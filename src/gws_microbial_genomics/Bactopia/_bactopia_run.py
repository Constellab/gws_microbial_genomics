#!/usr/bin/env python3
import argparse
import os
import shutil
import stat
import subprocess
from pathlib import Path
from typing import Dict, List


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--fastq-folder", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--max-cpus", default="")
    p.add_argument("--coverage", default="")
    return p.parse_args()


def _split_list(cell: str) -> List[str]:
    v = (cell or "").strip()
    if not v:
        return []
    return [s.strip() for s in v.split(",") if s.strip()]


def _resolve_cell(cell: str, fastq_folder: Path) -> str:
    parts = _split_list(cell)
    # portable samplesheet (basenames) -> abs paths under fastq_folder
    return ",".join(str((fastq_folder / Path(p).name).resolve()) for p in parts)


def _collect_missing_files(cell: str, fastq_folder: Path) -> List[Path]:
    missing: List[Path] = []
    for p in _split_list(cell):
        fp = fastq_folder / Path(p).name
        if not fp.is_file():
            missing.append(fp)
    return missing


def _require(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(msg)


def _write_executable(path: Path, content: str) -> None:
    path.write_text(content)
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def make_runtime(outdir: Path) -> Dict[str, Path]:
    base = outdir.parent / ".bactopia_runtime"
    d: Dict[str, Path] = {
        "base": base,
        "home": base / "home",
        "tmp": base / "tmp",
        "nxf_home": base / "nxf_home",
        "work": base / "work",
        "cache": base / "bactopia_cache",
        "conda": base / "bactopia_cache" / "conda",
        "conda_pkgs": base / "bactopia_cache" / "conda" / "pkgs",
        "datasets": base / "bactopia_cache" / "datasets",
        "singularity": base / "bactopia_cache" / "singularity",
        "bin": base / "bin",
    }
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    return d


def install_conda_shim(runtime: Dict[str, Path], real_conda: str) -> None:
    """
    Some conda builds crash on: conda create --mkdir ...
    Shim removes only the '--mkdir' argument.
    (Python shim, no embedded bash.)
    """
    shim_py = """#!/usr/bin/env python3
import os, sys
real = os.environ.get("REAL_CONDA", "conda")
args = [a for a in sys.argv[1:] if a != "--mkdir"]
os.execv(real, [real, *args])
"""
    _write_executable(runtime["bin"] / "conda", shim_py)


def write_nfconfig_ignore_dumpsoftwareversions(path: Path) -> None:
    """
    Your environment sometimes fails in DUMPSOFTWAREVERSIONS (versions.yml YAML parse).
    We do NOT disable the process; we ignore its error so the run continues.
    """
    path.write_text(
        "\n".join(
            [
                "process {",
                "  withName: /.*DUMPSOFTWAREVERSIONS.*/ {",
                "    errorStrategy = 'ignore'",
                "    maxRetries    = 0",
                "  }",
                "}",
                "",
            ]
        )
    )


# Official-ish runtype values used by Bactopia prepare/docs in samplesheet mode
VALID_RUNTYPES = {
    "paired-end",
    "single-end",
    "ont",
    "hybrid",
    "short_polish",
    "merge-pe",
    "hybrid-merge-pe",
    "short_polish-merge-pe",
    "merge-se",
    "assembly",
}

# Accept user/legacy variants, normalize to prepare-style naming
ALIASES = {
    "paired": "paired-end",
    "pe": "paired-end",
    "paired_end": "paired-end",
    "paired end": "paired-end",
    "single": "single-end",
    "se": "single-end",
    "single_end": "single-end",
    "single end": "single-end",
    "short-polish": "short_polish",
    "short polish": "short_polish",
}


def normalize_runtype(rt: str) -> str:
    raw = (rt or "").strip().lower()
    _require(raw != "", "Samplesheet: runtype vide")
    norm = ALIASES.get(raw, raw)
    _require(
        norm in VALID_RUNTYPES,
        f"Samplesheet: runtype invalide '{rt}'. Valeurs supportées: {', '.join(sorted(VALID_RUNTYPES))}",
    )
    # IMPORTANT: we DO NOT accept ambiguous "paired-end and ont" in runtype.
    if "paired" in norm and "ont" in norm:
        raise RuntimeError(
            "Samplesheet: runtype ambigu (contient à la fois paired et ont). "
            "Pour un hybride, mets explicitement runtype=hybrid OU runtype=short_polish "
            "et place les reads ONT dans la colonne extra."
        )
    return norm


def resolve_samplesheet(in_tsv: Path, out_tsv: Path, fastq_folder: Path) -> None:
    lines = in_tsv.read_text().splitlines()
    _require(bool(lines), "Samplesheet vide")

    header = lines[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}

    _require("sample" in idx, "Samplesheet: colonne requise manquante: sample")
    _require("runtype" in idx, "Samplesheet: colonne requise manquante: runtype")
    _require("r1" in idx, "Samplesheet: colonne requise manquante: r1")

    has_r2 = "r2" in idx
    has_extra = "extra" in idx

    out_lines = [lines[0]]
    missing_all: List[Path] = []

    for li, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue

        cols = line.split("\t")
        if len(cols) < len(header):
            cols += [""] * (len(header) - len(cols))

        sample = (cols[idx["sample"]] or "").strip()
        _require(sample != "", f"Samplesheet: sample vide (ligne {li})")

        rt_raw = cols[idx["runtype"]]
        rt = normalize_runtype(rt_raw)
        cols[idx["runtype"]] = rt

        r1 = (cols[idx["r1"]] or "").strip()
        r2 = (cols[idx["r2"]] or "").strip() if has_r2 else ""
        extra = (cols[idx["extra"]] or "").strip() if has_extra else ""

        # ---- Validate by runtype (Bactopia rules in samplesheet mode) ----
        if rt in {"paired-end", "merge-pe"}:
            _require(r1 != "", f"Samplesheet: r1 manquant (ligne {li})")
            _require(has_r2 and r2 != "", f"Samplesheet: r2 manquant pour runtype '{rt}' (ligne {li})")
            # extra usually empty; we don't hard-fail, but it's suspicious
        elif rt in {"single-end", "merge-se"}:
            _require(r1 != "", f"Samplesheet: r1 manquant (ligne {li})")
            _require(r2 == "", f"Samplesheet: r2 doit être vide pour runtype '{rt}' (ligne {li})")
        elif rt == "ont":
            # ONT-only: reads are in r1 (per prepare example)
            _require(r1 != "", f"Samplesheet: r1 (ONT) manquant (ligne {li})")
            _require(r2 == "", f"Samplesheet: r2 doit être vide pour runtype 'ont' (ligne {li})")
            # extra should be empty for ont-only
        elif rt in {"hybrid", "short_polish", "hybrid-merge-pe", "short_polish-merge-pe"}:
            _require(r1 != "", f"Samplesheet: r1 manquant (ligne {li})")
            _require(has_r2 and r2 != "", f"Samplesheet: r2 manquant pour runtype '{rt}' (ligne {li})")
            _require(has_extra and extra != "", f"Samplesheet: extra (ONT) manquant pour runtype '{rt}' (ligne {li})")
        elif rt == "assembly":
            _require(has_extra and extra != "", f"Samplesheet: extra (assembly) manquant pour runtype 'assembly' (ligne {li})")

        # ---- Check missing files ----
        missing_all += _collect_missing_files(r1, fastq_folder)
        if has_r2:
            missing_all += _collect_missing_files(r2, fastq_folder)
        if has_extra:
            missing_all += _collect_missing_files(extra, fastq_folder)

        # ---- Resolve paths ----
        if r1:
            cols[idx["r1"]] = _resolve_cell(r1, fastq_folder)
        if has_r2 and r2:
            cols[idx["r2"]] = _resolve_cell(r2, fastq_folder)
        if has_extra and extra:
            cols[idx["extra"]] = _resolve_cell(extra, fastq_folder)

        out_lines.append("\t".join(cols[: len(header)]))

    if missing_all:
        uniq = sorted({str(p) for p in missing_all})
        preview = "\n".join(uniq[:30])
        more = "" if len(uniq) <= 30 else f"\n... +{len(uniq) - 30} autres"
        raise RuntimeError(f"FASTQ(s) introuvable(s) dans fastq-folder:\n{preview}{more}")

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_tsv.write_text("\n".join(out_lines) + "\n")


def run_bactopia(cmd: List[str], env: Dict[str, str]) -> None:
    res = subprocess.run(cmd, text=True, capture_output=True, env=env)
    if res.returncode != 0:
        out = ((res.stdout or "") + "\n" + (res.stderr or "")).strip()
        out = "\n".join(out.splitlines()[-300:]) if out else "(aucune sortie)"
        raise RuntimeError(f"Bactopia a échoué (code {res.returncode}).\n{out}")


def main():
    a = parse_args()

    samplesheet = Path(a.samplesheet)
    fastq_folder = Path(a.fastq_folder)
    outdir = Path(a.outdir)

    _require(samplesheet.is_file(), f"Samplesheet introuvable: {samplesheet}")
    _require(fastq_folder.is_dir(), f"FASTQ folder introuvable: {fastq_folder}")
    outdir.mkdir(parents=True, exist_ok=True)

    real_conda = shutil.which("conda")
    _require(bool(real_conda), "conda introuvable dans l'environnement")

    runtime = make_runtime(outdir)
    install_conda_shim(runtime, real_conda)

    env = dict(os.environ)
    env["REAL_CONDA"] = real_conda
    env["PATH"] = str(runtime["bin"]) + os.pathsep + env.get("PATH", "")
    env["HOME"] = str(runtime["home"])
    env["TMPDIR"] = str(runtime["tmp"])
    env["NXF_HOME"] = str(runtime["nxf_home"])
    env["NXF_ANSI_LOG"] = "false"
    env["CONDA_ALWAYS_YES"] = "true"
    env["CONDA_PKGS_DIRS"] = str(runtime["conda_pkgs"])
    env["NXF_CONDA_CACHEDIR"] = str(runtime["conda"])
    env["BACTOPIA_CACHEDIR"] = str(runtime["cache"])

    nfconfig = runtime["base"] / "local_overrides.config"
    write_nfconfig_ignore_dumpsoftwareversions(nfconfig)

    resolved = runtime["base"] / "samples.resolved.tsv"
    resolve_samplesheet(samplesheet, resolved, fastq_folder)

    cmd = [
        "bactopia",
        "--samples", str(resolved),
        "--outdir", str(outdir),
        "--condadir", str(runtime["conda"]),
        "--datasets_cache", str(runtime["datasets"]),
        "--singularity_cache", str(runtime["singularity"]),
        "--nfconfig", str(nfconfig),
    ]

    if str(a.max_cpus).strip():
        cmd += ["--max_cpus", str(a.max_cpus).strip()]
    if str(a.coverage).strip():
        cmd += ["--coverage", str(a.coverage).strip()]

    # NOTE: no --ont / no --short_polish / no --hybrid here.
    # Those are driven by per-row 'runtype' in the samplesheet.
    run_bactopia(cmd, env)


if __name__ == "__main__":
    main()
