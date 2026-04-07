#!/usr/bin/env python3
# bactopia_explore.py
#
# Curate Bactopia outputs into flat ResourceSets (no nested ResourceSet).
# - Decompress selected useful compressed outputs (.tsv.gz/.csv.gz/.txt.gz, NanoPlot.tar.gz)
# - Keep only relevant artefacts per module (QC html, key TSV tables, sketcher reports, etc.)
# - Import tabular artefacts as Table when possible.

import gzip
import shutil
import tarfile
from pathlib import Path
from typing import Final, Optional, Tuple

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    Folder,
    InputSpec,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    ResourceSet,
    Settings,
    TableImporter,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

# -----------------------------
# Skip noisy stuff
# -----------------------------
SKIP_DIRS = {"logs"}
SKIP_PREFIXES = ("nf-",)
SKIP_SUFFIXES = (".err", ".log", ".trace", ".begin", ".out", ".run", ".sh")
SKIP_FILENAMES = {"versions.yml", "software_versions.yml", "software_versions_mqc.yml"}

# -----------------------------
# Decompression policy (targeted)
# -----------------------------
# We only decompress "useful" compressed outputs.
DECOMPRESS_GZ_SUFFIXES = (".tsv.gz", ".csv.gz", ".txt.gz")
DECOMPRESS_TAR_GZ_NAMES = ("nanoplot.tar.gz", "nanoplot-results.tar.gz")  # be permissive


def _is_skipped_path(p: Path) -> bool:
    if p.name in SKIP_FILENAMES:
        return True
    if any(p.name.startswith(x) for x in SKIP_PREFIXES):
        return True
    if any(p.name.endswith(x) for x in SKIP_SUFFIXES):
        return True
    if any(part in SKIP_DIRS for part in p.parts):
        return True
    return False


def _path_contains(rel: Path, token: str) -> bool:
    s = "/".join(rel.parts).lower()
    return f"/{token.lower()}/" in f"/{s}/"


def _infer_sample_id(rel: Path) -> Optional[str]:
    if not rel.parts:
        return None
    top = rel.parts[0]
    if top in {"bactopia-runs", "merged-results", "nf-reports"}:
        return None
    return top


def _infer_run_id(rel: Path) -> Optional[str]:
    if len(rel.parts) >= 2 and rel.parts[0] == "bactopia-runs":
        return rel.parts[1]
    return None


def _infer_annotator_tool(rel: Path) -> Optional[str]:
    s = "/".join(rel.parts).lower()
    if "/annotator/prokka/" in f"/{s}/":
        return "prokka"
    if "/annotator/bakta/" in f"/{s}/":
        return "bakta"
    return None


def _classify(rel: Path) -> str:
    # run-level
    if rel.parts and rel.parts[0] == "bactopia-runs":
        return "bactopia_runs"

    # sample-level modules
    if _path_contains(rel, "qc") or _path_contains(rel, "summary"):
        return "qc"
    if _path_contains(rel, "assembler"):
        return "assembler"
    if _path_contains(rel, "annotator"):
        return "annotator"
    if _path_contains(rel, "sketcher"):
        return "sketcher"
    if _path_contains(rel, "gather"):
        return "gather"
    if _path_contains(rel, "tools/mlst"):
        return "mlst"
    if _path_contains(rel, "tools/amrfinderplus"):
        return "amrfinderplus"

    # merged-results / nf-reports can appear outside bactopia-runs
    if _path_contains(rel, "merged-results") or _path_contains(rel, "nf-reports"):
        return "bactopia_runs"

    return "other"


def _keep(module: str, rel: Path) -> bool:
    name = rel.name.lower()

    if module == "qc":
        # Keep ONLY HTML reports:
        # - Illumina: fastp + fastqc (original + final)
        # - ONT: NanoPlot report HTML (original + final)
        if not name.endswith(".html"):
            return False
        return ("fastqc" in name) or ("fastp" in name) or ("nanoplot-report" in name)

    if module == "assembler":
        return name.endswith(".tsv") or name.endswith(".tsv.gz")

    if module == "annotator":
        return name.endswith(".tsv") or name.endswith(".tsv.gz")

    if module == "sketcher":
        # Your two key reports (often .txt, sometimes .txt.gz)
        return (
            name.endswith("-mash-refseq88-k21.txt")
            or name.endswith("-mash-refseq88-k21.txt.gz")
            or name.endswith("-sourmash-gtdb-rs207-k31.txt")
            or name.endswith("-sourmash-gtdb-rs207-k31.txt.gz")
        )

    if module in {"mlst", "amrfinderplus", "gather"}:
        return name.endswith(".tsv") or name.endswith(".tsv.gz")

    if module == "bactopia_runs":
        # keep: tables + html/dot
        if name.endswith(".html") or name.endswith(".dot"):
            return True
        if name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".txt") or name.endswith(".txt.gz") or name.endswith(".csv") or name.endswith(".csv.gz"):
            return True
        return False

    return False


def _label(module: str, rel: Path) -> str:
    sample = _infer_sample_id(rel)
    run_id = _infer_run_id(rel)

    if module == "bactopia_runs":
        prefix = f"RUN {run_id}" if run_id else "RUN"
        return f"{prefix} | {rel}"

    sid = sample or "sample"

    if module == "annotator":
        tool = _infer_annotator_tool(rel)
        if tool:
            return f"{sid} | Annotator({tool}) | {rel}"
        return f"{sid} | Annotator | {rel}"

    pretty = {
        "qc": "QC",
        "assembler": "Assembler",
        "annotator": "Annotator",
        "sketcher": "Sketcher",
        "mlst": "Sequence Typing (mlst)",
        "amrfinderplus": "Antimicrobial Resistance (amrfinderplus)",
        "gather": "Gather",
    }.get(module, module)

    return f"{sid} | {pretty} | {rel}"


def _gunzip_to(dst: Path, src: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as f_in, open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def _extract_nanoplot_tar_gz(src: Path, dst_dir: Path) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(src, "r:gz") as tf:
        tf.extractall(dst_dir)


def _materialize(src: Path, outdir: Path, work: Path) -> Tuple[Path, Path]:
    """
    Copy/decompress/extract into workdir.
    Returns (material_path, rel_for_label)
    """
    rel = src.relative_to(outdir)
    name = src.name.lower()
    suffixes = "".join(src.suffixes).lower()

    # Default: copy file as-is
    dst = work / rel
    dst.parent.mkdir(parents=True, exist_ok=True)

    # Decompress small useful gz tables/text
    if suffixes in DECOMPRESS_GZ_SUFFIXES:
        # remove only the last ".gz"
        dst = (work / rel).with_suffix("")  # drops ".gz"
        _gunzip_to(dst, src)
        return dst, dst.relative_to(work)

    # Extract NanoPlot tarball (ONT QC)
    # Bactopia commonly produces: <SAMPLE>-{final|original}_NanoPlot.tar.gz
    if name.endswith(".tar.gz") and ("nanoplot" in name):
        # Extract under a folder next to the tarball
        extract_dir = (work / rel).with_suffix("")  # drop .gz
        extract_dir = extract_dir.with_suffix("")   # drop .tar
        _extract_nanoplot_tar_gz(src, extract_dir)

        # Prefer the NanoPlot HTML report if present
        # We'll expose the report file itself (and keep rules will include it in QC anyway)
        candidates = list(extract_dir.rglob("*NanoPlot-report.html"))
        if candidates:
            report = candidates[0]
            return report, report.relative_to(work)

        # If no report, keep the extracted folder content isn't a Resource type here;
        # fall back to copying the tar.gz itself.
        shutil.copy2(src, dst)
        return dst, dst.relative_to(work)

    # Otherwise simple copy
    shutil.copy2(src, dst)
    return dst, dst.relative_to(work)


def _import_tsv(path: Path):
    # Force real tab separator
    return TableImporter.call(
        File(str(path)),
        {"delimiter": "tab", "header": 0, "file_format": "tsv"},
    )


def _looks_like_single_record_no_header_tsv(path: Path) -> bool:
    """
    Some TSV outputs (notably MLST) can contain a single result line without header.
    In that case, importing with header=0 produces a broken / empty table.
    """
    try:
        lines = [line for line in path.read_text(errors="ignore").splitlines() if line.strip()]
    except Exception:
        return False

    if len(lines) != 1:
        return False

    # One non-empty line with at least one tab strongly suggests a single TSV record
    # rather than a proper headered table.
    return "\t" in lines[0]


def _maybe_import_as_table(module: str, material: Path) -> object:
    name = material.name.lower()

    # MLST files are sometimes one-line TSV outputs without header.
    # Keep them as File in that case to avoid a bad table import.
    if module == "mlst" and name.endswith(".tsv"):
        if _looks_like_single_record_no_header_tsv(material):
            return File(str(material))
        try:
            return _import_tsv(material)
        except Exception:
            return File(str(material))

    if name.endswith(".tsv"):
        try:
            return _import_tsv(material)
        except Exception:
            return File(str(material))

    # sketcher txt are TSV-like (mash/sourmash reports)
    if module == "sketcher" and name.endswith(".txt"):
        try:
            return _import_tsv(material)
        except Exception:
            return File(str(material))

    # run-level txt/csv might be tabular; try TSV import
    if module == "bactopia_runs" and (name.endswith(".txt") or name.endswith(".csv")):
        try:
            return _import_tsv(material)
        except Exception:
            return File(str(material))

    return File(str(material))


@task_decorator(
    "BactopiaExplore",
    human_name="Bactopia_explore",
    short_description="Curated Bactopia outputs into flat ResourceSets (QC/Assembler/Annotator/...).",
)
class BactopiaExplore(Task):
    input_specs: Final[InputSpecs] = InputSpecs({
        "bactopia_results": InputSpec(Folder, human_name="Bactopia output folder"),
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "gather": OutputSpec(ResourceSet, human_name="Gather"),
        "qc": OutputSpec(ResourceSet, human_name="QC"),
        "assembler": OutputSpec(ResourceSet, human_name="Assembler"),
        "annotator": OutputSpec(ResourceSet, human_name="Annotator"),
        "sketcher": OutputSpec(ResourceSet, human_name="Sketcher"),
        "mlst": OutputSpec(ResourceSet, human_name="Sequence Typing (mlst)"),
        "amrfinderplus": OutputSpec(ResourceSet, human_name="Antimicrobial Resistance (amrfinderplus)"),
        "bactopia_runs": OutputSpec(ResourceSet, human_name="bactopia_runs"),
    })

    # No config params: this task is meant to be deterministic & curated.
    config_specs: Final[ConfigSpecs] = ConfigSpecs({})

    def run(self, p: ConfigParams, ins: TaskInputs) -> TaskOutputs:
        outdir = Path(ins["bactopia_results"].path)

        work = Path(Settings.make_temp_dir()) / "bactopia_explore"
        work.mkdir(parents=True, exist_ok=True)

        qc = ResourceSet()
        assembler = ResourceSet()
        annotator = ResourceSet()
        sketcher = ResourceSet()
        mlst = ResourceSet()
        amrfinderplus = ResourceSet()
        gather = ResourceSet()
        bactopia_runs = ResourceSet()

        for src in outdir.rglob("*"):
            if src.is_dir():
                continue
            if _is_skipped_path(src):
                continue

            rel = src.relative_to(outdir)
            module = _classify(rel)

            if module == "other":
                continue
            if not _keep(module, rel):
                continue

            material, rel_for_label = _materialize(src, outdir, work)
            label = _label(module, rel_for_label)
            res = _maybe_import_as_table(module, material)

            if module == "qc":
                qc.add_resource(res, label)
            elif module == "assembler":
                assembler.add_resource(res, label)
            elif module == "annotator":
                annotator.add_resource(res, label)
            elif module == "sketcher":
                sketcher.add_resource(res, label)
            elif module == "mlst":
                mlst.add_resource(res, label)
            elif module == "amrfinderplus":
                amrfinderplus.add_resource(res, label)
            elif module == "gather":
                gather.add_resource(res, label)
            elif module == "bactopia_runs":
                bactopia_runs.add_resource(res, label)

        return {
            "gather": gather,
            "qc": qc,
            "assembler": assembler,
            "annotator": annotator,
            "sketcher": sketcher,
            "mlst": mlst,
            "amrfinderplus": amrfinderplus,
            "bactopia_runs": bactopia_runs,
        }