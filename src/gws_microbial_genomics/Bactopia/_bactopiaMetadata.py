#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Generate portable Bactopia samplesheet (basenames only)")
    p.add_argument("--fastqs", required=True, help="Folder containing compressed FASTQ files")
    p.add_argument("--out", required=True, help="Output samplesheet path (TSV)")

    # Optional metadata fillers
    p.add_argument("--species", default="", help="Species value (optional)")
    p.add_argument("--genome-size", default="0", help="Genome size bp (optional)")

    # Paired-end naming rules
    p.add_argument("--fastq-separator", default="_", help="Separator used before mate token (default: _)")
    p.add_argument("--pe1-pattern", default="([Aa]|[Rr]1|1)", help="Pattern for mate 1 token")
    p.add_argument("--pe2-pattern", default="([Bb]|[Rr]2|2)", help="Pattern for mate 2 token")

    # Nanopore mode
    p.add_argument("--ont", action="store_true", help="Treat single-end reads as ONT reads")

    return p.parse_args()


def _basename_or_list(value: str) -> str:
    """Convert '/a/b/x.fastq.gz,/c/d/y.fastq.gz' -> 'x.fastq.gz,y.fastq.gz' """
    v = (value or "").strip()
    if not v:
        return ""
    parts = [s.strip() for s in v.split(",") if s.strip()]
    return ",".join(Path(s).name for s in parts)


def run_prepare(cmd):
    return subprocess.run(cmd, text=True, capture_output=True)


def main():
    a = parse_args()
    fastq_dir = Path(a.fastqs)
    out_path = Path(a.out)

    if not fastq_dir.is_dir():
        raise RuntimeError(f"FASTQ folder does not exist: {fastq_dir}")

    # sanity: fichiers compressés attendus
    if not list(fastq_dir.glob("*.gz")):
        raise RuntimeError(f"No compressed FASTQ files (*.gz) found in: {fastq_dir}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Variant A: options avec tirets (souvent montré dans la doc)
    cmd_a = [
        "bactopia", "prepare",
        "--path", str(fastq_dir),
        "--fastq-separator", a.fastq_separator,
        "--pe1-pattern", a.pe1_pattern,
        "--pe2-pattern", a.pe2_pattern,
    ]
    # Variant B: options avec underscores (vu dans certaines versions / aides)
    cmd_b = [
        "bactopia", "prepare",
        "--path", str(fastq_dir),
        "--fastq_separator", a.fastq_separator,
        "--pe1_pattern", a.pe1_pattern,
        "--pe2_pattern", a.pe2_pattern,
    ]

    if a.ont:
        cmd_a.append("--ont")
        cmd_b.append("--ont")

    if a.species.strip():
        cmd_a += ["--species", a.species.strip()]
        cmd_b += ["--species", a.species.strip()]

    gs = str(a.genome_size).strip()
    if gs and gs != "0":
        cmd_a += ["--genome-size", gs]
        cmd_b += ["--genome-size", gs]

    res = run_prepare(cmd_a)
    if res.returncode != 0:
        # si erreur liée aux options, on tente l'autre variante
        err = (res.stderr or "")
        if ("unrecognized arguments" in err) or ("no such option" in err) or ("unknown option" in err):
            res = run_prepare(cmd_b)

    if res.returncode != 0:
        raise RuntimeError(f"bactopia prepare failed:\n{res.stderr}")

    if not res.stdout.strip():
        raise RuntimeError("bactopia prepare succeeded but produced empty output")

    # Parse TSV and replace r1/r2/extra with basenames
    lines = res.stdout.splitlines()
    header = lines[0].split("\t")
    col_idx = {name: i for i, name in enumerate(header)}

    for required in ("r1", "r2", "extra"):
        if required not in col_idx:
            raise RuntimeError(f"Unexpected samplesheet header, missing column: {required}")

    out_lines = [lines[0]]
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < len(header):
            cols += [""] * (len(header) - len(cols))

        cols[col_idx["r1"]] = _basename_or_list(cols[col_idx["r1"]])
        cols[col_idx["r2"]] = _basename_or_list(cols[col_idx["r2"]])
        cols[col_idx["extra"]] = _basename_or_list(cols[col_idx["extra"]])

        out_lines.append("\t".join(cols[:len(header)]))

    out_path.write_text("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    main()
