#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


class EggNOGPipeline:
    def __init__(self, input_fasta: str, output_dir: str, cpus: int, itype: str):
        self.input_fasta = Path(input_fasta)
        self.output_dir = Path(output_dir)
        self.cpus = int(cpus)
        self.itype = itype
        self.prefix = self.input_fasta.stem

        data_dir = os.environ.get("EGGNOG_DATA_DIR")
        if not data_dir:
            raise RuntimeError("EGGNOG_DATA_DIR not set (should point to decompressed DB folder).")
        self.db_dir = Path(data_dir)

        self.required = [
            self.db_dir / "eggnog.db",
            self.db_dir / "eggnog.taxa.db",
            self.db_dir / "eggnog.taxa.db.traverse.pkl",
            self.db_dir / "eggnog_proteins.dmnd",
        ]

    @staticmethod
    def which(exe: str) -> Optional[str]:
        return shutil.which(exe)

    def run_cmd_stream(self, cmd: List[str]) -> None:
        print(f"[CMD] {' '.join(cmd)}", flush=True)
        proc = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr, text=True)
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)

    def ensure_db_present(self) -> None:
        missing = [p.name for p in self.required if not p.exists()]
        if missing:
            raise RuntimeError(f"DB missing in {self.db_dir}: {missing}")

    def run_emapper(self) -> Path:
        if not self.which("emapper.py"):
            raise RuntimeError("emapper.py not found in PATH. eggNOG-mapper installed?")

        self.output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            "emapper.py",
            "-i", str(self.input_fasta),
            "-o", self.prefix,
            "--output_dir", str(self.output_dir),
            "--data_dir", str(self.db_dir),
            "--itype", self.itype,
            "-m", "diamond",
            "--cpu", str(self.cpus),
            "--override",
        ]
        self.run_cmd_stream(cmd)

        out = self.output_dir / f"{self.prefix}.emapper.annotations"
        if not out.exists():
            raise FileNotFoundError(f"Expected emapper output not found: {out}")
        return out

    def clean_output(self, annotation_file: Path) -> Path:
        lines = annotation_file.read_text().splitlines()
        header = next((l for l in lines if l.startswith("#") and "query" in l.lower()), None)
        if header is None:
            raise ValueError("Missing header line containing 'query' in emapper output.")

        header = header.lstrip("#").strip()
        data = [l for l in lines if l and not l.startswith("#")]

        out_file = self.output_dir / "annotation.tsv"
        out_file.write_text("\n".join([header] + data) + "\n")
        return out_file

    def run_pipeline(self) -> Path:
        if not self.input_fasta.exists():
            raise FileNotFoundError(f"Input FASTA not found: {self.input_fasta}")
        self.ensure_db_present()
        ann = self.run_emapper()
        cleaned = self.clean_output(ann)
        print(f"[INFO] Done: {cleaned}", flush=True)
        return cleaned


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("input_fasta")
    p.add_argument("output_dir")
    p.add_argument("cpus", type=int)
    p.add_argument("itype", choices=["proteins", "CDS", "genome", "metagenome"])
    args = p.parse_args()

    EggNOGPipeline(args.input_fasta, args.output_dir, args.cpus, args.itype).run_pipeline()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
