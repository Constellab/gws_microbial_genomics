#!/usr/bin/env python3
from pathlib import Path
from typing import Final

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    Folder,
    InputSpec,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    ShellProxy,
    StrParam,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from .bactopiaMamba_env import BactopiaMambaShellProxyHelper


def _as_bool(v: str) -> bool:
    return str(v or "").strip().lower() in {"1", "true", "yes", "y", "on"}


@task_decorator(
    "BactopiaMetadata",
    human_name="Bactopia_metadata",
    short_description="Generate portable Bactopia samplesheet (filenames only) using `bactopia prepare`",
)
class BactopiaMetadata(Task):
    input_specs: Final[InputSpecs] = InputSpecs(
        {"fastq_folder": InputSpec(Folder, human_name="FASTQ folder")}
    )

    output_specs: Final[OutputSpecs] = OutputSpecs(
        {"samplesheet": OutputSpec(File, human_name="Bactopia samplesheet")}
    )

    config_specs: Final[ConfigSpecs] = ConfigSpecs(
        {
            "species": StrParam(default_value="", short_description="Species value (optional)"),
            "genome_size": StrParam(
                default_value="0", short_description="Genome size bp (optional)"
            ),
            "fastq_separator": StrParam(
                default_value="_", short_description="Separator used before mate token"
            ),
            "pe1_pattern": StrParam(default_value="", short_description="Pattern for mate 1 token (1,R1,r1)"),
            "pe2_pattern": StrParam(default_value="", short_description="Pattern for mate 2 token (2,r2,R2)"),
            "ont": StrParam(
                default_value="false",
                allowed_values=["false", "true"],
                short_description="Treat single-end reads as Oxford Nanopore reads",
            ),
        }
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        fastq_folder: Folder = inputs["fastq_folder"]
        shell: ShellProxy = BactopiaMambaShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "bactopia_prepare"
        out_dir.mkdir(parents=True, exist_ok=True)

        samplesheet_path = out_dir / "samples.tsv"
        wrapper_stderr = out_dir / "bactopia_prepare.wrapper.stderr.log"
        worker = Path(__file__).with_name("_bactopiaMetadata.py")

        species = (params["species"] or "").strip()
        genome_size = (params["genome_size"] or "0").strip() or "0"
        fastq_separator = (params["fastq_separator"] or "_").strip() or "_"
        pe1_pattern = (params["pe1_pattern"] or "([Aa]|[Rr]1|1)").strip() or "([Aa]|[Rr]1|1)"
        pe2_pattern = (params["pe2_pattern"] or "([Bb]|[Rr]2|2)").strip() or "([Bb]|[Rr]2|2)"
        ont = _as_bool(params["ont"])

        cmd_parts = [
            "python3",
            str(worker),
            "--fastqs",
            f'"{fastq_folder.path}"',
            "--out",
            f'"{samplesheet_path}"',
            "--fastq-separator",
            f'"{fastq_separator}"',
            "--pe1-pattern",
            f'"{pe1_pattern}"',
            "--pe2-pattern",
            f'"{pe2_pattern}"',
        ]
        if species:
            cmd_parts += ["--species", f'"{species}"']
        if genome_size and genome_size != "0":
            cmd_parts += ["--genome-size", f'"{genome_size}"']
        if ont:
            cmd_parts += ["--ont"]

        cmd_str = " ".join(cmd_parts) + f' 2> "{wrapper_stderr}"'
        print("[DEBUG] BactopiaPrepare cmd:", cmd_str)

        rc = shell.run(cmd_str, shell_mode=True)
        if rc != 0:
            raise RuntimeError(f"BactopiaPrepare failed (rc={rc}). See: {wrapper_stderr}")

        if not samplesheet_path.is_file() or samplesheet_path.stat().st_size == 0:
            raise RuntimeError(f"Samplesheet missing/empty: {samplesheet_path}")

        return {"samplesheet": File(str(samplesheet_path))}
