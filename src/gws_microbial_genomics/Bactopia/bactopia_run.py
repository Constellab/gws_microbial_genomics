#!/usr/bin/env python3
import shlex
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


@task_decorator(
    "BactopiaRun",
    human_name="Bactopia_run",
    short_description="Run Bactopia from a portable samplesheet (basenames) and a FASTQ folder",
)
class BactopiaRun(Task):
    """
    **Bactopia** is a **Nextflow-based** pipeline for end-to-end analysis of **bacterial isolates** from sequencing reads.
    It typically performs:
    - Read **QC** (trimming, stats, reports)
    - **Assembly** (short-read, long-read, hybrid)
    - **Annotation** (e.g., Prokka)
    - Typing / characterization tools (e.g., **MLST**, AMRfinderPlus, sketching/Mash, etc.)
    - Report aggregation (e.g., MultiQC)

    > Bactopia is designed for **one isolate per sample**.
    > It is **not** a metagenomics pipeline for complex multi-species mixtures.

    ---

    ## Inputs

    This task runs Bactopia from:
    - a **portable samplesheet**
    - a **FASTQ folder**

    The samplesheet format documented below is the **wrapper-specific format expected by this task**.
    It is converted internally into a Bactopia-compatible `--samples` file.

    ### Supported sequencing inputs
    - **Illumina Paired-End (PE)**: `r1` + `r2`
    - **Illumina Single-End (SE)**: `r1` only
    - **Oxford Nanopore (ONT)**: `r1` only
    - **Hybrid Illumina + ONT**: `r1` + `r2` + `extra`

    ---

    ## Samplesheet format

    ### Required columns
    - `sample`: sample / isolate ID
    - `runtype`: one of `paired-end`, `single-end`, `ont`, `hybrid`, `short_polish`
    - `r1`: read 1 FASTQ filename (required for all supported runtype values)

    ### Conditionally required columns
    - `r2`: required for `paired-end`, `hybrid`, `short_polish`
    - `extra`: required for `hybrid`, `short_polish` and contains ONT long reads

    ---

    ## `runtype` values and meaning

    `runtype` selects how the row is interpreted.

    - `paired-end`
    Illumina paired-end input. Uses `r1` and `r2`.

    - `single-end`
    Illumina single-end input. Uses `r1` only. `r2` must be empty.

    - `ont`
    ONT-only input. Uses `r1` only. `r2` must be empty.

    - `hybrid`
    Illumina paired-end + ONT input. Uses `r1` + `r2` for Illumina reads and `extra` for ONT reads.

    - `short_polish`
    Illumina paired-end + ONT input. Uses `r1` + `r2` for Illumina reads and `extra` for ONT reads, with short-read polishing workflow.

    > Hybrid choice is explicit.
    > The user must set `runtype=hybrid` or `runtype=short_polish`.

    ---

    ## Column semantics

    ### `r1`
    - Always required
    - `paired-end`: Illumina R1 FASTQ
    - `single-end`: Illumina single-end FASTQ
    - `ont`: ONT FASTQ
    - `hybrid` / `short_polish`: Illumina R1 FASTQ

    ### `r2`
    - Required only for `paired-end`, `hybrid`, `short_polish`
    - Must be empty for `single-end` and `ont`

    ### `extra`
    - Required only for `hybrid` and `short_polish`
    - Contains the ONT long reads associated with the same sample

    ---

    ## Recommended usage: `hybrid` vs `short_polish`

    ### Use `hybrid` when:
    - ONT data is noisy and/or low coverage
    - You want a classic short-read-first hybrid approach

    ### Use `short_polish` when:
    - ONT data is modern and high coverage
    - You want long-read-first assembly with short-read polishing

    ---

    ## Example samplesheet

    sample	runtype	r1	r2	extra
    sample1	paired-end	sample1_R1.fastq.gz	sample1_R2.fastq.gz
    sample2	single-end	sample2.fastq.gz
    sample3	ont	sample3.fastq.gz
    EcoliHybrid	hybrid	Ecoli_R1.fastq.gz	Ecoli_R2.fastq.gz	Ecoli_ont.fastq.gz
    EcoliPolish	short_polish	EcoliPolish_R1.fastq.gz	EcoliPolish_R2.fastq.gz	EcoliPolish_ont.fastq.gz
    """

    input_specs: Final[InputSpecs] = InputSpecs(
        {
            "samplesheet": InputSpec(File, human_name="Portable samplesheet (TSV)"),
            "fastq_folder": InputSpec(Folder, human_name="FASTQ folder"),
        }
    )

    output_specs: Final[OutputSpecs] = OutputSpecs(
        {"bactopia_results": OutputSpec(Folder, human_name="Bactopia output folder")}
    )

    config_specs: Final[ConfigSpecs] = ConfigSpecs(
        {
            "max_cpus": StrParam(default_value="", short_description="bactopia --max_cpus (optional)"),
            "coverage": StrParam(default_value="", short_description="bactopia --coverage (optional)"),
        }
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        samplesheet: File = inputs["samplesheet"]
        fastq_folder: Folder = inputs["fastq_folder"]

        shell: ShellProxy = BactopiaMambaShellProxyHelper.create_proxy(self.message_dispatcher)

        work_dir = Path(shell.working_dir)
        out_dir = work_dir / "bactopia_run"
        out_dir.mkdir(parents=True, exist_ok=True)

        bactopia_outdir = out_dir / "bactopia_out"
        worker = Path(__file__).with_name("_bactopia_run.py")

        max_cpus = (params.get("max_cpus", "") or "").strip()
        coverage = (params.get("coverage", "") or "").strip()

        args = [
            "python3",
            str(worker),
            "--samplesheet", str(samplesheet.path),
            "--fastq-folder", str(fastq_folder.path),
            "--outdir", str(bactopia_outdir),
        ]
        if max_cpus:
            args += ["--max-cpus", max_cpus]
        if coverage:
            args += ["--coverage", coverage]

        cmd_str = " ".join(shlex.quote(a) for a in args)
        print("[DEBUG] BactopiaRun cmd:", cmd_str)

        rc = shell.run(cmd_str, shell_mode=True)
        if rc != 0:
            raise RuntimeError("BactopiaRun failed (rc=1).")

        if not bactopia_outdir.is_dir():
            raise RuntimeError(f"Bactopia output folder not found: {bactopia_outdir}")

        return {"bactopia_results": Folder(str(bactopia_outdir))}
