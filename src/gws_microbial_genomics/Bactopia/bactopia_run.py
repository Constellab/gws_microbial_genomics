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
    > It is **not** a metagenomics pipeline for complex multi-species mixtures (e.g., raw soil samples).

    ---

    ## Inputs

    ### 1) FASTQ reads
    Bactopia supports the following sequencing inputs:
    - **Illumina Paired-End (PE)**: `R1` + `R2`
    - **Illumina Single-End (SE)**: `R1` only
    - **Oxford Nanopore (ONT)**: long reads (usually one FASTQ)

    ### 2) Samplesheet (TSV)
    For multi-sample runs, Bactopia takes a **tab-separated samplesheet** describing one sample per row.

    ---

    ## Samplesheet format

    ### Required columns (minimum)
    - `sample`: sample / isolate ID
    - `runtype`: how Bactopia should interpret the row
    - `r1`: read 1 FASTQ path (or comma-separated list)
    - `r2`: read 2 FASTQ path (or comma-separated list, required for PE/hybrid modes)

    ### Optional column
    - `extra`: **either** long reads (ONT) **or** a precomputed assembly, depending on `runtype`

    ---

    ## `runtype` values and meaning

    `runtype` is the **key** field: it selects the input interpretation and assembly strategy.

    Common supported values:
    - `paired-end`
    Illumina Paired-End. Uses `r1` and `r2`.

    - `single-end`
    Illumina Single-End. Uses `r1` only.

    - `ont`
    ONT long reads only. In many wrappers, the ONT FASTQ is stored in `r1` and `r2` is empty.

    - `hybrid`
    Illumina PE + ONT. Uses `r1`+`r2` (Illumina) and `extra` (ONT).
    Hybrid assembly approach: **short reads first**, then long reads bridge gaps (commonly via Unicycler).

    - `short_polish`
    Illumina PE + ONT. Uses `r1`+`r2` (Illumina) and `extra` (ONT).
    Hybrid assembly with short-read polishing: **long reads first**, then polish with short reads (often better for modern ONT data).

    > Important: **Hybrid choice is explicit**.
    > The user selects `hybrid` or `short_polish` in the samplesheet. No auto-selection.

    ---

    ## Column semantics (`r1`, `r2`, `extra`)

    ### `r1`
    - Always required.
    - For `paired-end`: Illumina R1 FASTQ(s).
    - For `single-end`: Illumina SE FASTQ(s).
    - For `ont`: ONT FASTQ(s), if your wrapper uses `r1` as the ONT container.

    ### `r2`
    - Required only for `paired-end`, `hybrid`, `short_polish`.
    - Must be empty for `single-end` and typically empty for `ont`.

    ### `extra`
    Per Bactopia documentation:
    > **extra**: Either the assembly or long reads associated with a sample

    Practical usage:
    - `hybrid` / `short_polish`: `extra` contains the **ONT long reads** for that same sample.
    - `assembly`: `extra` may contain a pre-built assembly FASTA (if your workflow supports it).

    ---

    ## Recommended usage: `hybrid` vs `short_polish`

    ### Use `hybrid` when:
    - ONT data is noisy and/or low coverage
    - You want a classic short-read-first hybrid approach

    ### Use `short_polish` when:
    - ONT data is **modern** and **high coverage**
    - You want long-read-first assembly with short-read polishing (often faster/better)

    ---

    ## Example samplesheet (mixed run)

    ```tsv
    sample	runtype	r1	r2	extra
    sample1	paired-end	sample1_1.fastq.gz	sample1_2.fastq.gz
    sample2	paired-end	sample2_1.fastq.gz	sample2_2.fastq.gz
    sample3	ont	sample3.fastq.gz
    EcoliHybrid	short_polish	Ecoli_R1.fastq.gz	Ecoli_R2.fastq.gz	Ecoli_ont.fastq.gz

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
