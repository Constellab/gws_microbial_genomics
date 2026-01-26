#!/usr/bin/env python3
import os
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
    IntParam,
    OutputSpec,
    OutputSpecs,
    ShellProxy,
    StrParam,
    Table,
    TableImporter,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from .eggnog_env_task import EggnogShellProxyHelper


@task_decorator(
    "Emapper",
    human_name="eggNOG Mapper",
    short_description="eggNOG-mapper using pre-downloaded DB folder (no download, no extract).",
)
class EggnogMapperTask(Task):

    input_specs: Final[InputSpecs] = InputSpecs({
        "fasta": InputSpec(File, human_name="FASTA file",
                           short_description="Protein or transcript FASTA to annotate"),
        "db_dir": InputSpec(Folder, human_name="eggNOG DB folder",
                            short_description="Output of EggnogDBDownloadTask (decompressed DB dir)"),
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "annotation_table": OutputSpec(Table, human_name="eggNOG annotations",
                                       short_description="Cleaned TSV output from eggNOG-mapper"),
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "itype": StrParam(default_value="proteins",
                          allowed_values=["proteins", "CDS", "genome", "metagenome"],
                          short_description="Sequence type for eggNOG input"),
        "cpus": IntParam(default_value=25, min_value=1,
                         short_description="Number of CPU threads to use"),
    })

    python_file_path: Final[str] = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), "_eggnog_emapper.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        fasta: File = inputs["fasta"]
        db_dir: Folder = inputs["db_dir"]

        itype: str = params["itype"]
        cpus: int = int(params["cpus"])

        shell: ShellProxy = EggnogShellProxyHelper.create_proxy(self.message_dispatcher)
        work_dir = Path(shell.working_dir)

        prefix = Path(fasta.path).stem
        output_dir = work_dir / f"{prefix}_output"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Validate DB folder content BEFORE running emapper
        db_path = Path(db_dir.path)
        required = [
            db_path / "eggnog.db",
            db_path / "eggnog.taxa.db",
            db_path / "eggnog.taxa.db.traverse.pkl",
            db_path / "eggnog_proteins.dmnd",
        ]
        missing = [p.name for p in required if not p.exists()]
        if missing:
            raise RuntimeError(f"DB folder invalid: missing {missing} in {db_path}")

        export_env = f"export EGGNOG_DATA_DIR={shlex.quote(str(db_path))};"

        cmd = " ".join([
            export_env,
            "python3",
            shlex.quote(self.python_file_path),
            shlex.quote(str(fasta.path)),
            shlex.quote(str(output_dir)),
            shlex.quote(str(cpus)),
            shlex.quote(str(itype)),
        ])

        ret_code = shell.run(cmd, shell_mode=True)
        if ret_code != 0:
            raise RuntimeError("eggNOG-mapper pipeline failed. Check logs for more info.")

        cleaned_output_path = output_dir / "annotation.tsv"
        if not cleaned_output_path.exists():
            raise FileNotFoundError(f"Expected cleaned annotation file not found: {cleaned_output_path}")

        annotation_table = TableImporter.call(
            File(str(cleaned_output_path)),
            {
                "delimiter": "tab",
                "header": 0,
                "file_format": "tsv",
                "index_column": None,
                "comment": None,
            }
        )

        return {"annotation_table": annotation_table}
