#!/usr/bin/env python3
# bactopia_explore.py

import shutil
from pathlib import Path
from typing import Final, Optional

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
    StrParam,
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


def _copy_to_work(src: Path, outdir: Path, work: Path) -> tuple[Path, Path]:
    rel = src.relative_to(outdir)
    dst = work / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst, rel


def _import_tsv_table(path: Path):
    # Force TSV tab separator
    return TableImporter.call(
        File(str(path)),
        {"delimiter": "tab", "header": 0, "file_format": "tsv"},
    )


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


def _path_contains(rel: Path, token: str) -> bool:
    s = "/".join(rel.parts).lower()
    return f"/{token.lower()}/" in f"/{s}/"


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

    # merged-results / nf-reports often live under bactopia-runs anyway
    if _path_contains(rel, "merged-results") or _path_contains(rel, "nf-reports"):
        return "bactopia_runs"

    return "other"


# -----------------------------
# Per-module "keep" rules
# -----------------------------
def _keep(module: str, rel: Path) -> bool:
    name = rel.name.lower()
    suf = "".join(Path(name).suffixes)

    if module == "qc":
        # ONLY HTML: fastp html + fastqc html (original + final)
        return name.endswith(".html") and ("fastqc" in name or "fastp" in name)

    if module == "assembler":
        # keep assembler TSV summaries only (your request)
        return name.endswith(".tsv")

    if module == "annotator":
        # keep TSV only (your request) - you can add gff/gbk later if you want
        return name.endswith(".tsv")

    if module == "sketcher":
        # keep only the 2 key reports you mentioned
        return (
            name.endswith("-mash-refseq88-k21.txt")
            or name.endswith("-sourmash-gtdb-rs207-k31.txt")
        )

    if module == "mlst":
        return name.endswith(".tsv")

    if module == "amrfinderplus":
        return name.endswith(".tsv")

    if module == "gather":
        return name.endswith(".tsv")

    if module == "bactopia_runs":
        # everything to table except .html and .dot
        if name.endswith(".tsv"):
            return True
        if name.endswith(".html") or name.endswith(".dot"):
            return True
        return False

    return False


def _label(module: str, rel: Path) -> str:
    sample = _infer_sample_id(rel)
    run_id = _infer_run_id(rel)

    if module == "bactopia_runs":
        prefix = f"RUN {run_id}" if run_id else "RUN"
        return f"{prefix} | {rel}"
    else:
        sid = sample or "sample"
        return f"{sid} | {module} | {rel}"


@task_decorator(
    "BactopiaExplore",
    human_name="Bactopia_explore",
    short_description="Curated Bactopia outputs into flat ResourceSets (no logs, tables when possible).",
)
class BactopiaExplore(Task):

    input_specs: Final[InputSpecs] = InputSpecs({
        "bactopia_results": InputSpec(Folder, human_name="Bactopia output folder"),
    })

    # Ordered as you requested (logical appearance)
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "qc": OutputSpec(ResourceSet, human_name="QC"),
        "assembler": OutputSpec(ResourceSet, human_name="Assembler"),
        "annotator": OutputSpec(ResourceSet, human_name="Annotator"),
        "sketcher": OutputSpec(ResourceSet, human_name="Sketcher"),
        "mlst": OutputSpec(ResourceSet, human_name="Sequence Typing (mlst)"),
        "amrfinderplus": OutputSpec(ResourceSet, human_name="Antimicrobial Resistance (amrfinderplus)"),
        "gather": OutputSpec(ResourceSet, human_name="Gather"),
        "bactopia_runs": OutputSpec(ResourceSet, human_name="bactopia_runs"),
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "include_bactopia_runs": StrParam(
            default_value="true",
            allowed_values=["true", "false"],
            short_description="Include bactopia-runs reports + merged tables",
        ),
    })

    def run(self, p: ConfigParams, ins: TaskInputs) -> TaskOutputs:
        outdir = Path(ins["bactopia_results"].path)

        work = Path(Settings.make_temp_dir()) / "bactopia_explore"
        work.mkdir(parents=True, exist_ok=True)

        include_runs = str(p.get("include_bactopia_runs", "true")).lower() in {"true", "1", "yes", "y", "on"}

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
            if module == "bactopia_runs" and not include_runs:
                continue
            if not _keep(module, rel):
                continue

            material, rel_for_label = _copy_to_work(src, outdir, work)
            label = _label(module, rel_for_label)

            name = rel_for_label.name.lower()

            # Tables
            if module in {"assembler", "annotator", "mlst", "amrfinderplus", "gather"} and name.endswith(".tsv"):
                try:
                    res = _import_tsv_table(material)
                except Exception:
                    res = File(str(material))
            elif module == "sketcher" and name.endswith(".txt"):
                # these are TSV-like; import as TSV table
                try:
                    res = _import_tsv_table(material)
                except Exception:
                    res = File(str(material))
            elif module == "bactopia_runs" and name.endswith(".tsv"):
                try:
                    res = _import_tsv_table(material)
                except Exception:
                    res = File(str(material))
            else:
                # QC html, run html/dot
                res = File(str(material))

            # Dispatch
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
            "qc": qc,
            "assembler": assembler,
            "annotator": annotator,
            "sketcher": sketcher,
            "mlst": mlst,
            "amrfinderplus": amrfinderplus,
            "gather": gather,
            "bactopia_runs": bactopia_runs,
        }
