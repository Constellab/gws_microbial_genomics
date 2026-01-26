#!/usr/bin/env python3
import gzip
import os
import shutil
import tarfile
import time
import urllib.request
from pathlib import Path
from typing import Final, Optional

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    Folder,
    OutputSpec,
    OutputSpecs,
    StrParam,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)


def _download_stream(url: str, dest: Path, timeout: int = 60) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f)


def _download_with_retries(url: str, dest: Path, retries: int = 5, sleep_sec: int = 5) -> None:
    last_err: Optional[Exception] = None
    for i in range(1, retries + 1):
        try:
            dest.parent.mkdir(parents=True, exist_ok=True)
            # overwrite partials to avoid corrupt states
            if dest.exists():
                dest.unlink(missing_ok=True)
            print(f"[INFO] Downloading ({i}/{retries}): {url}", flush=True)
            _download_stream(url, dest)
            return
        except Exception as e:
            last_err = e
            print(f"[WARN] Download failed ({i}/{retries}): {e}", flush=True)
            time.sleep(sleep_sec)
    raise RuntimeError(f"Failed to download after {retries} retries: {url}\nLast error: {last_err}")


def _assert_min_size(path: Path, min_bytes: int, label: str) -> None:
    size = path.stat().st_size if path.exists() else 0
    if size < min_bytes:
        raise RuntimeError(f"{label} too small ({size} bytes), expected >= {min_bytes}. File: {path}")


def _gunzip(gz_path: Path, out_path: Path) -> None:
    print(f"[INFO] gunzip: {gz_path} -> {out_path}", flush=True)
    with gzip.open(gz_path, "rb") as fin, open(out_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def _untar_gz(tar_gz_path: Path, out_dir: Path) -> None:
    print(f"[INFO] untar: {tar_gz_path} -> {out_dir}", flush=True)
    with tarfile.open(tar_gz_path, "r:gz") as tar:
        tar.extractall(path=out_dir)


@task_decorator(
    "EggnogDBDownload",
    human_name="eggNOG DB download",
    short_description="Download + extract eggNOG-mapper DBs once (Python), output as Folder.",
)
class EggnogDBDownloadTask(Task):

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "db_dir": OutputSpec(
            Folder,
            human_name="eggNOG DB folder",
            short_description="Decompressed eggNOG-mapper data directory",
        ),
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "base_url": StrParam(
            default_value="http://eggnog6.embl.de/download/emapperdb-5.0.2",
            short_description="Base URL hosting eggNOG-mapper DB files (HTTP recommended).",
        ),
    })

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        base_url: str = params["base_url"]

        # Folder artifact location
        work_dir = Path(self.working_dir) if hasattr(self, "working_dir") else Path(os.getcwd())
        db_dir = work_dir / "eggnog_data"
        db_dir.mkdir(parents=True, exist_ok=True)

        required = [
            db_dir / "eggnog.db",
            db_dir / "eggnog.taxa.db",
            db_dir / "eggnog.taxa.db.traverse.pkl",
            db_dir / "eggnog_proteins.dmnd",
        ]
        if all(p.exists() for p in required):
            print("[INFO] DB already present, skipping download.", flush=True)
            return {"db_dir": Folder(str(db_dir))}

        eggnog_db_gz = db_dir / "eggnog.db.gz"
        taxa_tar_gz = db_dir / "eggnog.taxa.tar.gz"
        dmnd_gz = db_dir / "eggnog_proteins.dmnd.gz"

        url_db = f"{base_url}/eggnog.db.gz"
        url_taxa = f"{base_url}/eggnog.taxa.tar.gz"
        url_dmnd = f"{base_url}/eggnog_proteins.dmnd.gz"

        # retries fixed here (not configurable)
        _download_with_retries(url_db, eggnog_db_gz, retries=5, sleep_sec=5)
        _download_with_retries(url_taxa, taxa_tar_gz, retries=5, sleep_sec=5)
        _download_with_retries(url_dmnd, dmnd_gz, retries=5, sleep_sec=5)

        # Sanity floors (catch empty/corrupt)
        _assert_min_size(taxa_tar_gz, 1_000_000, "eggnog.taxa.tar.gz")
        _assert_min_size(eggnog_db_gz, 1_000_000, "eggnog.db.gz")
        _assert_min_size(dmnd_gz, 1_000_000, "eggnog_proteins.dmnd.gz")

        # Extract once; delete archives
        if not (db_dir / "eggnog.db").exists():
            _gunzip(eggnog_db_gz, db_dir / "eggnog.db")
        eggnog_db_gz.unlink(missing_ok=True)

        if not (db_dir / "eggnog.taxa.db").exists() or not (db_dir / "eggnog.taxa.db.traverse.pkl").exists():
            if not tarfile.is_tarfile(taxa_tar_gz):
                raise RuntimeError(f"taxa archive is not a valid tar file: {taxa_tar_gz}")
            _untar_gz(taxa_tar_gz, db_dir)
        taxa_tar_gz.unlink(missing_ok=True)

        if not (db_dir / "eggnog_proteins.dmnd").exists():
            _gunzip(dmnd_gz, db_dir / "eggnog_proteins.dmnd")
        dmnd_gz.unlink(missing_ok=True)

        missing = [p.name for p in required if not p.exists()]
        if missing:
            raise RuntimeError(f"DB extraction finished but missing: {missing}")

        print(f"[INFO] eggNOG DB ready: {db_dir}", flush=True)
        return {"db_dir": Folder(str(db_dir))}
