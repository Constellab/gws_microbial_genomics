import os

import pandas as pd
from gws_core import BaseTestCase, File, Folder, TaskRunner

from gws_microbial_genomics.Bactopia.bactopiaMetadata import BactopiaMetadata

TESTDATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "../testdata"))
FASTQ_DIR = os.path.join(TESTDATA, "fastq")


class TestBactopiaMetadata(BaseTestCase):

    def test_prepare_pe_fastq(self):
        """BactopiaMetadata generates a samplesheet from a PE FASTQ folder."""
        self.print("Test BactopiaMetadata: sample1_R1/R2 paired-end")

        samplesheet: File = TaskRunner(
            task_type=BactopiaMetadata,
            inputs={"fastq_folder": Folder(FASTQ_DIR)},
            params={
                "species": "",
                "genome_size": "0",
                "fastq_separator": "_",
                "pe1_pattern": "R1",
                "pe2_pattern": "R2",
                "ont": "false",
            },
        ).run()["samplesheet"]

        self.assertIsNotNone(samplesheet)
        self.assertTrue(os.path.exists(samplesheet.path))

        df = pd.read_csv(samplesheet.path, sep="\t")
        self.print(f"Samplesheet:\n{df.to_string()}")

        self.assertIn("sample", df.columns)
        self.assertIn("runtype", df.columns)
        self.assertIn("r1", df.columns)
        self.assertGreater(len(df), 0, msg="Expected at least one row in samplesheet")

    def test_prepare_se_fastq(self):
        """BactopiaMetadata generates a samplesheet from a SE FASTQ folder."""
        self.print("Test BactopiaMetadata: sample2_SE single-end")

        # Use a subfolder with only SE file for clarity
        import tempfile, shutil, gzip
        with tempfile.TemporaryDirectory() as tmpdir:
            # copy only SE file
            src = os.path.join(FASTQ_DIR, "sample2.fastq.gz")
            dst = os.path.join(tmpdir, "sample2.fastq.gz")
            shutil.copy(src, dst)

            samplesheet: File = TaskRunner(
                task_type=BactopiaMetadata,
                inputs={"fastq_folder": Folder(tmpdir)},
                params={
                    "species": "",
                    "genome_size": "0",
                    "fastq_separator": "_",
                    "pe1_pattern": "R1",
                    "pe2_pattern": "R2",
                    "ont": "false",
                },
            ).run()["samplesheet"]

        self.assertIsNotNone(samplesheet)
        self.assertTrue(os.path.exists(samplesheet.path))

        df = pd.read_csv(samplesheet.path, sep="\t")
        self.print(f"SE Samplesheet:\n{df.to_string()}")
        self.assertGreater(len(df), 0, msg="Expected at least one row for SE sample")
