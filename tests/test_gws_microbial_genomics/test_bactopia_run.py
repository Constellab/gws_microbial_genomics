import os

from gws_core import BaseTestCase, File, Folder, ResourceSet, TaskRunner

from gws_microbial_genomics.Bactopia.bactopia_explore import BactopiaExplore
from gws_microbial_genomics.Bactopia.bactopiaMetadata import BactopiaMetadata
from gws_microbial_genomics.Bactopia.bactopia_run import BactopiaRun

TESTDATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "../testdata"))
FASTQ_DIR = os.path.join(TESTDATA, "fastq")


class TestBactopiaRun(BaseTestCase):

    def test_bactopia_run_pe(self):
        """Chain BactopiaMetadata → BactopiaRun on a PE sample (~10 min expected)."""
        self.print("Test BactopiaRun: BactopiaMetadata → BactopiaRun (PE)")

        # Step 1: generate samplesheet
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

        self.print(f"Samplesheet: {samplesheet.path}")

        # Step 2: run Bactopia
        bactopia_out: Folder = TaskRunner(
            task_type=BactopiaRun,
            inputs={
                "samplesheet": samplesheet,
                "fastq_folder": Folder(FASTQ_DIR),
            },
            params={
                "max_cpus": "4",
                "coverage": "100",
            },
        ).run()["bactopia_results"]

        self.assertIsNotNone(bactopia_out)
        self.assertTrue(os.path.isdir(bactopia_out.path), msg="Bactopia output folder should exist")
        output_files = list(os.walk(bactopia_out.path))
        self.assertGreater(len(output_files), 0, msg="Expected files in Bactopia output folder")
        self.print(f"Bactopia output top-level: {os.listdir(bactopia_out.path)}")

        # Step 3: explore the Bactopia output
        self.print("Step 3: BactopiaExplore on the output")
        explore_outputs = TaskRunner(
            task_type=BactopiaExplore,
            inputs={"bactopia_results": bactopia_out},
            params={},
        ).run()

        expected_keys = ("qc", "assembler", "annotator", "gather", "mlst",
                         "amrfinderplus", "bactopia_runs", "sketcher")
        for key in expected_keys:
            self.assertIn(key, explore_outputs, msg=f"Missing output key: {key}")
            rs: ResourceSet = explore_outputs[key]
            self.print(f"{key} resources: {list(rs.get_resources().keys())}")
