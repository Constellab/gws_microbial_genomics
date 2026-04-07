import os

from gws_core import BaseTestCase, Folder, ResourceSet, TaskRunner

from gws_microbial_genomics.Bactopia.bactopia_explore import BactopiaExplore

TESTDATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "../testdata"))
MOCK_BACTOPIA_OUT = os.path.join(TESTDATA, "mock_bactopia_out")


class TestBactopiaExplore(BaseTestCase):

    def test_explore_mock_output(self):
        """BactopiaExplore correctly classifies mock Bactopia outputs into ResourceSets."""
        self.print("Test BactopiaExplore: mock bactopia output folder")

        if not os.path.isdir(MOCK_BACTOPIA_OUT):
            self.skipTest(f"Mock folder not found: {MOCK_BACTOPIA_OUT} — run test_bactopia_run first")

        outputs = TaskRunner(
            task_type=BactopiaExplore,
            inputs={"bactopia_results": Folder(MOCK_BACTOPIA_OUT)},
            params={},
        ).run()

        # All 8 output keys must be present
        expected_keys = ("qc", "assembler", "annotator", "gather", "mlst",
                         "amrfinderplus", "bactopia_runs", "sketcher")
        for key in expected_keys:
            self.assertIn(key, outputs, msg=f"Missing output key: {key}")
            self.assertIsNotNone(outputs[key])

        # Log all outputs for diagnosis
        for key in expected_keys:
            rs: ResourceSet = outputs[key]
            self.print(f"{key} resources: {list(rs.get_resources().keys())}")

        # Gather: 2 TSV files (sample1 + sample2)
        gather: ResourceSet = outputs["gather"]
        self.assertGreaterEqual(len(gather.get_resources()), 2, msg="Expected ≥2 gather TSV files")

        # bactopia_runs: HTML reports + DAG + meta TSV
        br: ResourceSet = outputs["bactopia_runs"]
        self.assertGreaterEqual(len(br.get_resources()), 3, msg="Expected ≥3 resources in bactopia_runs")
