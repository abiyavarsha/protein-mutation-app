import os
import sys
import unittest

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from backend.main import health_check, list_proteins, validate_sequence, align_sequences, detect_mutations
from backend.models.schemas import ValidationRequest, AlignmentRequest, MutationDetectionRequest
from backend.database import SessionLocal, init_db


class TestFastAPIEndpoints(unittest.TestCase):

    def setUp(self):
        init_db()
        self.db = SessionLocal()

    def tearDown(self):
        self.db.close()

    def test_health_endpoint(self):
        res = health_check()
        self.assertEqual(res["status"], "healthy")
        self.assertEqual(res["database"], "connected")

    def test_list_proteins_endpoint(self):
        proteins = list_proteins(self.db)
        names = [p["name"] for p in proteins]
        self.assertIn("HBB", names)
        self.assertIn("CFTR", names)
        self.assertIn("EGFR", names)
        self.assertIn("TP53", names)

    def test_validate_endpoint_valid(self):
        req = ValidationRequest(sequence="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK")
        res = validate_sequence(req)
        self.assertTrue(res.is_valid)
        self.assertEqual(res.length, 60)

    def test_validate_endpoint_invalid(self):
        req = ValidationRequest(sequence="MVHLT123$#@")
        res = validate_sequence(req)
        self.assertFalse(res.is_valid)

    def test_align_endpoint(self):
        req = AlignmentRequest(
            ref_sequence="MVHLTPEEKS",
            mut_sequence="MVHLTPVEKS",
            protein_name="HBB"
        )
        res = align_sequences(req)
        self.assertAlmostEqual(res.sequence_identity, 90.0, places=1)
        self.assertEqual(len(res.mutations), 1)
        self.assertEqual(res.mutations[0].mature_notation, "E6V")

    def test_mutations_detection_endpoint(self):
        req = MutationDetectionRequest(
            ref_sequence="MQRSPLEKASVVSKLFFSWTRPILR",
            mut_sequence="MQRSPLEKASVVSKLFSWTRPILR",
            protein_name="CFTR"
        )
        res = detect_mutations(req)
        self.assertEqual(len(res["mutations_detected"]), 1)
        self.assertEqual(res["mutations_detected"][0]["mutation_type"], "DELETION")


if __name__ == "__main__":
    unittest.main()
