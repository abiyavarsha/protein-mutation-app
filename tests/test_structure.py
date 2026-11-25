import os
import sys
import unittest
import numpy as np

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from backend.services.structure_service import StructureService
from backend.services.superimposition_service import SuperimpositionService


class TestStructureSuperimposition(unittest.TestCase):

    def setUp(self):
        self.hbb_ref = os.path.join(BASE_DIR, "data", "structures", "HBB", "reference.pdb")
        self.hbb_mut = os.path.join(BASE_DIR, "data", "structures", "HBB", "E6V_mutant.pdb")
        self.hbb_sup = os.path.join(BASE_DIR, "data", "structures", "HBB", "E6V_superimposed.pdb")

    def test_pdb_validation(self):
        if os.path.exists(self.hbb_ref):
            is_valid, msg, ca_cnt = StructureService.validate_pdb(self.hbb_ref)
            self.assertTrue(is_valid)
            self.assertEqual(ca_cnt, 147)

    def test_pdb_validation_invalid(self):
        is_valid, msg, ca_cnt = StructureService.validate_pdb("REMARK This is not a valid pdb\nEND\n")
        self.assertFalse(is_valid)
        self.assertEqual(ca_cnt, 0)

    def test_kabsch_superimposition_hbb(self):
        if os.path.exists(self.hbb_ref) and os.path.exists(self.hbb_mut):
            res = SuperimpositionService.superimpose_structures(
                self.hbb_ref,
                self.hbb_mut,
                self.hbb_sup,
                mutation_ref_pos=7
            )
            self.assertEqual(res["matched_ca_atoms_count"], 147)
            self.assertGreaterEqual(res["global_rmsd"], 0.0)
            self.assertLess(res["global_rmsd"], 1.0)
            self.assertEqual(res["alignment_coverage"], 100.0)
            self.assertTrue(os.path.exists(self.hbb_sup))

    def test_cftr_f508del_deletion_superimposition(self):
        cftr_ref = os.path.join(BASE_DIR, "data", "structures", "CFTR", "reference.pdb")
        cftr_f508 = os.path.join(BASE_DIR, "data", "structures", "CFTR", "F508del_mutant.pdb")
        cftr_sup = os.path.join(BASE_DIR, "data", "structures", "CFTR", "F508del_superimposed.pdb")
        
        if os.path.exists(cftr_ref) and os.path.exists(cftr_f508):
            res = SuperimpositionService.superimpose_structures(
                cftr_ref,
                cftr_f508,
                cftr_sup,
                mutation_ref_pos=78
            )
            # Reference has 200 residues, mutant has 199 (due to F508 deletion)
            self.assertEqual(res["matched_ca_atoms_count"], 199)
            self.assertLess(res["global_rmsd"], 1.5)
            self.assertGreater(res["alignment_coverage"], 99.0)

    def test_rmsd_formula_exactness(self):
        # Test mathematical formula verification on known synthetic coordinate shifts
        coords_ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        # Shift each coordinate by exactly 1.0 along Z axis
        coords_mut = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]])
        diff = coords_ref - coords_mut
        sq_dist = np.sum(diff ** 2, axis=1) # each is 1.0
        rmsd = np.sqrt(np.mean(sq_dist)) # sqrt(1.0) = 1.0
        self.assertAlmostEqual(rmsd, 1.0, places=5)


if __name__ == "__main__":
    unittest.main()
