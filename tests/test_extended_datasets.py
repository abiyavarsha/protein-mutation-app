import unittest
import os
import sys

# Ensure app path is in sys.path
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from backend.dataset_registry import (
    DATASET_REGISTRY,
    get_all_proteins,
    get_protein,
    get_mutants_for_protein,
    get_mature_offsets,
    KNOWN_BIOLOGICAL_CONTEXT
)
from backend.services.sequence_service import SequenceService


class TestExtendedDatasets(unittest.TestCase):

    def test_registry_coverage(self):
        """Verifies that all 28 genes are present and well-formed in the registry."""
        self.assertEqual(len(DATASET_REGISTRY), 28)
        
        # Verify 4 benchmark proteins are flagged as Original Benchmark
        benchmark_genes = ["HBB", "CFTR", "EGFR", "TP53"]
        for bg in benchmark_genes:
            self.assertIn(bg, DATASET_REGISTRY)
            self.assertEqual(DATASET_REGISTRY[bg]["category"], "Original Benchmark")
            
        # Verify total registered mutations
        total_muts = sum(len(info.get("mutants", [])) for info in DATASET_REGISTRY.values())
        self.assertGreaterEqual(total_muts, 49)

    def test_mature_offsets(self):
        """Verifies mature offsets for proteins with cleavage."""
        offsets = get_mature_offsets()
        self.assertEqual(offsets["HBB"], -1)
        self.assertEqual(offsets["SOD1"], -1)
        self.assertEqual(offsets["F9"], -46)
        self.assertEqual(offsets["GBA1"], -39)
        self.assertEqual(offsets["EGFR"], 0)
        self.assertEqual(offsets["CFTR"], 0)

    def test_kras_g12d_alignment_and_mutation(self):
        """Tests alignment and mutation calling on extended dataset KRAS G12D."""
        fasta_dir = os.path.join(BASE_DIR, "data", "fasta")
        ref_path = os.path.join(fasta_dir, "KRAS_reference.fasta")
        mut_path = os.path.join(fasta_dir, "KRAS_G12D_mutant.fasta")
        
        self.assertTrue(os.path.exists(ref_path), "KRAS reference FASTA must exist")
        self.assertTrue(os.path.exists(mut_path), "KRAS G12D mutant FASTA must exist")
        
        ref_recs = SequenceService.parse_fasta(ref_path)
        mut_recs = SequenceService.parse_fasta(mut_path)
        
        res = SequenceService.align_sequences(
            ref_recs[0]["sequence"],
            mut_recs[0]["sequence"],
            "KRAS"
        )
        mutations = res["mutations"]
        self.assertEqual(len(mutations), 1)
        m = mutations[0]
        self.assertEqual(m["ref_aa"], "G")
        self.assertEqual(m["mut_aa"], "D")
        self.assertEqual(m["ref_pos"], 12)
        self.assertEqual(m["canonical_notation"], "G12D")

    def test_braf_v600e_alignment_and_mutation(self):
        """Tests alignment and mutation calling on extended dataset BRAF V600E."""
        fasta_dir = os.path.join(BASE_DIR, "data", "fasta")
        ref_path = os.path.join(fasta_dir, "BRAF_reference.fasta")
        mut_path = os.path.join(fasta_dir, "BRAF_V600E_mutant.fasta")
        
        ref_recs = SequenceService.parse_fasta(ref_path)
        mut_recs = SequenceService.parse_fasta(mut_path)
        
        res = SequenceService.align_sequences(
            ref_recs[0]["sequence"],
            mut_recs[0]["sequence"],
            "BRAF"
        )
        mutations = res["mutations"]
        self.assertEqual(len(mutations), 1)
        m = mutations[0]
        self.assertEqual(m["ref_aa"], "V")
        self.assertEqual(m["mut_aa"], "E")
        self.assertEqual(m["ref_pos"], 600)
        self.assertEqual(m["canonical_notation"], "V600E")

    def test_sod1_a4v_mature_numbering(self):
        """Tests mature chain numbering offset on SOD1 A4V (precursor Ala5Val)."""
        fasta_dir = os.path.join(BASE_DIR, "data", "fasta")
        ref_path = os.path.join(fasta_dir, "SOD1_reference.fasta")
        mut_path = os.path.join(fasta_dir, "SOD1_A4V_mutant.fasta")
        
        ref_recs = SequenceService.parse_fasta(ref_path)
        mut_recs = SequenceService.parse_fasta(mut_path)
        
        res = SequenceService.align_sequences(
            ref_recs[0]["sequence"],
            mut_recs[0]["sequence"],
            "SOD1"
        )
        mutations = res["mutations"]
        self.assertEqual(len(mutations), 1)
        m = mutations[0]
        self.assertEqual(m["ref_aa"], "A")
        self.assertEqual(m["mut_aa"], "V")
        self.assertEqual(m["ref_pos"], 5)   # Uniprot precursor position
        self.assertEqual(m["mature_pos"], 4)      # Mature protein position (offset -1)
        self.assertEqual(m["mature_notation"], "A4V")

    def test_ryr1_documented_discrepancy(self):
        """Verifies that RYR1 is flagged as unavailable and documented properly."""
        ryr1 = get_protein("RYR1")
        self.assertIsNotNone(ryr1)
        muts = ryr1.get("mutants", [])
        self.assertEqual(len(muts), 1)
        self.assertFalse(muts[0].get("available", True))
        self.assertIn("P21817", muts[0]["clinical_note"])


if __name__ == "__main__":
    unittest.main()
