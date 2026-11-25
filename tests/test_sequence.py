import os
import sys
import unittest

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from backend.services.sequence_service import SequenceService


class TestSequenceProcessing(unittest.TestCase):

    def test_fasta_validation_valid(self):
        seq = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        res = SequenceService.validate_sequence(seq)
        self.assertTrue(res["is_valid"])
        self.assertEqual(res["length"], 60)
        self.assertEqual(len(res["errors"]), 0)

    def test_fasta_validation_invalid_characters(self):
        seq = "MVHLTPEEK123$$$ZZZ"
        res = SequenceService.validate_sequence(seq)
        self.assertFalse(res["is_valid"])
        self.assertTrue(any("Invalid non-amino acid" in err for err in res["errors"]))

    def test_fasta_validation_empty(self):
        res = SequenceService.validate_sequence("   \n\t  ")
        self.assertFalse(res["is_valid"])
        self.assertIn("Sequence cannot be empty.", res["errors"])

    def test_single_substitution(self):
        ref = "MVHLTPEEKS"
        mut = "MVHLTPVEKS" # E7V (precursor)
        aln = SequenceService.align_sequences(ref, mut, protein_name="CUSTOM")
        self.assertEqual(len(aln["mutations"]), 1)
        m = aln["mutations"][0]
        self.assertEqual(m["mutation_type"], "SUBSTITUTION")
        self.assertEqual(m["ref_pos"], 7)
        self.assertEqual(m["ref_aa"], "E")
        self.assertEqual(m["mut_aa"], "V")
        self.assertEqual(m["canonical_notation"], "E7V")

    def test_mature_numbering_hbb(self):
        ref = "MVHLTPEEKS"
        mut = "MVHLTPVEKS"
        aln = SequenceService.align_sequences(ref, mut, protein_name="HBB")
        self.assertEqual(len(aln["mutations"]), 1)
        m = aln["mutations"][0]
        # In HBB, mature offset is -1
        self.assertEqual(m["mature_pos"], 6)
        self.assertEqual(m["mature_notation"], "E6V")

    def test_deletion_mutation(self):
        ref = "MQRSPLEKASVVSKLWFSWTRPILR"
        # Delete W at position 16: MQRSPLEKASVVSKL[W]FSWTRPILR -> MQRSPLEKASVVSKLFSWTRPILR
        mut = "MQRSPLEKASVVSKLFSWTRPILR"
        aln = SequenceService.align_sequences(ref, mut, protein_name="CFTR")
        self.assertEqual(len(aln["mutations"]), 1)
        m = aln["mutations"][0]
        self.assertEqual(m["mutation_type"], "DELETION")
        self.assertEqual(m["ref_pos"], 16)
        self.assertEqual(m["ref_aa"], "W")
        self.assertIsNone(m["mut_aa"])
        self.assertEqual(m["canonical_notation"], "W16del")

    def test_insertion_mutation(self):
        ref = "MVHLTPEEKS"
        # Insert G between pos 5 and 6: MVHLT[G]PEEKS
        mut = "MVHLTGPEEKS"
        aln = SequenceService.align_sequences(ref, mut, protein_name="CUSTOM")
        self.assertEqual(len(aln["mutations"]), 1)
        m = aln["mutations"][0]
        self.assertEqual(m["mutation_type"], "INSERTION")
        self.assertEqual(m["mut_aa"], "G")

    def test_multiple_substitutions(self):
        ref = "ACDEFGHIKL"
        mut = "ACREFGHIKM" # D3R and L10M
        aln = SequenceService.align_sequences(ref, mut, protein_name="CUSTOM")
        self.assertEqual(len(aln["mutations"]), 2)
        notations = [m["canonical_notation"] for m in aln["mutations"]]
        self.assertIn("D3R", notations)
        self.assertIn("L10M", notations)

    def test_benchmark_cftr_f508del(self):
        ref_p = os.path.join(BASE_DIR, "data", "fasta", "CFTR_reference.fasta")
        mut_p = os.path.join(BASE_DIR, "data", "fasta", "CFTR_F508del_mutant.fasta")
        if os.path.exists(ref_p) and os.path.exists(mut_p):
            r_seq = SequenceService.parse_fasta(ref_p)[0]["sequence"]
            m_seq = SequenceService.parse_fasta(mut_p)[0]["sequence"]
            aln = SequenceService.align_sequences(r_seq, m_seq, protein_name="CFTR")
            self.assertEqual(len(aln["mutations"]), 1)
            self.assertEqual(aln["mutations"][0]["canonical_notation"], "F508del")
            self.assertEqual(aln["mutations"][0]["mutation_type"], "DELETION")


if __name__ == "__main__":
    unittest.main()
