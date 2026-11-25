import re
import io
import csv
import json
from typing import List, Dict, Any, Tuple, Optional
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# Mature chain offsets (e.g. HBB initiator Met is cleaved in mature hemoglobin, mature E6 is precursor E7)
try:
    from backend.dataset_registry import get_mature_offsets
    MATURE_OFFSETS = get_mature_offsets()
except Exception:
    MATURE_OFFSETS = {
        "HBB": -1,
        "CFTR": 0,
        "EGFR": 0,
        "TP53": 0,
        "SOD1": -1,
        "F9": -46,
        "GBA1": -39
    }


class SequenceService:

    @staticmethod
    def validate_sequence(sequence_str: str) -> Dict[str, Any]:
        """
        Validates an amino acid sequence string.
        Removes whitespace and checks for standard amino acid characters.
        """
        errors = []
        warnings = []
        
        # Clean sequence
        cleaned = re.sub(r"\s+", "", sequence_str).upper()
        
        if not cleaned:
            errors.append("Sequence cannot be empty.")
            return {
                "is_valid": False,
                "length": 0,
                "cleaned_sequence": "",
                "errors": errors,
                "warnings": warnings
            }
            
        # Check characters
        invalid_chars = set(cleaned) - STANDARD_AMINO_ACIDS
        if invalid_chars:
            # Check if common non-standard or stop codon
            if "*" in invalid_chars:
                warnings.append("Sequence contains stop codon '*'.")
            if "X" in invalid_chars:
                warnings.append("Sequence contains unknown amino acid code 'X'.")
            if "U" in invalid_chars:
                warnings.append("Sequence contains selenocysteine 'U'.")
            other_invalid = invalid_chars - {"*", "X", "U"}
            if other_invalid:
                errors.append(f"Invalid non-amino acid characters found: {sorted(list(other_invalid))}")
                
        if len(cleaned) < 10:
            warnings.append(f"Sequence is very short ({len(cleaned)} residues). Typical proteins are >30 residues.")
            
        return {
            "is_valid": len(errors) == 0,
            "length": len(cleaned),
            "cleaned_sequence": cleaned,
            "errors": errors,
            "warnings": warnings
        }

    @staticmethod
    def parse_fasta(fasta_content_or_path: str) -> List[Dict[str, str]]:
        """
        Parses FASTA format from file path or text string.
        Returns list of dicts: [{'header': ..., 'sequence': ..., 'description': ...}]
        """
        records = []
        is_path = False
        try:
            import os
            if os.path.exists(fasta_content_or_path):
                is_path = True
        except Exception:
            pass
            
        if is_path:
            with open(fasta_content_or_path, "r") as f:
                for rec in SeqIO.parse(f, "fasta"):
                    records.append({
                        "header": rec.id,
                        "description": rec.description,
                        "sequence": str(rec.seq).upper().strip()
                    })
        else:
            handle = io.StringIO(fasta_content_or_path)
            for rec in SeqIO.parse(handle, "fasta"):
                records.append({
                    "header": rec.id,
                    "description": rec.description,
                    "sequence": str(rec.seq).upper().strip()
                })
                
        # If no FASTA headers found but raw sequence text provided
        if not records and fasta_content_or_path.strip():
            raw_seq = re.sub(r"\s+", "", fasta_content_or_path).upper()
            records.append({
                "header": "custom_sequence",
                "description": "User provided sequence",
                "sequence": raw_seq
            })
            
        return records

    @staticmethod
    def align_sequences(
        ref_seq: str,
        mut_seq: str,
        protein_name: str = "CUSTOM"
    ) -> Dict[str, Any]:
        """
        Performs global pairwise sequence alignment using Biopython PairwiseAligner.
        Returns alignment metrics, aligned strings, and identified mutations.
        """
        val_ref = SequenceService.validate_sequence(ref_seq)
        val_mut = SequenceService.validate_sequence(mut_seq)
        
        if not val_ref["is_valid"]:
            raise ValueError(f"Reference sequence invalid: {'; '.join(val_ref['errors'])}")
        if not val_mut["is_valid"]:
            raise ValueError(f"Mutant sequence invalid: {'; '.join(val_mut['errors'])}")
            
        r_seq = val_ref["cleaned_sequence"]
        m_seq = val_mut["cleaned_sequence"]
        
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -2.0
        aligner.extend_gap_score = -0.5
        
        alignments = aligner.align(r_seq, m_seq)
        if not alignments:
            raise RuntimeError("No alignment could be computed between sequences.")
            
        best_aln = alignments[0]
        aligned_ref = str(best_aln[0])
        aligned_mut = str(best_aln[1])
        
        # Calculate matches, mismatches, gaps
        matches = 0
        mismatches = 0
        gaps = 0
        
        for r_ch, m_ch in zip(aligned_ref, aligned_mut):
            if r_ch == '-' or m_ch == '-':
                gaps += 1
            elif r_ch == m_ch:
                matches += 1
            else:
                mismatches += 1
                
        total_aln_len = len(aligned_ref)
        seq_identity = (matches / max(len(r_seq), len(m_seq))) * 100.0 if max(len(r_seq), len(m_seq)) > 0 else 0.0
        similarity = (matches / (matches + mismatches)) * 100.0 if (matches + mismatches) > 0 else 0.0
        
        # Detect mutations
        mutations = SequenceService.detect_mutations_from_alignment(
            aligned_ref=aligned_ref,
            aligned_mut=aligned_mut,
            protein_name=protein_name
        )
        
        # Format human-readable visual alignment snippet
        alignment_repr = SequenceService.format_alignment_visual(aligned_ref, aligned_mut)
        
        return {
            "ref_length": len(r_seq),
            "mut_length": len(m_seq),
            "alignment_score": float(best_aln.score),
            "sequence_identity": round(seq_identity, 2),
            "similarity_percentage": round(similarity, 2),
            "gaps_count": gaps,
            "matches_count": matches,
            "mismatches_count": mismatches,
            "aligned_ref": aligned_ref,
            "aligned_mut": aligned_mut,
            "alignment_repr": alignment_repr,
            "mutations": mutations
        }

    @staticmethod
    def detect_mutations_from_alignment(
        aligned_ref: str,
        aligned_mut: str,
        protein_name: str = "CUSTOM"
    ) -> List[Dict[str, Any]]:
        """
        Traces residue-by-residue correspondence to detect:
        - Substitutions
        - Deletions
        - Insertions
        Applies mature-protein numbering (e.g. HBB E7V -> mature E6V).
        """
        mature_offset = MATURE_OFFSETS.get(protein_name.upper(), 0)
        
        mutations = []
        ref_idx = 0  # 1-based index in reference
        mut_idx = 0  # 1-based index in mutant
        
        i = 0
        n = len(aligned_ref)
        
        while i < n:
            r_ch = aligned_ref[i]
            m_ch = aligned_mut[i]
            
            if r_ch != '-':
                ref_idx += 1
            if m_ch != '-':
                mut_idx += 1
                
            if r_ch != m_ch:
                # 1. Insertion (gap in ref, character in mut)
                if r_ch == '-':
                    ins_start_ref = ref_idx
                    ins_start_mut = mut_idx
                    ins_seq = [m_ch]
                    j = i + 1
                    while j < n and aligned_ref[j] == '-':
                        mut_idx += 1
                        ins_seq.append(aligned_mut[j])
                        j += 1
                    ins_str = "".join(ins_seq)
                    canonical_notation = f"ins_{ins_start_ref}_{ins_start_ref+1}{ins_str}"
                    mutations.append({
                        "mutation_type": "INSERTION",
                        "ref_pos": ins_start_ref,
                        "mut_pos": ins_start_mut,
                        "ref_aa": None,
                        "mut_aa": ins_str,
                        "canonical_notation": canonical_notation,
                        "mature_pos": ins_start_ref + mature_offset if ins_start_ref + mature_offset > 0 else None,
                        "mature_notation": canonical_notation,
                        "description": f"Insertion of {ins_str} between reference residues {ins_start_ref} and {ins_start_ref+1}"
                    })
                    i = j
                    continue
                    
                # 2. Deletion (character in ref, gap in mut)
                elif m_ch == '-':
                    del_start_ref = ref_idx
                    del_seq = [r_ch]
                    j = i + 1
                    while j < n and aligned_mut[j] == '-':
                        ref_idx += 1
                        del_seq.append(aligned_ref[j])
                        j += 1
                    del_str = "".join(del_seq)
                    if len(del_seq) == 1:
                        canonical_notation = f"{del_str}{del_start_ref}del"
                        mature_pos = del_start_ref + mature_offset
                        mature_notation = f"{del_str}{mature_pos}del" if mature_offset != 0 else canonical_notation
                        desc = f"Deletion of {del_str} at reference position {del_start_ref}"
                    else:
                        end_del_ref = del_start_ref + len(del_seq) - 1
                        canonical_notation = f"{del_str[0]}{del_start_ref}_{del_str[-1]}{end_del_ref}del"
                        mature_pos = del_start_ref + mature_offset
                        mature_notation = f"{del_str[0]}{del_start_ref + mature_offset}_{del_str[-1]}{end_del_ref + mature_offset}del" if mature_offset != 0 else canonical_notation
                        desc = f"Deletion of {del_str} from position {del_start_ref} to {end_del_ref}"
                        
                    mutations.append({
                        "mutation_type": "DELETION",
                        "ref_pos": del_start_ref,
                        "mut_pos": mut_idx,
                        "ref_aa": del_str,
                        "mut_aa": None,
                        "canonical_notation": canonical_notation,
                        "mature_pos": mature_pos,
                        "mature_notation": mature_notation,
                        "description": desc
                    })
                    i = j
                    continue
                    
                # 3. Substitution
                else:
                    canonical_notation = f"{r_ch}{ref_idx}{m_ch}"
                    mature_pos = ref_idx + mature_offset
                    mature_notation = f"{r_ch}{mature_pos}{m_ch}" if mature_offset != 0 else canonical_notation
                    
                    desc = f"Substitution of {r_ch} to {m_ch} at position {ref_idx}"
                    if mature_offset != 0:
                        desc += f" (mature chain position: {mature_pos}, notation: {mature_notation})"
                        
                    mutations.append({
                        "mutation_type": "SUBSTITUTION",
                        "ref_pos": ref_idx,
                        "mut_pos": mut_idx,
                        "ref_aa": r_ch,
                        "mut_aa": m_ch,
                        "canonical_notation": canonical_notation,
                        "mature_pos": mature_pos,
                        "mature_notation": mature_notation,
                        "description": desc
                    })
            i += 1
            
        return mutations

    @staticmethod
    def format_alignment_visual(aligned_ref: str, aligned_mut: str, chunk_size: int = 60) -> str:
        """
        Formats visual pairwise alignment with match symbols (| for match, : for mismatch, space for gap).
        """
        match_line = []
        for r, m in zip(aligned_ref, aligned_mut):
            if r == '-' or m == '-':
                match_line.append(' ')
            elif r == m:
                match_line.append('|')
            else:
                match_line.append('.')
        match_str = "".join(match_line)
        
        blocks = []
        for start in range(0, len(aligned_ref), chunk_size):
            end = min(start + chunk_size, len(aligned_ref))
            ref_sub = aligned_ref[start:end]
            mid_sub = match_str[start:end]
            mut_sub = aligned_mut[start:end]
            blocks.append(f"Ref  {start+1:5d}  {ref_sub}  {end:5d}\n              {mid_sub}\nMut  {start+1:5d}  {mut_sub}  {end:5d}\n")
            
        return "\n".join(blocks)

    @staticmethod
    def generate_csv_summary(mutations: List[Dict[str, Any]], protein_name: str) -> str:
        """
        Generates CSV format string summarizing detected mutations.
        """
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow([
            "Protein",
            "Mutation_Type",
            "Canonical_Notation",
            "Mature_Notation",
            "Ref_Pos",
            "Mature_Pos",
            "Mut_Pos",
            "Ref_AA",
            "Mut_AA",
            "Description"
        ])
        
        for m in mutations:
            writer.writerow([
                protein_name,
                m.get("mutation_type"),
                m.get("canonical_notation"),
                m.get("mature_notation"),
                m.get("ref_pos"),
                m.get("mature_pos"),
                m.get("mut_pos"),
                m.get("ref_aa") or "-",
                m.get("mut_aa") or "-",
                m.get("description")
            ])
            
        return output.getvalue()

    @staticmethod
    def validate_mutation_label(expected_label: str, detected_mutations: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Checks if detected mutations match the expected mutation label (e.g. 'F508del', 'E6V', 'L858R').
        """
        if not detected_mutations:
            return {
                "matches": False,
                "message": f"No mutations detected, but expected '{expected_label}'."
            }
            
        label_clean = expected_label.strip()
        matched = False
        matched_item = None
        
        for m in detected_mutations:
            if label_clean in (m.get("canonical_notation", ""), m.get("mature_notation", "")):
                matched = True
                matched_item = m
                break
            # Also check without 'del' or case variations
            if label_clean.replace("del", "").lower() in m.get("canonical_notation", "").lower():
                matched = True
                matched_item = m
                break
                
        if matched:
            return {
                "matches": True,
                "message": f"Detected mutation {matched_item.get('canonical_notation')} (mature: {matched_item.get('mature_notation')}) successfully verified against expected label '{expected_label}'."
            }
        else:
            found_notations = [m.get("canonical_notation") for m in detected_mutations]
            return {
                "matches": False,
                "message": f"Expected mutation '{expected_label}' did not directly match detected mutations: {', '.join(found_notations)}."
            }
