import os
import math
import numpy as np
from typing import Dict, Any, List, Tuple, Optional
from Bio.PDB import PDBParser, PDBIO, Superimposer, Structure, Model, Chain, Residue, Atom
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner


class SuperimpositionService:

    @staticmethod
    def load_structure(pdb_path: str, structure_id: str = "protein") -> Structure.Structure:
        """Loads a PDB structure using Bio.PDB.PDBParser."""
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"PDB file not found at: {pdb_path}")
        parser = PDBParser(QUIET=True)
        return parser.get_structure(structure_id, pdb_path)

    @staticmethod
    def extract_ca_atoms_and_sequence(structure: Structure.Structure) -> Tuple[List[Atom.Atom], List[Residue.Residue], str]:
        """
        Extracts C-alpha atoms, residues containing CA, and 1-letter amino acid sequence from structure.
        """
        ca_atoms = []
        ca_residues = []
        seq_chars = []
        
        # Take the first model
        model = structure[0]
        for chain in model:
            for residue in chain:
                # Check for standard amino acid and CA atom
                if residue.has_id("CA") and residue.get_id()[0] == " ":
                    res_name = residue.get_resname()
                    try:
                        one_letter = seq1(res_name)
                    except Exception:
                        one_letter = "X"
                    ca_atoms.append(residue["CA"])
                    ca_residues.append(residue)
                    seq_chars.append(one_letter)
                    
        return ca_atoms, ca_residues, "".join(seq_chars)

    @staticmethod
    def map_residues_by_sequence_alignment(
        ref_residues: List[Residue.Residue],
        mut_residues: List[Residue.Residue],
        ref_seq: str,
        mut_seq: str
    ) -> List[Tuple[Residue.Residue, Residue.Residue, int, int]]:
        """
        Aligns the actual sequences extracted from reference and mutant PDB structures
        using Biopython PairwiseAligner to find exact residue-to-residue correspondences.
        Excludes gaps and deletion/insertion sites.
        Returns list of (ref_res, mut_res, ref_idx_in_pdb, mut_idx_in_pdb).
        """
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -2.0
        aligner.extend_gap_score = -0.5
        
        alignments = aligner.align(ref_seq, mut_seq)
        if not alignments:
            raise RuntimeError("Could not align sequence from PDB files for residue mapping.")
            
        best_aln = alignments[0]
        aln_ref = str(best_aln[0])
        aln_mut = str(best_aln[1])
        
        matched_pairs = []
        ref_ptr = 0
        mut_ptr = 0
        
        for r_ch, m_ch in zip(aln_ref, aln_mut):
            r_has = (r_ch != '-')
            m_has = (m_ch != '-')
            
            if r_has and m_has:
                # Corresponding residue pair (could be match or substitution)
                matched_pairs.append((
                    ref_residues[ref_ptr],
                    mut_residues[mut_ptr],
                    ref_ptr,
                    mut_ptr
                ))
                
            if r_has:
                ref_ptr += 1
            if m_has:
                mut_ptr += 1
                
        return matched_pairs

    @staticmethod
    def superimpose_structures(
        ref_pdb_path: str,
        mut_pdb_path: str,
        output_superimposed_path: str,
        mutation_ref_pos: Optional[int] = None,
        local_window_size: int = 10
    ) -> Dict[str, Any]:
        """
        Performs structural superimposition using Python's Bio.PDB.Superimposer (Kabsch algorithm).
        1. Loads reference and mutant PDB structures.
        2. Sequence-based alignment of CA residues (handles deletions such as CFTR F508del seamlessly).
        3. Excludes unmatched residues and gaps.
        4. Calculates Kabsch optimal rotation & translation.
        5. Applies transformation to all atoms in mutant structure.
        6. Saves superimposed mutant PDB.
        7. Computes exact Global RMSD and Local RMSD using the formula:
           RMSD = sqrt[(1/N) * sum(||X_i - Y_i||^2)]
        8. Computes TM-score approximation and structural alignment coverage.
        """
        ref_struct = SuperimpositionService.load_structure(ref_pdb_path, "ref")
        mut_struct = SuperimpositionService.load_structure(mut_pdb_path, "mut")
        
        ref_ca_atoms, ref_residues, ref_seq = SuperimpositionService.extract_ca_atoms_and_sequence(ref_struct)
        mut_ca_atoms, mut_residues, mut_seq = SuperimpositionService.extract_ca_atoms_and_sequence(mut_struct)
        
        if not ref_ca_atoms:
            raise ValueError(f"No C-alpha atoms found in reference structure: {ref_pdb_path}")
        if not mut_ca_atoms:
            raise ValueError(f"No C-alpha atoms found in mutant structure: {mut_pdb_path}")
            
        # Map equivalent residues via sequence alignment
        matched_pairs = SuperimpositionService.map_residues_by_sequence_alignment(
            ref_residues=ref_residues,
            mut_residues=mut_residues,
            ref_seq=ref_seq,
            mut_seq=mut_seq
        )
        
        if len(matched_pairs) < 3:
            raise ValueError(f"Insufficient matched residues ({len(matched_pairs)}) for 3D superimposition (need >= 3).")
            
        # Extract matched CA atoms
        matched_ref_atoms = [pair[0]["CA"] for pair in matched_pairs]
        matched_mut_atoms = [pair[1]["CA"] for pair in matched_pairs]
        
        # Superimpose using Bio.PDB.Superimposer (Kabsch algorithm)
        superimposer = Superimposer()
        superimposer.set_atoms(matched_ref_atoms, matched_mut_atoms)
        
        # Apply transformation (rotation and translation) to the entire mutant structure
        superimposer.apply(mut_struct.get_atoms())
        
        # Save the superimposed mutant structure
        os.makedirs(os.path.dirname(os.path.abspath(output_superimposed_path)), exist_ok=True)
        io_pdb = PDBIO()
        io_pdb.set_structure(mut_struct)
        io_pdb.save(output_superimposed_path)
        
        # Compute exact Global C-alpha RMSD from coordinates
        ref_coords = np.array([atom.get_coord() for atom in matched_ref_atoms], dtype=np.float64)
        mut_coords = np.array([atom.get_coord() for atom in matched_mut_atoms], dtype=np.float64)
        
        diff = ref_coords - mut_coords
        sq_dist = np.sum(diff ** 2, axis=1) # ||X_i - Y_i||^2
        global_rmsd = math.sqrt(np.mean(sq_dist))
        
        # Coverage
        total_ref_res = len(ref_residues)
        total_mut_res = len(mut_residues)
        coverage = (len(matched_pairs) / max(total_ref_res, total_mut_res)) * 100.0
        
        # Local RMSD around mutation site
        local_rmsd = None
        local_window_residues = []
        target_idx = None
        
        # 1. Try finding explicit residue number matching mutation_ref_pos
        if mutation_ref_pos is not None:
            for idx, (r_res, m_res, r_idx, m_idx) in enumerate(matched_pairs):
                if r_res.get_id()[1] == mutation_ref_pos:
                    target_idx = idx
                    break

        # 2. If not found by residue number (e.g. domain slicing or None), auto-detect from sequence mismatch
        if target_idx is None:
            for idx, (r_res, m_res, r_idx, m_idx) in enumerate(matched_pairs):
                try:
                    r_code = seq1(r_res.get_resname())
                    m_code = seq1(m_res.get_resname())
                    if r_code != m_code:
                        target_idx = idx
                        break
                except Exception:
                    pass

        # 3. If still not found, check for insertion/deletion gap between ref_seq and mut_seq
        if target_idx is None:
            aligner = PairwiseAligner()
            aligner.mode = "global"
            aln_res = aligner.align(ref_seq, mut_seq)
            if aln_res:
                best_a = aln_res[0]
                ref_count = 0
                for r_ch, m_ch in zip(str(best_a[0]), str(best_a[1])):
                    if m_ch == '-' or r_ch == '-':
                        target_idx = min(ref_count, len(matched_pairs) - 1)
                        break
                    if r_ch != '-':
                        ref_count += 1

        # 4. If target_idx determined, calculate Local RMSD in a +/- local_window_size residue window
        if target_idx is not None and 0 <= target_idx < len(matched_pairs):
            target_res_num = matched_pairs[target_idx][0].get_id()[1]
            local_sq_dists = []
            for r_res, m_res, r_idx, m_idx in matched_pairs:
                res_num = r_res.get_id()[1]
                if abs(res_num - target_res_num) <= local_window_size:
                    c_ref = r_res["CA"].get_coord()
                    c_mut = m_res["CA"].get_coord()
                    d2 = np.sum((c_ref - c_mut) ** 2)
                    local_sq_dists.append(d2)
                    local_window_residues.append(res_num)
                    
            if local_sq_dists:
                local_rmsd = round(math.sqrt(float(np.mean(local_sq_dists))), 4)
                
        # TM-score approximation (Levitt-Gerstein formula)
        # TM = (1 / L_target) * sum(1 / (1 + (d_i / d0)^2))
        L_target = float(max(total_ref_res, total_mut_res))
        if L_target > 15:
            d0 = 1.24 * ((L_target - 15) ** (1.0 / 3.0)) - 1.8
            d0 = max(d0, 0.5) # Avoid zero or negative for small L
            distances = np.sqrt(sq_dist)
            tm_terms = 1.0 / (1.0 + (distances / d0) ** 2)
            tm_score = round(float(np.sum(tm_terms) / L_target), 4)
        else:
            tm_score = None
            
        return {
            "matched_residues_count": len(matched_pairs),
            "matched_ca_atoms_count": len(matched_pairs),
            "global_rmsd": round(global_rmsd, 4),
            "local_rmsd": local_rmsd,
            "local_window_residues": local_window_residues,
            "alignment_coverage": round(coverage, 2),
            "tm_score": tm_score,
            "superimposed_pdb_path": output_superimposed_path,
            "ref_pdb_path": ref_pdb_path,
            "mut_pdb_path": mut_pdb_path,
            "total_ref_residues": total_ref_res,
            "total_mut_residues": total_mut_res
        }
