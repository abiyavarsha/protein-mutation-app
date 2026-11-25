import numpy as np
from typing import Dict, Any, List, Optional
from Bio.PDB import PDBParser


# Amino acid biophysical properties
AA_PROPERTIES = {
    'A': {'name': 'Alanine', 'hydropathy': 1.8, 'mw': 89.1, 'charge': 0, 'volume': 88.6, 'class': 'Aliphatic hydrophobic'},
    'R': {'name': 'Arginine', 'hydropathy': -4.5, 'mw': 174.2, 'charge': +1, 'volume': 173.4, 'class': 'Basic positively charged'},
    'N': {'name': 'Asparagine', 'hydropathy': -3.5, 'mw': 132.1, 'charge': 0, 'volume': 114.1, 'class': 'Polar uncharged'},
    'D': {'name': 'Aspartate', 'hydropathy': -3.5, 'mw': 133.1, 'charge': -1, 'volume': 111.1, 'class': 'Acidic negatively charged'},
    'C': {'name': 'Cysteine', 'hydropathy': 2.5, 'mw': 121.2, 'charge': 0, 'volume': 108.5, 'class': 'Polar uncharged / Disulfide'},
    'E': {'name': 'Glutamate', 'hydropathy': -3.5, 'mw': 147.1, 'charge': -1, 'volume': 138.4, 'class': 'Acidic negatively charged'},
    'Q': {'name': 'Glutamine', 'hydropathy': -3.5, 'mw': 146.2, 'charge': 0, 'volume': 143.8, 'class': 'Polar uncharged'},
    'G': {'name': 'Glycine', 'hydropathy': -0.4, 'mw': 75.1, 'charge': 0, 'volume': 60.1, 'class': 'Conformationally flexible'},
    'H': {'name': 'Histidine', 'hydropathy': -3.2, 'mw': 155.2, 'charge': 0.1, 'volume': 153.2, 'class': 'Basic / Aromatic'},
    'I': {'name': 'Isoleucine', 'hydropathy': 4.5, 'mw': 131.2, 'charge': 0, 'volume': 166.7, 'class': 'Aliphatic hydrophobic'},
    'L': {'name': 'Leucine', 'hydropathy': 3.8, 'mw': 131.2, 'charge': 0, 'volume': 166.7, 'class': 'Aliphatic hydrophobic'},
    'K': {'name': 'Lysine', 'hydropathy': -3.9, 'mw': 146.2, 'charge': +1, 'volume': 168.6, 'class': 'Basic positively charged'},
    'M': {'name': 'Methionine', 'hydropathy': 1.9, 'mw': 149.2, 'charge': 0, 'volume': 162.9, 'class': 'Hydrophobic / Sulfur-containing'},
    'F': {'name': 'Phenylalanine', 'hydropathy': 2.8, 'mw': 165.2, 'charge': 0, 'volume': 189.9, 'class': 'Aromatic hydrophobic'},
    'P': {'name': 'Proline', 'hydropathy': -1.6, 'mw': 115.1, 'charge': 0, 'volume': 112.7, 'class': 'Conformationally rigid (imino)'},
    'S': {'name': 'Serine', 'hydropathy': -0.8, 'mw': 105.1, 'charge': 0, 'volume': 89.0, 'class': 'Polar uncharged'},
    'T': {'name': 'Threonine', 'hydropathy': -0.7, 'mw': 119.1, 'charge': 0, 'volume': 116.1, 'class': 'Polar uncharged'},
    'W': {'name': 'Tryptophan', 'hydropathy': -0.9, 'mw': 204.2, 'charge': 0, 'volume': 227.8, 'class': 'Aromatic hydrophobic'},
    'Y': {'name': 'Tyrosine', 'hydropathy': -1.3, 'mw': 181.2, 'charge': 0, 'volume': 193.6, 'class': 'Aromatic polar'},
    'V': {'name': 'Valine', 'hydropathy': 4.2, 'mw': 117.1, 'charge': 0, 'volume': 140.0, 'class': 'Aliphatic hydrophobic'}
}

try:
    from backend.dataset_registry import KNOWN_BIOLOGICAL_CONTEXT
except Exception:
    KNOWN_BIOLOGICAL_CONTEXT = {
        "HBB_E6V": (
            "Sickle cell anemia driver mutation: Glutamate (acidic, negatively charged, hydrophilic) "
            "is replaced by Valine (hydrophobic, aliphatic) at residue 6 of the mature beta-globin chain. "
            "This exposes a hydrophobic patch on the deoxygenated tetramer surface that interlocks with a complementary "
            "hydrophobic pocket (Phe85 and Leu88) on an adjacent beta-globin subunit, driving pathological polymerisation."
        ),
        "CFTR_F508del": (
            "Primary cause of Cystic Fibrosis (~70% of mutant alleles): In-frame deletion of Phenylalanine 508 "
            "in the first nucleotide-binding domain (NBD1). This prevents proper domain folding and disrupts the interface "
            "between NBD1 and intracellular loops, causing retention in the endoplasmic reticulum and degradation by ERAD."
        ),
        "CFTR_G551D": (
            "Cystic Fibrosis class III gating mutation: Glycine 551 at the highly conserved ABC transporter signature sequence "
            "(LSGGQ motif) in NBD1 is mutated to negatively charged Aspartate. The protein reaches the apical membrane but exhibits "
            "severely diminished open probability, treatable by the CFTR potentiator ivacaftor."
        ),
        "EGFR_L858R": (
            "Oncogenic driver in non-small cell lung cancer (NSCLC): Leucine 858 in the activation loop (A-loop) of the kinase domain "
            "is replaced by bulky, basic Arginine. This substitution destabilizes the auto-inhibited inactive conformation and stabilizes "
            "the active kinase state, imparting constitutive signaling and high sensitivity to EGFR tyrosine kinase inhibitors."
        ),
        "TP53_R175H": (
            "p53 conformational / structural hotspot mutation: Arg175 is located in the L2 loop near the zinc-coordination center. "
            "Substitution with Histidine disrupts zinc binding and local tertiary packing, causing partial unfolding of the core domain "
            "and loss of sequence-specific DNA-binding transcriptional capability."
        ),
        "TP53_R273H": (
            "p53 DNA-contact hotspot mutation: Arg273 is located in the S10-H2 loop and directly inserts into the DNA minor groove, "
            "forming electrostatic salt bridges with phosphate oxygens. Replacement with Histidine abolishes direct DNA contact without "
            "grossly perturbing overall core domain tertiary folding (evidenced by low global RMSD)."
        )
    }


class AnalysisService:

    @staticmethod
    def extract_plddt_scores(pdb_path: str) -> Dict[str, Any]:
        """
        Extracts pLDDT scores from the B-factor column of ESMFold / AlphaFold predicted PDB files.
        Returns mean pLDDT, per-residue pLDDT dict, and confidence categories.
        """
        parser = PDBParser(QUIET=True)
        try:
            struct = parser.get_structure("plddt", pdb_path)
        except Exception as e:
            return {"error": str(e), "mean_plddt": None, "plddt_by_res": {}}
            
        plddts = []
        plddt_by_res = {}
        for res in struct.get_residues():
            if res.has_id("CA") and res.get_id()[0] == " ":
                b = res["CA"].get_bfactor()
                res_num = res.get_id()[1]
                plddts.append(b)
                plddt_by_res[res_num] = round(float(b), 2)
                
        if not plddts:
            return {"mean_plddt": None, "plddt_by_res": {}}
            
        arr = np.array(plddts)
        very_high = float(np.sum(arr >= 90)) / len(arr) * 100.0
        confident = float(np.sum((arr >= 70) & (arr < 90))) / len(arr) * 100.0
        low = float(np.sum((arr >= 50) & (arr < 70))) / len(arr) * 100.0
        very_low = float(np.sum(arr < 50)) / len(arr) * 100.0
        
        return {
            "mean_plddt": round(float(np.mean(arr)), 2),
            "median_plddt": round(float(np.median(arr)), 2),
            "min_plddt": round(float(np.min(arr)), 2),
            "max_plddt": round(float(np.max(arr)), 2),
            "confidence_bands": {
                "very_high_plddt_gt90": round(very_high, 1),
                "confident_plddt_70_90": round(confident, 1),
                "low_confidence_plddt_50_70": round(low, 1),
                "very_low_confidence_plddt_lt50": round(very_low, 1)
            },
            "plddt_by_res": plddt_by_res
        }

    @staticmethod
    def calculate_biophysical_deltas(mutations: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Calculates changes in hydropathy, molecular weight, charge, and volume for each mutation.
        """
        results = []
        for m in mutations:
            m_type = m.get("mutation_type")
            ref_aa = m.get("ref_aa")
            mut_aa = m.get("mut_aa")
            notation = m.get("canonical_notation")
            mature_notation = m.get("mature_notation")
            
            if m_type == "SUBSTITUTION" and ref_aa in AA_PROPERTIES and mut_aa in AA_PROPERTIES:
                r_prop = AA_PROPERTIES[ref_aa]
                m_prop = AA_PROPERTIES[mut_aa]
                
                d_hydro = m_prop["hydropathy"] - r_prop["hydropathy"]
                d_mw = m_prop["mw"] - r_prop["mw"]
                d_charge = m_prop["charge"] - r_prop["charge"]
                d_vol = m_prop["volume"] - r_prop["volume"]
                
                results.append({
                    "notation": notation,
                    "mature_notation": mature_notation,
                    "type": "SUBSTITUTION",
                    "ref_residue": f"{r_prop['name']} ({ref_aa})",
                    "mut_residue": f"{m_prop['name']} ({mut_aa})",
                    "delta_hydropathy": round(d_hydro, 2),
                    "delta_molecular_weight": round(d_mw, 2),
                    "delta_charge": round(d_charge, 2),
                    "delta_volume": round(d_vol, 2),
                    "ref_class": r_prop["class"],
                    "mut_class": m_prop["class"],
                    "interpretation": AnalysisService._interpret_delta(d_hydro, d_charge, d_vol, r_prop, m_prop)
                })
            elif m_type == "DELETION":
                del_desc = []
                for char in (ref_aa or ""):
                    if char in AA_PROPERTIES:
                        del_desc.append(AA_PROPERTIES[char]["name"])
                results.append({
                    "notation": notation,
                    "mature_notation": mature_notation,
                    "type": "DELETION",
                    "ref_residue": f"{', '.join(del_desc)} ({ref_aa})",
                    "mut_residue": "None (Deleted)",
                    "delta_hydropathy": None,
                    "delta_molecular_weight": None,
                    "delta_charge": None,
                    "delta_volume": None,
                    "ref_class": AA_PROPERTIES.get(ref_aa[0], {}).get("class", "Unknown") if ref_aa else "Unknown",
                    "mut_class": "N/A",
                    "interpretation": f"Deletion of {ref_aa} creates a backbone gap, shortening the polypeptide chain and potentially altering local secondary structure and packing."
                })
            elif m_type == "INSERTION":
                results.append({
                    "notation": notation,
                    "mature_notation": mature_notation,
                    "type": "INSERTION",
                    "ref_residue": "None",
                    "mut_residue": mut_aa,
                    "delta_hydropathy": None,
                    "delta_molecular_weight": None,
                    "delta_charge": None,
                    "delta_volume": None,
                    "ref_class": "N/A",
                    "mut_class": "Inserted sequence",
                    "interpretation": f"Insertion of {mut_aa} increases chain length and may introduce new conformational flexibility or steric clashes."
                })
        return results

    @staticmethod
    def _interpret_delta(d_hydro: float, d_charge: float, d_vol: float, r_prop: dict, m_prop: dict) -> str:
        parts = []
        if abs(d_charge) > 0.5:
            if d_charge > 0:
                parts.append("Introduces positive electrostatic charge")
            else:
                parts.append("Introduces negative electrostatic charge / loss of positive charge")
        if abs(d_hydro) >= 2.0:
            if d_hydro > 0:
                parts.append("Substantial increase in hydrophobicity (may favor hydrophobic core burial or surface aggregation)")
            else:
                parts.append("Substantial decrease in hydrophobicity (polar/charged residue introduced)")
        if abs(d_vol) >= 30.0:
            if d_vol > 0:
                parts.append("Significant increase in steric volume (potential steric clash)")
            else:
                parts.append("Significant decrease in side-chain volume (potential loss of packing cavity)")
        if not parts:
            parts.append(f"Conservative biochemical change ({r_prop['class']} -> {m_prop['class']})")
        return "; ".join(parts)

    @staticmethod
    def generate_structural_report(
        protein_name: str,
        mutations: List[Dict[str, Any]],
        global_rmsd: float,
        local_rmsd: Optional[float],
        coverage: float,
        ref_plddt: Optional[Dict[str, Any]] = None,
        mut_plddt: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Synthesizes calculated structural metrics, biophysical property changes,
        and scientific interpretations while maintaining rigorous distinctions between
        calculated data and biological claims.
        """
        biophysical = AnalysisService.calculate_biophysical_deltas(mutations)
        
        # Check for known biological context
        bio_context = []
        for m in mutations:
            c_not = m.get("canonical_notation", "")
            m_not = m.get("mature_notation", "")
            for key, text in KNOWN_BIOLOGICAL_CONTEXT.items():
                k_prot, k_mut = key.split("_", 1)
                if k_prot.upper() == protein_name.upper():
                    if k_mut in c_not or k_mut in m_not:
                        bio_context.append(text)
                        
        # Formulate scientific discussion
        discussion_lines = [
            f"**Structural Superimposition Analysis for {protein_name}:**",
            f"- Structural superimposition yielded a global C-alpha RMSD of {global_rmsd:.4f} Å across {coverage:.1f}% alignment coverage."
        ]
        if local_rmsd is not None:
            discussion_lines.append(f"- Local C-alpha RMSD in the vicinity (±10 residues) of the mutation site was measured at {local_rmsd:.4f} Å.")
            
        if global_rmsd < 1.0:
            discussion_lines.append(
                "- The low global RMSD (< 1.0 Å) indicates that the overall tertiary fold of the mutant remains closely congruent with the reference backbone. "
                "Any functional alterations are likely driven by local side-chain electrostatic or steric effects rather than large-scale domain unfolding."
            )
        elif global_rmsd < 2.5:
            discussion_lines.append(
                "- A moderate global RMSD (1.0 - 2.5 Å) reflects detectable backbone conformational shifts, which may indicate altered domain orientation or flexible loop reorganizations."
            )
        else:
            discussion_lines.append(
                "- A significant global RMSD (> 2.5 Å) suggests substantial tertiary conformational rearrangement between reference and mutant predicted models."
            )
            
        if bio_context:
            discussion_lines.append("\n**Literature-Documented Biological Context:**")
            for ctx in bio_context:
                discussion_lines.append(f"- {ctx}")
                
        # Mandatory Scientific Limitations & Disclaimer
        disclaimer = (
            "IMPORTANT SCIENTIFIC LIMITATIONS:\n"
            "1. Calculated measurements (RMSD, sequence identity, coverage, and pLDDT) quantify geometric differences between predicted models.\n"
            "2. RMSD alone DOES NOT prove protein destabilization, altered folding kinetics, or loss of physiological function.\n"
            "3. Pathogenicity or clinical significance cannot be inferred solely from predicted computational structures.\n"
            "4. Free energy changes (ΔΔG) cannot be derived directly from RMSD; independent rigorous computational thermodynamics (such as FoldX, Rosetta, or FEP) or experimental biophysics (CD, ITC, DSC) are required."
        )
        
        return {
            "protein_name": protein_name,
            "global_rmsd": global_rmsd,
            "local_rmsd": local_rmsd,
            "alignment_coverage": coverage,
            "biophysical_changes": biophysical,
            "ref_plddt": ref_plddt,
            "mut_plddt": mut_plddt,
            "discussion": "\n".join(discussion_lines),
            "disclaimer": disclaimer
        }
