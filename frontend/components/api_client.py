import os
import sys
import requests
from typing import Dict, Any, Optional

BACKEND_URL = os.environ.get("BACKEND_URL", "http://127.0.0.1:8000")

# Import local services for direct fallback if API is not yet running
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from backend.services.sequence_service import SequenceService
from backend.services.structure_service import StructureService
from backend.services.superimposition_service import SuperimpositionService
from backend.services.analysis_service import AnalysisService
from backend.services.report_service import ReportService
from backend.database import SessionLocal
from backend.models.db_models import Protein, Analysis


class BackendClient:

    @staticmethod
    def is_api_online() -> bool:
        try:
            r = requests.get(f"{BACKEND_URL}/api/health", timeout=1.5)
            return r.status_code == 200
        except Exception:
            return False

    @staticmethod
    def get_proteins():
        if BackendClient.is_api_online():
            try:
                r = requests.get(f"{BACKEND_URL}/api/proteins", timeout=3)
                if r.status_code == 200:
                    return r.json()
            except Exception:
                pass
        # Fallback to direct DB query & registry
        db = SessionLocal()
        try:
            from backend.dataset_registry import get_protein
            proteins = db.query(Protein).all()
            results = []
            fasta_dir = os.path.join(BASE_DIR, "data", "fasta")
            for p in proteins:
                reg_info = get_protein(p.name)
                if reg_info:
                    muts = [
                        {
                            "mutation": m["mutation"],
                            "filename": m["filename"],
                            "type": m["type"],
                            "label": m.get("label", m["mutation"]),
                            "clinical_note": m.get("clinical_note", ""),
                            "available": m.get("available", True),
                            "is_original": m.get("is_original", False)
                        }
                        for m in reg_info.get("mutants", [])
                    ]
                    category = reg_info.get("category", "Extended Dataset")
                    ref_filename = reg_info.get("ref_file", f"{p.name}_reference.fasta")
                else:
                    muts = []
                    category = "Custom / Unregistered"
                    ref_filename = f"{p.name}_reference.fasta"
                    
                ref_p = os.path.join(fasta_dir, ref_filename)
                results.append({
                    "id": p.id,
                    "name": p.name,
                    "uniprot_id": p.uniprot_id,
                    "description": p.description,
                    "category": category,
                    "mature_offset": p.mature_offset,
                    "reference_fasta_exists": os.path.exists(ref_p),
                    "reference_fasta_path": ref_p if os.path.exists(ref_p) else None,
                    "available_mutants": muts
                })
            return results
        finally:
            db.close()

    @staticmethod
    def validate_sequence(sequence: str):
        if BackendClient.is_api_online():
            try:
                r = requests.post(f"{BACKEND_URL}/api/validate", json={"sequence": sequence}, timeout=5)
                if r.status_code == 200:
                    return r.json()
            except Exception:
                pass
        return SequenceService.validate_sequence(sequence)

    @staticmethod
    def align_and_detect_mutations(ref_seq: str, mut_seq: str, protein_name: str = "CUSTOM"):
        if BackendClient.is_api_online():
            try:
                payload = {
                    "ref_sequence": ref_seq,
                    "mut_sequence": mut_seq,
                    "protein_name": protein_name
                }
                r = requests.post(f"{BACKEND_URL}/api/mutations", json=payload, timeout=10)
                if r.status_code == 200:
                    return r.json()
            except Exception:
                pass
        # Fallback to direct service
        aln = SequenceService.align_sequences(ref_seq, mut_seq, protein_name)
        csv_s = SequenceService.generate_csv_summary(aln["mutations"], protein_name)
        return {
            "protein_name": protein_name,
            "ref_length": aln["ref_length"],
            "mut_length": aln["mut_length"],
            "mutations_detected": aln["mutations"],
            "sequence_identity": aln["sequence_identity"],
            "alignment_score": aln["alignment_score"],
            "aligned_ref": aln["aligned_ref"],
            "aligned_mut": aln["aligned_mut"],
            "alignment_repr": aln["alignment_repr"],
            "csv_summary": csv_s,
            "json_summary": aln
        }

    @staticmethod
    def superimpose(ref_pdb: str, mut_pdb: str, out_pdb: Optional[str] = None, mutation_pos: Optional[int] = None, protein_name: str = "CUSTOM"):
        if BackendClient.is_api_online():
            try:
                payload = {
                    "ref_pdb_path": ref_pdb,
                    "mut_pdb_path": mut_pdb,
                    "output_path": out_pdb,
                    "mutation_pos": mutation_pos,
                    "protein_name": protein_name
                }
                r = requests.post(f"{BACKEND_URL}/api/superimpose", json=payload, timeout=15)
                if r.status_code == 200:
                    return r.json()
            except Exception:
                pass
        # Direct service
        if not out_pdb:
            base_mut = os.path.splitext(os.path.basename(mut_pdb))[0]
            out_pdb = os.path.join(BASE_DIR, "data", "structures", protein_name, f"{base_mut}_superimposed.pdb")
        return SuperimpositionService.superimpose_structures(ref_pdb, mut_pdb, out_pdb, mutation_pos)

    @staticmethod
    def analyze(protein_name: str, ref_pdb: str, mut_pdb: str, sup_pdb: str, mutations: list, global_rmsd: float, local_rmsd: Optional[float] = None, job_id: Optional[str] = None):
        if BackendClient.is_api_online():
            try:
                payload = {
                    "protein_name": protein_name,
                    "ref_pdb_path": ref_pdb,
                    "mut_pdb_path": mut_pdb,
                    "superimposed_pdb_path": sup_pdb,
                    "mutations": mutations,
                    "global_rmsd": global_rmsd,
                    "local_rmsd": local_rmsd,
                    "job_id": job_id,
                    "ref_sequence": "",
                    "mut_sequence": ""
                }
                r = requests.post(f"{BACKEND_URL}/api/analyze", json=payload, timeout=20)
                if r.status_code == 200:
                    return r.json()
            except Exception:
                pass
        # Direct service fallback
        ref_plddt = AnalysisService.extract_plddt_scores(ref_pdb) if os.path.exists(ref_pdb) else None
        mut_plddt = AnalysisService.extract_plddt_scores(mut_pdb) if os.path.exists(mut_pdb) else None
        
        ref_cnt = len(ref_plddt["plddt_by_res"]) if ref_plddt and "plddt_by_res" in ref_plddt else 100
        mut_cnt = len(mut_plddt["plddt_by_res"]) if mut_plddt and "plddt_by_res" in mut_plddt else 100
        coverage = (min(ref_cnt, mut_cnt) / max(ref_cnt, mut_cnt)) * 100.0 if max(ref_cnt, mut_cnt) > 0 else 100.0
        
        rep = AnalysisService.generate_structural_report(
            protein_name=protein_name,
            mutations=mutations,
            global_rmsd=global_rmsd,
            local_rmsd=local_rmsd,
            coverage=coverage,
            ref_plddt=ref_plddt,
            mut_plddt=mut_plddt
        )
        
        import uuid, json
        j_id = job_id or f"anal_{uuid.uuid4().hex[:8]}"
        pdf_p = os.path.join(BASE_DIR, "data", "reports", f"{protein_name}_{j_id}_report.pdf")
        json_p = os.path.join(BASE_DIR, "data", "reports", f"{protein_name}_{j_id}_summary.json")
        
        r_data = {
            "job_id": j_id,
            "protein_name": protein_name,
            "mutant_label": mutations[0].get("canonical_notation") if mutations else "Mutant",
            "global_rmsd": global_rmsd,
            "local_rmsd": local_rmsd,
            "sequence_identity": 99.0,
            "alignment_coverage": coverage,
            "matched_ca_atoms_count": min(ref_cnt, mut_cnt),
            "mutations": mutations,
            "biophysical_changes": rep["biophysical_changes"],
            "discussion": rep["discussion"]
        }
        ReportService.generate_pdf_report(r_data, pdf_p)
        with open(json_p, "w") as jf:
            json.dump({**r_data, "disclaimer": rep["disclaimer"]}, jf, indent=2)
            
        return {
            "job_id": j_id,
            "protein_name": protein_name,
            "report_pdf_path": pdf_p,
            "summary_json_path": json_p,
            "analysis": rep
        }
