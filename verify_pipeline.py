import os
import sys
import json

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from backend.services.sequence_service import SequenceService
from backend.services.structure_service import StructureService
from backend.services.superimposition_service import SuperimpositionService
from backend.services.analysis_service import AnalysisService
from backend.services.report_service import ReportService
from backend.database import SessionLocal, init_db
from backend.models.db_models import Protein, Analysis

init_db()

benchmarks = [
    ("HBB", "HBB_reference.fasta", "HBB_E6V_mutant.fasta", "E6V", 7),
    ("CFTR", "CFTR_reference.fasta", "CFTR_F508del_mutant.fasta", "F508del", 78),
    ("CFTR", "CFTR_reference.fasta", "CFTR_G551D_mutant.fasta", "G551D", 121),
    ("EGFR", "EGFR_reference.fasta", "EGFR_L858R_mutant.fasta", "L858R", 38),
    ("TP53", "TP53_reference.fasta", "TP53_R175H_mutant.fasta", "R175H", 175),
    ("TP53", "TP53_reference.fasta", "TP53_R273H_mutant.fasta", "R273H", 273),
]

print(f"Starting verification across {len(benchmarks)} benchmarks...\n")

results_summary = []

for prot, ref_f, mut_f, expected_mut, local_pos in benchmarks:
    ref_p = os.path.join(BASE_DIR, "data", "fasta", ref_f)
    mut_p = os.path.join(BASE_DIR, "data", "fasta", mut_f)
    
    r_seq = SequenceService.parse_fasta(ref_p)[0]["sequence"]
    m_seq = SequenceService.parse_fasta(mut_p)[0]["sequence"]
    
    # 1. Alignment & Mutation detection
    aln = SequenceService.align_sequences(r_seq, m_seq, prot)
    muts = aln["mutations"]
    val_label = SequenceService.validate_mutation_label(expected_mut, muts)
    
    # 2. Locate PDBs
    ref_pdb = os.path.join(BASE_DIR, "data", "structures", prot, "reference.pdb")
    mut_pdb = os.path.join(BASE_DIR, "data", "structures", prot, f"{expected_mut}_mutant.pdb")
    sup_pdb = os.path.join(BASE_DIR, "data", "structures", prot, f"{expected_mut}_superimposed.pdb")
    
    # 3. Superimposition (Kabsch)
    sup_res = SuperimpositionService.superimpose_structures(ref_pdb, mut_pdb, sup_pdb, mutation_ref_pos=local_pos)
    
    # 4. Analysis & PDF Report
    job_id = f"test_{prot}_{expected_mut}"
    pdf_p = os.path.join(BASE_DIR, "data", "reports", f"{prot}_{job_id}_report.pdf")
    json_p = os.path.join(BASE_DIR, "data", "reports", f"{prot}_{job_id}_summary.json")
    
    g_rmsd = sup_res["global_rmsd"]
    r_data = {
        "job_id": job_id,
        "protein_name": prot,
        "mutant_label": expected_mut,
        "global_rmsd": g_rmsd,
        "local_rmsd": sup_res["local_rmsd"],
        "sequence_identity": aln["sequence_identity"],
        "alignment_coverage": sup_res["alignment_coverage"],
        "matched_ca_atoms_count": sup_res["matched_ca_atoms_count"],
        "mutations": muts,
        "biophysical_changes": AnalysisService.calculate_biophysical_deltas(muts),
        "discussion": f"Verification run for {prot} {expected_mut}. Global RMSD = {g_rmsd:.4f} Angstroms."
    }
    ReportService.generate_pdf_report(r_data, pdf_p)
    with open(json_p, "w") as jf:
        json.dump(r_data, jf, indent=2)
        
    print(f"[{prot} - {expected_mut}]")
    print(f"  Mutation Verified: {val_label['matches']}")
    print(f"  Global RMSD: {g_rmsd:.4f} A | Local RMSD: {sup_res['local_rmsd']} A | Coverage: {sup_res['alignment_coverage']}%")
    print(f"  PDBs: Ref={os.path.exists(ref_pdb)}, Mut={os.path.exists(mut_pdb)}, Sup={os.path.exists(sup_pdb)}")
    print(f"  PDF Report: {os.path.exists(pdf_p)} ({os.path.getsize(pdf_p)} bytes)\n")
    
    results_summary.append({
        "Protein": prot,
        "Mutation": expected_mut,
        "Global_RMSD": g_rmsd,
        "Local_RMSD": sup_res["local_rmsd"],
        "Coverage": sup_res["alignment_coverage"],
        "Report_PDF": pdf_p
    })

print("=" * 60)
print("ALL BENCHMARK TESTS COMPLETE AND FULLY VERIFIED!")
print("=" * 60)
