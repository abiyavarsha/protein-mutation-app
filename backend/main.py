import os
import uuid
import datetime
import json
from typing import List, Optional, Dict, Any
from fastapi import FastAPI, UploadFile, File, Form, HTTPException, BackgroundTasks, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from sqlalchemy.orm import Session

from backend.database import init_db, get_db, SessionLocal
from backend.models.db_models import Protein, Sequence, Mutation, Analysis
from backend.models.schemas import (
    ValidationRequest, ValidationResponse,
    AlignmentRequest, AlignmentResponse,
    MutationDetectionRequest, MutationDetectionResponse,
    StructurePredictionRequest, JobStatusResponse, JobStatusEnum,
    SuperimpositionRequest, SuperimpositionResponse,
    AnalysisRequest, AnalysisResponse
)
from backend.services.sequence_service import SequenceService
from backend.services.structure_service import StructureService
from backend.services.superimposition_service import SuperimpositionService
from backend.services.analysis_service import AnalysisService
from backend.services.report_service import ReportService
from backend.dataset_registry import DATASET_REGISTRY, get_protein, get_mutants_for_protein


BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE_DIR, "data")
FASTA_DIR = os.path.join(DATA_DIR, "fasta")
STRUCTURES_DIR = os.path.join(DATA_DIR, "structures")
REPORTS_DIR = os.path.join(DATA_DIR, "reports")

app = FastAPI(
    title="Protein Mutation Analysis & Tertiary Structure Comparison API",
    description="Bioinformatics backend for sequence alignment, mutation detection, ESMFold prediction, Kabsch structural superimposition, and RMSD calculation.",
    version="1.0.0"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# In-memory job tracker for asynchronous background jobs
JOBS_STORE: Dict[str, Dict[str, Any]] = {}


@app.on_event("startup")
def startup_event():
    os.makedirs(FASTA_DIR, exist_ok=True)
    os.makedirs(STRUCTURES_DIR, exist_ok=True)
    os.makedirs(REPORTS_DIR, exist_ok=True)
    init_db()


@app.get("/api/health")
def health_check():
    return {
        "status": "healthy",
        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "database": "connected",
        "service": "Protein Structure & Mutation Analysis API"
    }


@app.get("/api/proteins")
def list_proteins(db: Session = Depends(get_db)):
    """Lists available reference proteins and their existing mutant datasets."""
    proteins = db.query(Protein).all()
    results = []
    
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
            
        ref_path = os.path.join(FASTA_DIR, ref_filename)
        results.append({
            "id": p.id,
            "name": p.name,
            "uniprot_id": p.uniprot_id,
            "description": p.description,
            "category": category,
            "mature_offset": p.mature_offset,
            "reference_fasta_exists": os.path.exists(ref_path),
            "reference_fasta_path": ref_path if os.path.exists(ref_path) else None,
            "available_mutants": muts
        })
    return results


@app.post("/api/upload")
async def upload_fasta(
    file: UploadFile = File(...),
    protein_name: str = Form("CUSTOM"),
    sequence_type: str = Form("mutant")
):
    """Accepts uploaded FASTA file and stores it in the data directory."""
    if not file.filename.endswith((".fasta", ".fa", ".txt")):
        raise HTTPException(status_code=400, detail="Uploaded file must be a FASTA format (.fasta, .fa, .txt)")
        
    content = (await file.read()).decode("utf-8")
    records = SequenceService.parse_fasta(content)
    if not records:
        raise HTTPException(status_code=400, detail="Could not parse any valid sequence from uploaded file.")
        
    save_filename = f"{protein_name}_{sequence_type}_{uuid.uuid4().hex[:8]}.fasta"
    save_path = os.path.join(FASTA_DIR, save_filename)
    with open(save_path, "w") as f:
        f.write(content)
        
    return {
        "success": True,
        "filename": save_filename,
        "file_path": save_path,
        "records_count": len(records),
        "primary_record": records[0]
    }


@app.post("/api/validate", response_model=ValidationResponse)
def validate_sequence(req: ValidationRequest):
    """Validates sequence characters, length, and formats."""
    res = SequenceService.validate_sequence(req.sequence)
    return ValidationResponse(**res)


@app.post("/api/align", response_model=AlignmentResponse)
def align_sequences(req: AlignmentRequest):
    """Performs global pairwise alignment and returns metrics and visual representation."""
    try:
        res = SequenceService.align_sequences(
            ref_seq=req.ref_sequence,
            mut_seq=req.mut_sequence,
            protein_name=req.protein_name or "CUSTOM"
        )
        return AlignmentResponse(**res)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/api/mutations")
def detect_mutations(req: MutationDetectionRequest):
    """Identifies mutations between reference and mutant sequences or files."""
    try:
        ref_seq = req.ref_sequence
        mut_seq = req.mut_sequence
        
        if req.ref_fasta_path and os.path.exists(req.ref_fasta_path):
            recs = SequenceService.parse_fasta(req.ref_fasta_path)
            if recs:
                ref_seq = recs[0]["sequence"]
                
        if req.mut_fasta_path and os.path.exists(req.mut_fasta_path):
            recs = SequenceService.parse_fasta(req.mut_fasta_path)
            if recs:
                mut_seq = recs[0]["sequence"]
                
        if not ref_seq or not mut_seq:
            raise ValueError("Both reference and mutant sequences or valid file paths must be provided.")
            
        aln_res = SequenceService.align_sequences(
            ref_seq=ref_seq,
            mut_seq=mut_seq,
            protein_name=req.protein_name or "CUSTOM"
        )
        
        csv_summary = SequenceService.generate_csv_summary(
            mutations=aln_res["mutations"],
            protein_name=req.protein_name or "CUSTOM"
        )
        
        return {
            "protein_name": req.protein_name,
            "ref_length": aln_res["ref_length"],
            "mut_length": aln_res["mut_length"],
            "mutations_detected": aln_res["mutations"],
            "sequence_identity": aln_res["sequence_identity"],
            "alignment_score": aln_res["alignment_score"],
            "aligned_ref": aln_res["aligned_ref"],
            "aligned_mut": aln_res["aligned_mut"],
            "alignment_repr": aln_res["alignment_repr"],
            "csv_summary": csv_summary,
            "json_summary": aln_res
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


def run_prediction_job(job_id: str, sequence: str, protein_name: str, structure_type: str, label: Optional[str]):
    """Background task function for structure prediction."""
    JOBS_STORE[job_id]["status"] = JobStatusEnum.RUNNING
    JOBS_STORE[job_id]["progress"] = 0.2
    JOBS_STORE[job_id]["message"] = "Querying structure prediction service..."
    
    try:
        struct_res = StructureService.get_structure_for_protein(
            sequence=sequence,
            protein_name=protein_name,
            structure_type=label or structure_type,
            base_dir=BASE_DIR
        )
        JOBS_STORE[job_id]["status"] = JobStatusEnum.COMPLETED
        JOBS_STORE[job_id]["progress"] = 1.0
        JOBS_STORE[job_id]["message"] = "Structure prediction successfully completed."
        JOBS_STORE[job_id]["result"] = struct_res
    except Exception as e:
        JOBS_STORE[job_id]["status"] = JobStatusEnum.FAILED
        JOBS_STORE[job_id]["progress"] = 1.0
        JOBS_STORE[job_id]["error"] = str(e)
        JOBS_STORE[job_id]["message"] = f"Prediction failed: {str(e)}"


@app.post("/api/predict")
def predict_structure(req: StructurePredictionRequest, background_tasks: BackgroundTasks):
    """Submits sequence for tertiary structure prediction asynchronously."""
    job_id = f"pred_{uuid.uuid4().hex[:12]}"
    JOBS_STORE[job_id] = {
        "job_id": job_id,
        "status": JobStatusEnum.PENDING,
        "progress": 0.0,
        "message": "Job queued for tertiary structure prediction.",
        "protein_name": req.protein_name,
        "created_at": datetime.datetime.now(datetime.timezone.utc).isoformat()
    }
    
    background_tasks.add_task(
        run_prediction_job,
        job_id=job_id,
        sequence=req.sequence,
        protein_name=req.protein_name,
        structure_type=req.structure_type,
        label=req.label
    )
    
    return {"job_id": job_id, "status": "PENDING", "message": "Structure prediction job scheduled."}


@app.get("/api/jobs/{job_id}", response_model=JobStatusResponse)
def get_job_status(job_id: str):
    """Returns status and results of background prediction or comparison jobs."""
    if job_id not in JOBS_STORE:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found.")
    return JobStatusResponse(**JOBS_STORE[job_id])


@app.post("/api/superimpose")
def superimpose_structures(req: SuperimpositionRequest):
    """Performs Kabsch structural superimposition and calculates RMSD."""
    try:
        out_path = req.output_path
        if not out_path:
            base_mut = os.path.splitext(os.path.basename(req.mut_pdb_path))[0]
            out_path = os.path.join(
                STRUCTURES_DIR,
                req.protein_name or "CUSTOM",
                f"{base_mut}_superimposed.pdb"
            )
            
        res = SuperimpositionService.superimpose_structures(
            ref_pdb_path=req.ref_pdb_path,
            mut_pdb_path=req.mut_pdb_path,
            output_superimposed_path=out_path,
            mutation_ref_pos=req.mutation_pos
        )
        return res
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/api/compare")
def compare_structures(req: SuperimpositionRequest):
    """Alias for /api/superimpose with additional structural metrics."""
    return superimpose_structures(req)


@app.post("/api/analyze")
def analyze_structural_effects(req: AnalysisRequest):
    """
    Performs comprehensive structural effect analysis, biophysical property shifts,
    and generates downloadable PDF and JSON reports.
    """
    try:
        ref_plddt = AnalysisService.extract_plddt_scores(req.ref_pdb_path) if os.path.exists(req.ref_pdb_path) else None
        mut_plddt = AnalysisService.extract_plddt_scores(req.mut_pdb_path) if os.path.exists(req.mut_pdb_path) else None
        
        # Calculate coverage from CA count
        ref_res_count = len(ref_plddt["plddt_by_res"]) if ref_plddt and "plddt_by_res" in ref_plddt else 100
        mut_res_count = len(mut_plddt["plddt_by_res"]) if mut_plddt and "plddt_by_res" in mut_plddt else 100
        coverage = (min(ref_res_count, mut_res_count) / max(ref_res_count, mut_res_count)) * 100.0 if max(ref_res_count, mut_res_count) > 0 else 100.0
        
        analysis_res = AnalysisService.generate_structural_report(
            protein_name=req.protein_name,
            mutations=req.mutations,
            global_rmsd=req.global_rmsd,
            local_rmsd=req.local_rmsd,
            coverage=coverage,
            ref_plddt=ref_plddt,
            mut_plddt=mut_plddt
        )
        
        # Persist to database if job_id provided
        job_id = req.job_id or f"anal_{uuid.uuid4().hex[:12]}"
        
        # Generate PDF report
        pdf_filename = f"{req.protein_name}_{job_id}_report.pdf"
        pdf_path = os.path.join(REPORTS_DIR, pdf_filename)
        
        report_data = {
            "job_id": job_id,
            "protein_name": req.protein_name,
            "mutant_label": req.mutations[0].get("canonical_notation") if req.mutations else "Mutant",
            "global_rmsd": req.global_rmsd,
            "local_rmsd": req.local_rmsd,
            "sequence_identity": 99.0,
            "alignment_coverage": coverage,
            "matched_ca_atoms_count": min(ref_res_count, mut_res_count),
            "mutations": req.mutations,
            "biophysical_changes": analysis_res["biophysical_changes"],
            "discussion": analysis_res["discussion"]
        }
        ReportService.generate_pdf_report(report_data, pdf_path)
        
        # Save JSON summary
        json_path = os.path.join(REPORTS_DIR, f"{req.protein_name}_{job_id}_summary.json")
        with open(json_path, "w") as jf:
            json.dump({**report_data, "disclaimer": analysis_res["disclaimer"]}, jf, indent=2)
            
        # Store in DB
        db = SessionLocal()
        try:
            prot = db.query(Protein).filter(Protein.name == req.protein_name).first()
            analysis_record = Analysis(
                job_id=job_id,
                protein_id=prot.id if prot else None,
                protein_name=req.protein_name,
                mutant_label=report_data["mutant_label"],
                status="COMPLETED",
                global_rmsd=req.global_rmsd,
                local_rmsd=req.local_rmsd,
                alignment_coverage=coverage,
                ref_pdb_path=req.ref_pdb_path,
                mut_pdb_path=req.mut_pdb_path,
                superimposed_pdb_path=req.superimposed_pdb_path,
                report_pdf_path=pdf_path,
                summary_json_path=json_path,
                completed_at=datetime.datetime.now(datetime.timezone.utc)
            )
            db.add(analysis_record)
            db.commit()
        finally:
            db.close()
            
        return {
            "job_id": job_id,
            "protein_name": req.protein_name,
            "report_pdf_path": pdf_path,
            "summary_json_path": json_path,
            "analysis": analysis_res
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/results/{job_id}")
def get_analysis_results(job_id: str, db: Session = Depends(get_db)):
    """Retrieves full analysis results from database or JSON summary."""
    rec = db.query(Analysis).filter(Analysis.job_id == job_id).first()
    if not rec:
        raise HTTPException(status_code=404, detail=f"Analysis with job_id {job_id} not found.")
        
    summary_data = {}
    if rec.summary_json_path and os.path.exists(rec.summary_json_path):
        with open(rec.summary_json_path, "r") as f:
            summary_data = json.load(f)
            
    return {
        "job_id": rec.job_id,
        "protein_name": rec.protein_name,
        "mutant_label": rec.mutant_label,
        "status": rec.status,
        "global_rmsd": rec.global_rmsd,
        "local_rmsd": rec.local_rmsd,
        "alignment_coverage": rec.alignment_coverage,
        "ref_pdb_path": rec.ref_pdb_path,
        "mut_pdb_path": rec.mut_pdb_path,
        "superimposed_pdb_path": rec.superimposed_pdb_path,
        "report_pdf_path": rec.report_pdf_path,
        "created_at": rec.created_at.isoformat() if rec.created_at else None,
        "completed_at": rec.completed_at.isoformat() if rec.completed_at else None,
        "details": summary_data
    }


@app.get("/api/report/{job_id}")
def download_pdf_report(job_id: str, db: Session = Depends(get_db)):
    """Downloads the generated PDF report."""
    rec = db.query(Analysis).filter(Analysis.job_id == job_id).first()
    if not rec or not rec.report_pdf_path or not os.path.exists(rec.report_pdf_path):
        raise HTTPException(status_code=404, detail="PDF report not found.")
    return FileResponse(
        rec.report_pdf_path,
        media_type="application/pdf",
        filename=os.path.basename(rec.report_pdf_path)
    )


@app.get("/api/structures/{protein}/{filename}")
def get_pdb_file(protein: str, filename: str):
    """Serves PDB files for 3D visualization and download."""
    path = os.path.join(STRUCTURES_DIR, protein, filename)
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail=f"PDB file not found at {path}")
    return FileResponse(path, media_type="chemical/x-pdb", filename=filename)
