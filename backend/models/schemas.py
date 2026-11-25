from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from enum import Enum


class ProteinNameEnum(str, Enum):
    HBB = "HBB"
    CFTR = "CFTR"
    EGFR = "EGFR"
    TP53 = "TP53"
    CUSTOM = "CUSTOM"


class MutationTypeEnum(str, Enum):
    SUBSTITUTION = "SUBSTITUTION"
    DELETION = "DELETION"
    INSERTION = "INSERTION"


class JobStatusEnum(str, Enum):
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"


class ValidationRequest(BaseModel):
    sequence: str
    protein_name: Optional[str] = None
    header: Optional[str] = None


class ValidationResponse(BaseModel):
    is_valid: bool
    length: int
    cleaned_sequence: str
    errors: List[str] = []
    warnings: List[str] = []


class AlignmentRequest(BaseModel):
    ref_sequence: str
    mut_sequence: str
    protein_name: Optional[str] = "CUSTOM"
    ref_header: Optional[str] = "Reference"
    mut_header: Optional[str] = "Mutant"


class MutationItem(BaseModel):
    mutation_type: MutationTypeEnum
    ref_pos: int
    mut_pos: int
    ref_aa: Optional[str] = None
    mut_aa: Optional[str] = None
    canonical_notation: str
    mature_pos: Optional[int] = None
    mature_notation: Optional[str] = None
    description: str


class AlignmentResponse(BaseModel):
    ref_length: int
    mut_length: int
    alignment_score: float
    sequence_identity: float
    similarity_percentage: float
    gaps_count: int
    aligned_ref: str
    aligned_mut: str
    alignment_repr: str
    mutations: List[MutationItem] = []


class MutationDetectionRequest(BaseModel):
    ref_fasta_path: Optional[str] = None
    mut_fasta_path: Optional[str] = None
    ref_sequence: Optional[str] = None
    mut_sequence: Optional[str] = None
    protein_name: Optional[str] = "CUSTOM"


class MutationDetectionResponse(BaseModel):
    protein_name: str
    ref_length: int
    mut_length: int
    mutations_detected: List[MutationItem]
    sequence_identity: float
    alignment_score: float
    csv_summary: str
    json_summary: Dict[str, Any]
    label_consistency: Optional[str] = None


class StructurePredictionRequest(BaseModel):
    sequence: str
    protein_name: str
    structure_type: str = "mutant" # "reference" or "mutant"
    label: Optional[str] = None


class JobStatusResponse(BaseModel):
    job_id: str
    status: JobStatusEnum
    progress: float = 0.0
    message: str = ""
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


class SuperimpositionRequest(BaseModel):
    ref_pdb_path: str
    mut_pdb_path: str
    protein_name: Optional[str] = "CUSTOM"
    output_path: Optional[str] = None
    mutation_pos: Optional[int] = None


class SuperimpositionResponse(BaseModel):
    matched_residues_count: int
    matched_ca_atoms_count: int
    global_rmsd: float
    local_rmsd: Optional[float] = None
    local_window_residues: Optional[List[int]] = None
    alignment_coverage: float
    tm_score: Optional[float] = None
    superimposed_pdb_path: str
    ref_pdb_path: str
    mut_pdb_path: str


class AnalysisRequest(BaseModel):
    job_id: Optional[str] = None
    protein_name: str
    ref_sequence: str
    mut_sequence: str
    ref_pdb_path: str
    mut_pdb_path: str
    superimposed_pdb_path: str
    mutations: List[Dict[str, Any]]
    global_rmsd: float
    local_rmsd: Optional[float] = None


class AnalysisResponse(BaseModel):
    protein_name: str
    mutation_summary: List[Dict[str, Any]]
    global_rmsd: float
    local_rmsd: Optional[float] = None
    sequence_identity: float
    alignment_coverage: float
    ref_plddt_mean: Optional[float] = None
    mut_plddt_mean: Optional[float] = None
    local_plddt_mean: Optional[float] = None
    biophysical_changes: List[Dict[str, Any]]
    structural_interpretation: str
    disclaimer: str
