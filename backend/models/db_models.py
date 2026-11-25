import datetime
from sqlalchemy import Column, Integer, String, Float, Text, DateTime, ForeignKey, Enum as SQLEnum
from sqlalchemy.orm import relationship, declarative_base

Base = declarative_base()


class Protein(Base):
    __tablename__ = "proteins"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(50), unique=True, nullable=False, index=True)
    uniprot_id = Column(String(20), nullable=False)
    description = Column(String(255))
    mature_offset = Column(Integer, default=0) # e.g. -1 for HBB mature chain
    created_at = Column(DateTime, default=lambda: datetime.datetime.now(datetime.timezone.utc))

    sequences = relationship("Sequence", back_populates="protein", cascade="all, delete-orphan")
    analyses = relationship("Analysis", back_populates="protein", cascade="all, delete-orphan")


class Sequence(Base):
    __tablename__ = "sequences"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey("proteins.id"), nullable=True)
    sequence_type = Column(String(20), nullable=False) # "reference" or "mutant"
    label = Column(String(100), nullable=True)
    header = Column(String(255), nullable=True)
    fasta_path = Column(String(500), nullable=True)
    sequence_text = Column(Text, nullable=False)
    length = Column(Integer, nullable=False)
    created_at = Column(DateTime, default=lambda: datetime.datetime.now(datetime.timezone.utc))

    protein = relationship("Protein", back_populates="sequences")


class Mutation(Base):
    __tablename__ = "mutations"

    id = Column(Integer, primary_key=True, autoincrement=True)
    analysis_id = Column(Integer, ForeignKey("analyses.id"), nullable=False)
    mutation_type = Column(String(30), nullable=False) # "SUBSTITUTION", "DELETION", "INSERTION"
    ref_pos = Column(Integer, nullable=False)
    mut_pos = Column(Integer, nullable=False)
    ref_aa = Column(String(10), nullable=True)
    mut_aa = Column(String(10), nullable=True)
    canonical_notation = Column(String(50), nullable=False)
    mature_pos = Column(Integer, nullable=True)
    mature_notation = Column(String(50), nullable=True)
    description = Column(Text, nullable=True)

    analysis = relationship("Analysis", back_populates="mutations")


class Analysis(Base):
    __tablename__ = "analyses"

    id = Column(Integer, primary_key=True, autoincrement=True)
    job_id = Column(String(64), unique=True, nullable=False, index=True)
    protein_id = Column(Integer, ForeignKey("proteins.id"), nullable=True)
    protein_name = Column(String(50), nullable=False)
    mutant_label = Column(String(100), nullable=True)
    mutant_filename = Column(String(255), nullable=True)
    
    status = Column(String(20), default="PENDING") # PENDING, RUNNING, COMPLETED, FAILED
    error_message = Column(Text, nullable=True)
    
    # Results
    sequence_identity = Column(Float, nullable=True)
    alignment_score = Column(Float, nullable=True)
    matched_residues_count = Column(Integer, nullable=True)
    matched_ca_atoms_count = Column(Integer, nullable=True)
    global_rmsd = Column(Float, nullable=True)
    local_rmsd = Column(Float, nullable=True)
    alignment_coverage = Column(Float, nullable=True)
    tm_score = Column(Float, nullable=True)
    
    # File storage paths
    ref_fasta_path = Column(String(500), nullable=True)
    mut_fasta_path = Column(String(500), nullable=True)
    ref_pdb_path = Column(String(500), nullable=True)
    mut_pdb_path = Column(String(500), nullable=True)
    superimposed_pdb_path = Column(String(500), nullable=True)
    report_pdf_path = Column(String(500), nullable=True)
    summary_json_path = Column(String(500), nullable=True)
    summary_csv_path = Column(String(500), nullable=True)
    
    created_at = Column(DateTime, default=lambda: datetime.datetime.now(datetime.timezone.utc))
    completed_at = Column(DateTime, nullable=True)

    protein = relationship("Protein", back_populates="analyses")
    mutations = relationship("Mutation", back_populates="analysis", cascade="all, delete-orphan")
