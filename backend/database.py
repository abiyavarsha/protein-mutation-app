import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from backend.models.db_models import Base, Protein

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE_DIR, "data")
DB_PATH = os.path.join(DATA_DIR, "app.db")
SQLALCHEMY_DATABASE_URL = f"sqlite:///{DB_PATH}"

os.makedirs(DATA_DIR, exist_ok=True)

engine = create_engine(
    SQLALCHEMY_DATABASE_URL,
    connect_args={"check_same_thread": False}
)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def init_db():
    Base.metadata.create_all(bind=engine)
    
    # Seed known proteins if not already present
    db = SessionLocal()
    try:
        from backend.dataset_registry import get_all_proteins
        seed_proteins = [
            {
                "name": p["name"],
                "uniprot_id": p["uniprot_id"],
                "description": p["description"],
                "mature_offset": p.get("mature_offset", 0)
            }
            for p in get_all_proteins()
        ]
    except Exception:
        # Fallback to original 4 benchmark proteins if registry cannot be loaded
        seed_proteins = [
            {
                "name": "HBB",
                "uniprot_id": "P68871",
                "description": "Hemoglobin subunit beta (Homo sapiens)",
                "mature_offset": -1
            },
            {
                "name": "CFTR",
                "uniprot_id": "P13569",
                "description": "Cystic fibrosis transmembrane conductance regulator (Homo sapiens)",
                "mature_offset": 0
            },
            {
                "name": "EGFR",
                "uniprot_id": "P00533",
                "description": "Epidermal growth factor receptor (Homo sapiens)",
                "mature_offset": 0
            },
            {
                "name": "TP53",
                "uniprot_id": "P04637",
                "description": "Cellular tumor antigen p53 (Homo sapiens)",
                "mature_offset": 0
            }
        ]
        
    try:
        for p_data in seed_proteins:
            existing = db.query(Protein).filter(Protein.name == p_data["name"]).first()
            if not existing:
                protein = Protein(
                    name=p_data["name"],
                    uniprot_id=p_data["uniprot_id"],
                    description=p_data["description"],
                    mature_offset=p_data["mature_offset"]
                )
                db.add(protein)
        db.commit()
    finally:
        db.close()
