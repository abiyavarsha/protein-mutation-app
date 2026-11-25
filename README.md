# Protein Mutation Analysis & Tertiary Structure Comparison Application

A complete, self-contained Python bioinformatics web application with a **Streamlit frontend** and **FastAPI backend** for:
- Reference and mutant protein FASTA sequence validation and pairwise alignment
- Mutation calling: substitutions, in-frame deletions (e.g., CFTR F508del), and insertions
- Mature-protein chain numbering (e.g. human beta-globin cleaved Met1; E6V in mature chain vs E7V in precursor)
- Automated 3D tertiary structure prediction via hosted ESMFold
- Structural superimposition via Python using the **Kabsch algorithm** (`Bio.PDB.Superimposer` and NumPy)
- Exact C-alpha RMSD quantification (Global RMSD and Local RMSD around mutation sites)
- Interactive 3D WebGL protein visualization (`3Dmol.js`) with pLDDT confidence coloring and mutation highlighting
- Structural effect and biophysical alteration analysis (hydropathy, charge, volume, molecular weight)
- Publication-quality PDF report generation (`ReportLab`)
- Persistent SQLite metadata database (`SQLAlchemy`)

---

## 📦 What Is Included in This Package

```
protein_app/
│
├── data/
│   ├── fasta/                         # All Reference and Mutant FASTA datasets (77 validated files)
│   │   ├── HBB_reference.fasta        # UniProt P68871 (147 aa)
│   │   ├── CFTR_reference.fasta       # UniProt P13569 (1480 aa)
│   │   ├── EGFR_reference.fasta       # UniProt P00533 (1210 aa)
│   │   ├── TP53_reference.fasta       # UniProt P04637 (393 aa)
│   │   ├── [24 Extended Reference FASTA files: KRAS, BRAF, PIK3CA, SOD1, SNCA, APP, ...]
│   │   ├── HBB_E6V_mutant.fasta       # Sickle cell driver mutation
│   │   ├── CFTR_F508del_mutant.fasta  # In-frame deletion mutation
│   │   ├── CFTR_G551D_mutant.fasta    # Class III gating mutation
│   │   ├── EGFR_L858R_mutant.fasta    # NSCLC oncogenic kinase driver
│   │   ├── TP53_R175H_mutant.fasta    # Conformational zinc-binding hotspot
│   │   ├── TP53_R273H_mutant.fasta    # Direct DNA minor-groove contact hotspot
│   │   └── [44 Extended Mutant FASTA files across 24 human disease genes]
│   ├── structures/                    # Pre-generated and validated PDB structures
│   │   ├── HBB/ (reference.pdb, E6V_mutant.pdb, E6V_superimposed.pdb)
│   │   ├── CFTR/ (reference.pdb, F508del_mutant.pdb, F508del_superimposed.pdb, G551D_mutant.pdb, G551D_superimposed.pdb)
│   │   ├── EGFR/ (reference.pdb, L858R_mutant.pdb, L858R_superimposed.pdb)
│   │   └── TP53/ (reference.pdb, R175H_mutant.pdb, R175H_superimposed.pdb, R273H_mutant.pdb, R273H_superimposed.pdb)
│   ├── reports/                       # Generated PDF and JSON reports
│   └── app.db                         # SQLite database for analysis history
│
├── backend/
│   ├── dataset_registry.py            # Centralized registry for 28 genes, mature offsets, and 50+ mutations
│   ├── models/
│   │   ├── db_models.py               # SQLAlchemy models (Protein, Sequence, Mutation, Analysis)
│   │   └── schemas.py                 # Pydantic request & response schemas
│   ├── services/
│   │   ├── sequence_service.py        # Biopython PairwiseAligner & mutation caller
│   │   ├── structure_service.py       # ESMFold client with validation & local cache
│   │   ├── superimposition_service.py # Kabsch algorithm & sequence-based CA matching
│   │   ├── analysis_service.py        # Biophysical deltas & pLDDT stats
│   │   └── report_service.py          # ReportLab PDF report generator
│   ├── database.py                    # SQLite engine & session management
│   └── main.py                        # FastAPI endpoints & background worker
│
├── frontend/
│   ├── app.py                         # 10-page Streamlit web application with tiered dataset selector
│   └── components/
│       ├── api_client.py              # Resilient API communication client
│       └── visualizer_3d.py           # 3Dmol.js WebGL viewer with pLDDT coloring
│
├── tests/
│   ├── test_sequence.py               # Unit tests: FASTA validation, alignments, mature numbering
│   ├── test_structure.py              # Unit tests: PDB validation, Kabsch alignment, RMSD formula
│   ├── test_api.py                    # Unit tests: FastAPI route verification
│   └── test_extended_datasets.py      # Unit tests: 28-protein registry coverage, mutation calling, RYR1 note
│
├── requirements.txt                   # Complete dependencies
├── run.bat                            # Windows Command Prompt 1-click launcher
├── run.ps1                            # Windows PowerShell 1-click launcher
├── start_servers.py                   # Python cross-platform launcher (Win/Mac/Linux)
├── verify_pipeline.py                 # Automated pipeline verification script
└── README.md                          # Documentation and instructions
```

---

## 🚀 Quick Start Guide (How to Run on Another Laptop)

### Prerequisites
1. **Python 3.10+** (Python 3.10, 3.11, 3.12, 3.13, 3.14).
   - Download from [python.org](https://www.python.org/downloads/).
   - ⚠️ **Important during installation on Windows:** Check the box **"Add python.exe to PATH"**.
2. Any standard laptop with **8 GB RAM** or higher.

---

### Method 1: One-Click Startup (Easiest on Windows)

1. Unzip `protein_mutation_app.zip` to any folder on your laptop.
2. Open the unzipped `protein_app` folder.
3. Double-click **`run.bat`** (or right-click **`run.ps1`** -> *Run with PowerShell*).
4. The script will automatically:
   - Check and install required packages from `requirements.txt`.
   - Start the FastAPI backend on `http://127.0.0.1:8000`.
   - Start the Streamlit frontend on `http://localhost:8501`.
   - Open the web application directly in your browser!

---

### Method 2: Universal Cross-Platform Startup (Windows / macOS / Linux)

Open a terminal or command prompt in the `protein_app` directory and run:

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Run both servers simultaneously
python start_servers.py
```

Then open your browser to:
👉 **[http://localhost:8501](http://localhost:8501)**

---

### Method 3: Manual Startup in Two Separate Terminals

#### Terminal 1: Start FastAPI Backend
```bash
cd protein_app
uvicorn backend.main:app --host 127.0.0.1 --port 8000
```
*(Backend Swagger API docs will be live at `http://127.0.0.1:8000/docs`)*

#### Terminal 2: Start Streamlit Frontend
```bash
cd protein_app
streamlit run frontend/app.py --server.port 8501
```
*(Streamlit frontend will open automatically at `http://localhost:8501`)*

---

## 🧪 Running Automated Unit Tests

To verify that all modules (FASTA parsing, mutation calling, mature numbering, Kabsch superimposition, RMSD mathematical formula, and API endpoints) are functioning properly on the new laptop:

```bash
cd protein_app
python -m unittest discover -s tests -p "test_*.py"
```

All 26 unit tests will run and report:
```
Ran 26 tests in 0.56s — OK
```

---

## 📊 End-to-End Pipeline Verification Script

You can also run the pre-built verification script which processes all six benchmark datasets, calculates RMSDs, and outputs PDF reports:

```bash
cd protein_app
python verify_pipeline.py
```

Expected benchmark metrics:
- **HBB E6V:** Global RMSD = `0.0111 Å` | Coverage = `100.0%` | TM-Score = `1.0000`
- **CFTR F508del:** Global RMSD = `0.3593 Å` | Coverage = `99.5%` | TM-Score = `0.9909`
- **CFTR G551D:** Global RMSD = `0.1252 Å` | Coverage = `100.0%` | TM-Score = `0.9994`
- **EGFR L858R:** Global RMSD = `0.2432 Å` | Coverage = `100.0%` | TM-Score = `0.9943`
- **TP53 R175H:** Global RMSD = `0.7624 Å` | Coverage = `100.0%` | TM-Score = `0.9892`
- **TP53 R273H:** Global RMSD = `0.8151 Å` | Coverage = `100.0%` | TM-Score = `0.9878`

---

## 🔬 Scientific Methodology & Disclaimers

1. **Structural Superimposition**: Implemented via Biopython's `Bio.PDB.Superimposer` using the Kabsch algorithm. Mismatched residues, gaps, and deletions (e.g. CFTR F508del) are excluded via sequence alignment before coordinate transformation.
2. **RMSD Formula**:
   $$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^N \|X_i - Y_i\|^2}$$
3. **Mandatory Scientific Limitations**:
   - RMSD quantifies geometric differences between models; it does not prove protein destabilization or loss of physiological function.
   - Clinical pathogenicity cannot be inferred solely from in silico predicted structures.
   - Thermodynamic stability ($\Delta\Delta G$) cannot be derived directly from RMSD; independent computational methods (e.g. FoldX, Rosetta, or FEP) or experimental biophysics are required.
