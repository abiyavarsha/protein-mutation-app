import os
import sys
import json
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# Ensure backend modules can be imported
APP_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if APP_DIR not in sys.path:
    sys.path.append(APP_DIR)

from frontend.components.api_client import BackendClient
from frontend.components.visualizer_3d import render_3d_viewer
from backend.services.sequence_service import SequenceService
from backend.services.structure_service import StructureService
from backend.services.superimposition_service import SuperimpositionService
from backend.services.analysis_service import AnalysisService
from backend.services.report_service import ReportService
from backend.database import SessionLocal, init_db
from backend.models.db_models import Protein, Analysis
from backend.dataset_registry import DATASET_REGISTRY, get_dataset_map, get_protein


def get_mutant_pdb_paths(struct_dir, selected_mutant_file):
    """
    Resolves the exact reference, mutant, and superimposed PDB paths for the selected mutant dataset.
    Prevents cross-contamination when multiple mutants exist for the same protein.
    """
    ref_pdb = os.path.join(struct_dir, "reference.pdb")
    mut_pdb = None
    sup_pdb = None
    
    if not os.path.exists(struct_dir):
        return ref_pdb, None, None
        
    mut_tag = None
    if selected_mutant_file:
        parts = os.path.basename(selected_mutant_file).replace(".fasta", "").split("_")
        if len(parts) >= 2:
            mut_tag = parts[1]
            
    if mut_tag:
        p1 = os.path.join(struct_dir, f"{mut_tag}_mutant.pdb")
        s1 = os.path.join(struct_dir, f"{mut_tag}_mutant_superimposed.pdb")
        s2 = os.path.join(struct_dir, f"{mut_tag}_superimposed.pdb")
        if os.path.exists(p1):
            mut_pdb = p1
            sup_pdb = s1 if os.path.exists(s1) else (s2 if os.path.exists(s2) else None)
            
    if not mut_pdb:
        for f in sorted(os.listdir(struct_dir)):
            if f.endswith(".pdb") and "mutant" in f and "superimposed" not in f:
                if mut_tag and mut_tag.lower() in f.lower():
                    mut_pdb = os.path.join(struct_dir, f)
                    break
                elif not mut_pdb:
                    mut_pdb = os.path.join(struct_dir, f)
                    
    if mut_pdb and not sup_pdb:
        base = os.path.splitext(os.path.basename(mut_pdb))[0]
        c1 = os.path.join(struct_dir, f"{base}_superimposed.pdb")
        base_clean = base.replace("_mutant", "")
        c2 = os.path.join(struct_dir, f"{base_clean}_superimposed.pdb")
        if os.path.exists(c1):
            sup_pdb = c1
        elif os.path.exists(c2):
            sup_pdb = c2
            
    return ref_pdb, mut_pdb, sup_pdb

# Initialize DB
init_db()

# Page config
st.set_page_config(
    page_title="Protein Mutation & Tertiary Structure Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Unique Scientific Styling
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;500;600&display=swap');

    html, body, [class*="css"] {
        font-family: 'Plus Jakarta Sans', system-ui, -apple-system, sans-serif;
    }

    /* Unique Brand Hero Banner */
    .brand-banner {
        background: linear-gradient(135deg, #0B132B 0%, #1C2541 60%, #1E1B4B 100%);
        border: 1px solid rgba(56, 189, 248, 0.25);
        border-radius: 16px;
        padding: 24px 28px;
        margin-bottom: 22px;
        box-shadow: 0 10px 30px -5px rgba(11, 19, 43, 0.35);
    }
    .brand-title {
        font-size: 2.15rem;
        font-weight: 800;
        letter-spacing: -0.025em;
        background: linear-gradient(135deg, #38BDF8 0%, #818CF8 50%, #C084FC 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 6px;
    }
    .brand-subtitle {
        color: #94A3B8;
        font-size: 0.96rem;
        margin-bottom: 14px;
        line-height: 1.5;
    }
    .tag-chip {
        display: inline-block;
        padding: 3px 10px;
        border-radius: 20px;
        font-size: 0.74rem;
        font-weight: 600;
        margin-right: 6px;
        background: rgba(99, 102, 241, 0.18);
        color: #A5B4FC;
        border: 1px solid rgba(129, 140, 248, 0.3);
    }

    /* Interactive Pipeline Stepper Navigation */
    .stepper-bar {
        display: flex;
        align-items: center;
        justify-content: space-between;
        background: #F8FAFC;
        border: 1px solid #E2E8F0;
        border-radius: 12px;
        padding: 8px 14px;
        margin-bottom: 24px;
        overflow-x: auto;
    }
    .stepper-step {
        font-size: 0.8rem;
        font-weight: 600;
        color: #64748B;
        display: flex;
        align-items: center;
        gap: 6px;
        white-space: nowrap;
    }
    .stepper-step.step-active {
        color: #4F46E5;
        font-weight: 700;
    }
    .stepper-step.step-done {
        color: #059669;
    }
    .step-num {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        width: 22px;
        height: 22px;
        border-radius: 50%;
        background: #E2E8F0;
        color: #475569;
        font-size: 0.7rem;
    }
    .step-active .step-num {
        background: linear-gradient(135deg, #4F46E5, #06B6D4);
        color: white;
        box-shadow: 0 0 8px rgba(79, 70, 229, 0.4);
    }
    .step-done .step-num {
        background: #10B981;
        color: white;
    }
    .stepper-arrow {
        color: #CBD5E1;
        font-size: 0.72rem;
        margin: 0 4px;
    }

    /* Distinctive Glassmorphism Metric Cards */
    .metric-card-pro {
        background: linear-gradient(135deg, #FFFFFF, #F8FAFC);
        border: 1px solid #E2E8F0;
        border-top: 3px solid #6366F1;
        border-radius: 12px;
        padding: 16px;
        text-align: center;
        box-shadow: 0 4px 14px -2px rgba(15, 23, 42, 0.05);
        transition: transform 0.2s ease, box-shadow 0.2s ease;
    }
    .metric-card-pro:hover {
        transform: translateY(-2px);
        box-shadow: 0 8px 22px -3px rgba(99, 102, 241, 0.15);
    }
    .metric-value-pro {
        font-size: 1.85rem;
        font-weight: 800;
        color: #1E293B;
        line-height: 1.2;
    }
    .metric-label-pro {
        font-size: 0.75rem;
        color: #64748B;
        text-transform: uppercase;
        font-weight: 700;
        letter-spacing: 0.07em;
        margin-top: 4px;
    }

    /* Active Dataset Badge */
    .dataset-info-box {
        background: #F8FAFC;
        border: 1px solid #CBD5E1;
        border-radius: 10px;
        padding: 10px 14px;
        margin: 8px 0 14px 0;
    }

    /* Disclaimer alert */
    .disclaimer-box {
        background: #FEF2F2;
        border-left: 4px solid #EF4444;
        padding: 14px 18px;
        border-radius: 6px;
        color: #991B1B;
        font-size: 0.88rem;
        margin-top: 18px;
    }

    /* Styled buttons */
    .stButton>button {
        border-radius: 8px;
        font-weight: 600;
        transition: all 0.2s ease;
    }
</style>
""", unsafe_allow_html=True)


# Session state initialization
if "selected_protein" not in st.session_state:
    st.session_state.selected_protein = "HBB"
if "selected_mutant_file" not in st.session_state:
    st.session_state.selected_mutant_file = "HBB_E6V_mutant.fasta"
if "current_analysis" not in st.session_state:
    st.session_state.current_analysis = None
if "ref_seq" not in st.session_state:
    st.session_state.ref_seq = ""
if "mut_seq" not in st.session_state:
    st.session_state.mut_seq = ""
if "ref_pdb_path" not in st.session_state:
    st.session_state.ref_pdb_path = None
if "mut_pdb_path" not in st.session_state:
    st.session_state.mut_pdb_path = None
if "sup_pdb_path" not in st.session_state:
    st.session_state.sup_pdb_path = None
if "current_page_idx" not in st.session_state:
    st.session_state.current_page_idx = 0


# Streamlined 8-Step Bioinformatics Pipeline (Combined 5 & 6; Removed 10)
PAGES = [
    "1. Overview & Benchmark Hub",
    "2. FASTA Sequence Input",
    "3. Alignment & Mutation Calling",
    "4. Tertiary Structure Prediction",
    "5. Superimposition & RMSD Analysis",
    "6. Interactive 3D Visualization",
    "7. Biophysical Impact & Dashboard",
    "8. Publication PDF Report"
]


# Sidebar Header
st.sidebar.markdown("""
<div style="display: flex; align-items: center; gap: 12px; margin-bottom: 6px;">
    <div style="background: linear-gradient(135deg, #4F46E5, #06B6D4); border-radius: 12px; width: 44px; height: 44px; display: flex; align-items: center; justify-content: center; font-size: 24px; box-shadow: 0 4px 10px rgba(79, 70, 229, 0.3);">
        🧬
    </div>
    <div>
        <div style="font-size: 1.15rem; font-weight: 800; line-height: 1.2;">Protein Structure Suite</div>
        <div style="font-size: 0.74rem; color: #94A3B8;">Biopython · ESMFold · Kabsch · 3Dmol</div>
    </div>
</div>
""", unsafe_allow_html=True)

st.sidebar.markdown("---")

# User-friendly workflow navigation
def on_nav_change():
    if st.session_state.nav_radio in PAGES:
        st.session_state.current_page_idx = PAGES.index(st.session_state.nav_radio)

page = st.sidebar.radio(
    "Workflow Navigation:",
    PAGES,
    index=st.session_state.current_page_idx,
    key="nav_radio",
    on_change=on_nav_change
)
# Always derive current_idx from session_state so Prev/Next buttons drive the page
current_idx = st.session_state.current_page_idx
# Keep page in sync with current_page_idx (handles Prev/Next button navigation)
page = PAGES[current_idx]


# Pipeline Stepper and Footer Helper Functions
def render_pipeline_stepper(active_idx: int):
    short_titles = [
        "1. Overview", "2. FASTA", "3. Alignment", "4. Prediction",
        "5. Kabsch & RMSD", "6. 3D Viewer", "7. Biophysics", "8. Report"
    ]
    parts = []
    for idx, name in enumerate(short_titles):
        if idx == active_idx:
            cls = "stepper-step step-active"
            badge = f'<span class="step-num">{idx+1}</span>'
        elif idx < active_idx:
            cls = "stepper-step step-done"
            badge = '<span class="step-num">✓</span>'
        else:
            cls = "stepper-step"
            badge = f'<span class="step-num">{idx+1}</span>'
        parts.append(f'<span class="{cls}">{badge} {name.split(". ")[1]}</span>')
        if idx < len(short_titles) - 1:
            parts.append('<span class="stepper-arrow">➔</span>')
    st.markdown(f'<div class="stepper-bar">{"".join(parts)}</div>', unsafe_allow_html=True)


def render_step_footer(active_idx: int):
    st.markdown("<div style='margin-top: 35px;'></div>", unsafe_allow_html=True)
    c_prev, c_center, c_next = st.columns([2, 3, 2])
    with c_prev:
        if active_idx > 0:
            prev_label = PAGES[active_idx - 1].split(". ")[1]
            if st.button(f"⬅ Previous: {prev_label}", key=f"btn_prev_{active_idx}", use_container_width=True):
                st.session_state.current_page_idx = active_idx - 1
                st.rerun()
    with c_center:
        st.caption(f"<div style='text-align: center; color: #64748B; font-weight: 600; padding-top: 8px;'>Stage {active_idx + 1} of {len(PAGES)}</div>", unsafe_allow_html=True)
    with c_next:
        if active_idx < len(PAGES) - 1:
            next_label = PAGES[active_idx + 1].split(". ")[1]
            if st.button(f"Next: {next_label} ➔", key=f"btn_next_{active_idx}", type="primary", use_container_width=True):
                st.session_state.current_page_idx = active_idx + 1
                st.rerun()


st.sidebar.markdown("---")
st.sidebar.subheader("Select Dataset")

# Centralized dataset registry loading
DATASET_MAP = get_dataset_map()
DATASET_MAP["Custom Upload"] = {
    "uniprot": "Custom",
    "description": "User-uploaded Reference & Mutant sequences",
    "category": "Custom",
    "ref_file": None,
    "mutants": [],
    "has_mutants": False
}

benchmark_keys = ["HBB", "CFTR", "EGFR", "TP53"]
extended_keys = [k for k in DATASET_MAP.keys() if k not in benchmark_keys and k != "Custom Upload"]

# Category filter
cat_options = [
    "Original Benchmark Proteins (4)",
    "Extended Disease Datasets (24)",
    "Custom Upload"
]

# Preserve category matching selected_protein
if st.session_state.selected_protein in benchmark_keys:
    default_cat_idx = 0
elif st.session_state.selected_protein == "Custom Upload":
    default_cat_idx = 2
else:
    default_cat_idx = 1

category_choice = st.sidebar.radio("Dataset Tier / Group:", cat_options, index=default_cat_idx)

if category_choice.startswith("Original"):
    available_proteins = benchmark_keys
elif category_choice.startswith("Extended"):
    available_proteins = extended_keys
else:
    available_proteins = ["Custom Upload"]

default_p_idx = available_proteins.index(st.session_state.selected_protein) if st.session_state.selected_protein in available_proteins else 0
protein_choice = st.sidebar.selectbox("Target Protein:", available_proteins, index=default_p_idx)
st.session_state.selected_protein = protein_choice

if protein_choice != "Custom Upload":
    prot_info = DATASET_MAP[protein_choice]
    mut_opts = prot_info.get("mutants", [])
    
    if mut_opts:
        mut_labels = [m[0] for m in mut_opts]
        default_m_idx = 0
        curr_file = st.session_state.selected_mutant_file
        for idx, m in enumerate(mut_opts):
            if m[1] == curr_file:
                default_m_idx = idx
                break
        mut_choice = st.sidebar.selectbox("Mutant Dataset:", mut_labels, index=default_m_idx)
        selected_file = [m[1] for m in mut_opts if m[0] == mut_choice][0]
        st.session_state.selected_mutant_file = selected_file
        st.sidebar.caption(f"Reference: `{prot_info['ref_file']}`")
        st.sidebar.caption(f"Mutant: `{selected_file}`")
        st.sidebar.caption(f"UniProt: `{prot_info['uniprot']}`")
    else:
        st.session_state.selected_mutant_file = None
        st.sidebar.caption(f"Reference: `{prot_info['ref_file']}`")
        st.sidebar.caption(f"UniProt: `{prot_info['uniprot']}`")
        if protein_choice == "RYR1":
            st.sidebar.warning(
                "⚠️ Notice on RYR1: In canonical UniProt sequence P21817, residue 615 is Serine (S) rather than Arginine (R) "
                "(manifest note: neighbours A612 V613 R614 S615 N616 Q617 D618). As instructed, no assumed sequence was fabricated."
            )
        else:
            st.sidebar.warning("No mutant datasets registered for this entry.")

# One-Click Full Pipeline Runner
st.sidebar.markdown("---")
run_all_btn = st.sidebar.button("⚡ Run Full Pipeline Now", type="primary", use_container_width=True)


def load_dataset_sequences(prot_name, mut_file):
    fasta_dir = os.path.join(APP_DIR, "data", "fasta")
    ref_file = DATASET_MAP.get(prot_name, {}).get("ref_file")
    ref_p = os.path.join(fasta_dir, ref_file) if ref_file else ""
    mut_p = os.path.join(fasta_dir, mut_file) if mut_file else ""
    
    r_recs = SequenceService.parse_fasta(ref_p) if (ref_p and os.path.exists(ref_p)) else []
    m_recs = SequenceService.parse_fasta(mut_p) if (mut_p and os.path.exists(mut_p)) else []
    
    r_seq = r_recs[0]["sequence"] if r_recs else ""
    m_seq = m_recs[0]["sequence"] if m_recs else ""
    return r_seq, m_seq, ref_p, mut_p


# Auto-load sequences when dataset changes
if protein_choice != "Custom Upload":
    r_s, m_s, r_p, m_p = load_dataset_sequences(protein_choice, st.session_state.selected_mutant_file)
    st.session_state.ref_seq = r_s
    st.session_state.mut_seq = m_s
    st.session_state.ref_fasta_path = r_p
    st.session_state.mut_fasta_path = m_p


# Helper function to execute complete pipeline
def execute_full_pipeline():
    with st.spinner("Executing end-to-end bioinformatics pipeline..."):
        p_name = st.session_state.selected_protein
        r_seq = st.session_state.ref_seq
        m_seq = st.session_state.mut_seq
        
        if not r_seq or not m_seq:
            st.error("Both reference and mutant sequences must be loaded before running the pipeline.")
            return
            
        # 1. Alignment & Mutations
        aln = SequenceService.align_sequences(r_seq, m_seq, p_name)
        muts = aln["mutations"]
        
        # 2. Locate / Verify Tertiary Structures
        struct_dir = os.path.join(APP_DIR, "data", "structures", p_name)
        os.makedirs(struct_dir, exist_ok=True)
        
        # Resolve specific mutant file tag
        mut_tag = "mutant"
        if st.session_state.selected_mutant_file:
            parts = os.path.basename(st.session_state.selected_mutant_file).replace(".fasta", "").split("_")
            if len(parts) >= 2:
                mut_tag = parts[1]
        elif muts:
            mut_tag = muts[0].get('mature_notation') or muts[0].get('canonical_notation') or "mutant"
            
        ref_pdb = os.path.join(struct_dir, "reference.pdb")
        mut_pdb = os.path.join(struct_dir, f"{mut_tag}_mutant.pdb")
        sup_pdb = os.path.join(struct_dir, f"{mut_tag}_mutant_superimposed.pdb")
        
        # Check existing structures
        if not os.path.exists(ref_pdb) or not os.path.exists(mut_pdb):
            # Windowing for long sequences (> 400 aa) to comply with hosted ESMFold limits
            if len(r_seq) > 400 or len(m_seq) > 400:
                mut_pos = muts[0].get("ref_pos", 200) if muts else 200
                half_win = 100
                start_r = max(0, mut_pos - half_win)
                end_r = min(len(r_seq), start_r + 200)
                start_m = max(0, mut_pos - half_win)
                end_m = min(len(m_seq), start_m + 200)
                r_fold = r_seq[start_r:end_r]
                m_fold = m_seq[start_m:end_m]
                superimpose_mut_pos = mut_pos - start_r
            else:
                r_fold = r_seq
                m_fold = m_seq
                superimpose_mut_pos = muts[0].get("ref_pos") if muts else None
                
            try:
                if not os.path.exists(ref_pdb):
                    StructureService.predict_with_esmfold(r_fold, ref_pdb)
                if not os.path.exists(mut_pdb):
                    StructureService.predict_with_esmfold(m_fold, mut_pdb)
            except Exception as pred_err:
                st.error(f"❌ Tertiary structure prediction service encountered an issue: {pred_err}")
                st.info(
                    "💡 **Troubleshooting & Guidance:**\n"
                    "- The hosted ESMFold public API (`api.esmatlas.com`) can encounter rate limits or temporary `HTTP 504: Gateway Timeout` under heavy community traffic.\n"
                    "- Pre-computed, validated tertiary structures are ready-to-run for **HBB (E6V)**, **CFTR (F508del, G551D)**, **EGFR (L858R)**, and **TP53 (R175H, R273H)**.\n"
                    "- You may retry prediction in a few moments, or upload your own reference and mutant `.pdb` files directly on Page 4 (Tertiary Structure Prediction)."
                )
                return
        else:
            # Handle mature offset if applicable
            superimpose_mut_pos = muts[0].get("ref_pos") if muts else None
            
        # 3. Superimposition (Kabsch)
        sup_res = SuperimpositionService.superimpose_structures(
            ref_pdb_path=ref_pdb,
            mut_pdb_path=mut_pdb,
            output_superimposed_path=sup_pdb,
            mutation_ref_pos=superimpose_mut_pos
        )
        
        # 4. Analysis & PDF Report
        import uuid
        job_id = f"job_{uuid.uuid4().hex[:8]}"
        analysis_res = BackendClient.analyze(
            protein_name=p_name,
            ref_pdb=ref_pdb,
            mut_pdb=mut_pdb,
            sup_pdb=sup_pdb,
            mutations=muts,
            global_rmsd=sup_res["global_rmsd"],
            local_rmsd=sup_res["local_rmsd"],
            job_id=job_id
        )
        
        st.session_state.current_analysis = {
            "job_id": job_id,
            "alignment": aln,
            "superimposition": sup_res,
            "analysis": analysis_res,
            "ref_pdb_path": ref_pdb,
            "mut_pdb_path": mut_pdb,
            "sup_pdb_path": sup_pdb
        }
        st.session_state.ref_pdb_path = ref_pdb
        st.session_state.mut_pdb_path = mut_pdb
        st.session_state.sup_pdb_path = sup_pdb
        
        st.success(f"Pipeline finished! Global RMSD: {sup_res['global_rmsd']:.4f} Å | Coverage: {sup_res['alignment_coverage']:.1f}%")

if run_all_btn:
    execute_full_pipeline()


# ==========================================
# PAGE 1: HOME / PROJECT OVERVIEW
# ==========================================
if page.startswith("1."):
    render_pipeline_stepper(0)
    st.markdown('''
<div class="brand-banner">
  <div class="brand-title">🧬 ProteinMutSuite</div>
  <div class="brand-subtitle">Automated protein mutation analysis — sequence alignment, ESMFold structure prediction, Kabsch superimposition, RMSD quantification &amp; 3D visualization.</div>
  <span class="tag-chip">ESMFold</span>
  <span class="tag-chip">Kabsch Algorithm</span>
  <span class="tag-chip">28 Proteins</span>
  <span class="tag-chip">50+ Mutations</span>
  <span class="tag-chip">3Dmol.js</span>
</div>
''', unsafe_allow_html=True)
    col1, col2, col3, col4, col5 = st.columns(5)
    with col1:
        st.markdown('<div class="metric-card-pro"><div class="metric-value-pro">4</div><div class="metric-label-pro">Benchmark Proteins</div></div>', unsafe_allow_html=True)
    with col2:
        st.markdown('<div class="metric-card-pro"><div class="metric-value-pro">24</div><div class="metric-label-pro">Extended Proteins</div></div>', unsafe_allow_html=True)
    with col3:
        st.markdown('<div class="metric-card-pro"><div class="metric-value-pro">50+</div><div class="metric-label-pro">Mutant Datasets</div></div>', unsafe_allow_html=True)
    with col4:
        st.markdown('<div class="metric-card-pro"><div class="metric-value-pro">Kabsch</div><div class="metric-label-pro">Superimposition Engine</div></div>', unsafe_allow_html=True)
    with col5:
        st.markdown('<div class="metric-card-pro"><div class="metric-value-pro">ESMFold</div><div class="metric-label-pro">Tertiary Predictor</div></div>', unsafe_allow_html=True)
        
    st.markdown("### System Architecture & Pipeline Workflow")
    
    st.markdown("""
    ```mermaid
    flowchart TD
        A["FASTA Input (Reference & Mutant)"] --> B["Sequence Validation (IUPAC Standards)"]
        B --> C["Pairwise Sequence Alignment (Biopython PairwiseAligner)"]
        C --> D["Mutation Detection (Substitutions, In-frame Deletions, Insertions)"]
        D --> E["Tertiary Structure Prediction (ESMFold Hosted API)"]
        E --> F["Sequence-Based CA Atom Mapping"]
        F --> G["Optimal Structural Superimposition (Kabsch Algorithm / Bio.PDB)"]
        G --> H["RMSD Calculation (Global RMSD & Local RMSD Around Mutation)"]
        H --> I["Interactive 3D WebGL Visualization (3Dmol.js)"]
        H --> J["Biophysical & Structural Effect Analysis"]
        J --> K["Downloadable Publication PDF Report (ReportLab)"]
    ```
    """)
    
    st.markdown("### 1. Existing Validated Benchmark Datasets (Tiers 1-4)")
    benchmark_df = pd.DataFrame([
        {"Protein": "HBB (Haemoglobin beta)", "UniProt": "P68871", "Mutation": "E6V (Mature) / E7V", "Type": "Substitution", "FASTA File": "HBB_E6V_mutant.fasta", "Clinical Relevance": "Sickle Cell Anemia"},
        {"Protein": "CFTR", "UniProt": "P13569", "Mutation": "F508del", "Type": "In-frame Deletion", "FASTA File": "CFTR_F508del_mutant.fasta", "Clinical Relevance": "Cystic Fibrosis (~70% alleles)"},
        {"Protein": "CFTR", "UniProt": "P13569", "Mutation": "G551D", "Type": "Substitution", "FASTA File": "CFTR_G551D_mutant.fasta", "Clinical Relevance": "Cystic Fibrosis Gating Defect"},
        {"Protein": "EGFR", "UniProt": "P00533", "Mutation": "L858R", "Type": "Substitution", "FASTA File": "EGFR_L858R_mutant.fasta", "Clinical Relevance": "NSCLC Oncogenic Driver (TKI-sensitive)"},
        {"Protein": "TP53 (p53)", "UniProt": "P04637", "Mutation": "R175H", "Type": "Substitution", "FASTA File": "TP53_R175H_mutant.fasta", "Clinical Relevance": "Conformational Zinc-loss Hotspot"},
        {"Protein": "TP53 (p53)", "UniProt": "P04637", "Mutation": "R273H", "Type": "Substitution", "FASTA File": "TP53_R273H_mutant.fasta", "Clinical Relevance": "Direct DNA-contact Hotspot"}
    ])
    st.dataframe(benchmark_df, use_container_width=True, hide_index=True)

    with st.expander("🔬 2. View Extended Disease Datasets (24 Additional Proteins, 44+ Mutations)"):
        extended_rows = []
        for g, info in DATASET_REGISTRY.items():
            if info.get("category") == "Original Benchmark" and g in ["HBB", "CFTR", "EGFR", "TP53"]:
                # Also list additional mutants if any
                add_muts = [m for m in info.get("mutants", []) if not m.get("is_original", False)]
                for m in add_muts:
                    extended_rows.append({
                        "Gene": g,
                        "UniProt": info["uniprot_id"],
                        "Mutation": m["mutation"],
                        "Type": m["type"],
                        "Filename": m["filename"],
                        "Clinical Annotation": m.get("clinical_note", "")
                    })
            else:
                for m in info.get("mutants", []):
                    extended_rows.append({
                        "Gene": g,
                        "UniProt": info["uniprot_id"],
                        "Mutation": m["mutation"],
                        "Type": m["type"],
                        "Filename": m["filename"],
                        "Clinical Annotation": m.get("clinical_note", "")
                    })
        if extended_rows:
            st.dataframe(pd.DataFrame(extended_rows), use_container_width=True, hide_index=True)

    st.markdown('<div class="disclaimer-box"><b>Important Scientific Disclaimer:</b> RMSD metrics quantify geometric differences between predicted models. RMSD alone does NOT prove protein destabilization, altered folding kinetics, or loss of physiological function. Thermodynamic stability changes (ΔΔG) cannot be derived directly from RMSD.</div>', unsafe_allow_html=True)

    render_step_footer(0)

# ==========================================
# PAGE 2: FASTA INPUT
# ==========================================
elif page.startswith("2."):
    render_pipeline_stepper(1)
    st.markdown("### Reference and Mutant FASTA Input")
    
    if protein_choice == "Custom Upload":
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Reference Sequence")
            up_ref = st.file_uploader("Upload Reference FASTA (.fasta, .fa)", type=["fasta", "fa", "txt"], key="ref_up")
            if up_ref:
                recs = SequenceService.parse_fasta(up_ref.read().decode("utf-8"))
                if recs:
                    st.session_state.ref_seq = recs[0]["sequence"]
            ref_input = st.text_area("Reference Sequence (Paste)", value=st.session_state.ref_seq, height=220)
            st.session_state.ref_seq = ref_input
            
        with col2:
            st.subheader("Mutant Sequence")
            up_mut = st.file_uploader("Upload Mutant FASTA (.fasta, .fa)", type=["fasta", "fa", "txt"], key="mut_up")
            if up_mut:
                recs = SequenceService.parse_fasta(up_mut.read().decode("utf-8"))
                if recs:
                    st.session_state.mut_seq = recs[0]["sequence"]
            mut_input = st.text_area("Mutant Sequence (Paste)", value=st.session_state.mut_seq, height=220)
            st.session_state.mut_seq = mut_input
    else:
        st.info(f"Loaded benchmark dataset: **{protein_choice}** with mutant file **{st.session_state.selected_mutant_file}**")
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Reference Sequence")
            st.caption(f"UniProt ID: `{DATASET_MAP[protein_choice]['uniprot']}` | Length: {len(st.session_state.ref_seq)} residues")
            st.text_area("FASTA Content", value=st.session_state.ref_seq, height=260, disabled=True)
        with col2:
            st.subheader("Mutant Sequence")
            st.caption(f"Filename: `{st.session_state.selected_mutant_file}` | Length: {len(st.session_state.mut_seq)} residues")
            st.text_area("FASTA Content", value=st.session_state.mut_seq, height=260, disabled=True)
            
    # Sequence Validation Section
    st.markdown("### Sequence Validation Results")
    val_c1, val_c2 = st.columns(2)
    with val_c1:
        if st.session_state.ref_seq:
            v_ref = BackendClient.validate_sequence(st.session_state.ref_seq)
            if v_ref["is_valid"]:
                st.success(f"✅ Reference valid: {v_ref['length']} standard amino acid residues.")
            else:
                st.error(f"❌ Reference invalid: {', '.join(v_ref['errors'])}")
            if v_ref.get("warnings"):
                st.warning(f"⚠️ {', '.join(v_ref['warnings'])}")
    with val_c2:
        if st.session_state.mut_seq:
            v_mut = BackendClient.validate_sequence(st.session_state.mut_seq)
            if v_mut["is_valid"]:
                st.success(f"✅ Mutant valid: {v_mut['length']} standard amino acid residues.")
            else:
                st.error(f"❌ Mutant invalid: {', '.join(v_mut['errors'])}")
            if v_mut.get("warnings"):
                st.warning(f"⚠️ {', '.join(v_mut['warnings'])}")

    render_step_footer(1)

# ==========================================
# PAGE 3: SEQUENCE ALIGNMENT & MUTATIONS
# ==========================================
elif page.startswith("3."):
    render_pipeline_stepper(2)
    st.markdown("### Global Sequence Alignment & Mutation Detection")
    
    if not st.session_state.ref_seq or not st.session_state.mut_seq:
        st.warning("Please load or select reference and mutant sequences in Page 2.")
    else:
        align_res = BackendClient.align_and_detect_mutations(
            st.session_state.ref_seq,
            st.session_state.mut_seq,
            st.session_state.selected_protein
        )
        
        # Metric indicators
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Sequence Identity", f"{align_res['sequence_identity']:.2f}%")
        c2.metric("Alignment Score", f"{align_res['alignment_score']:.1f}")
        c3.metric("Reference Length", f"{align_res['ref_length']} aa")
        c4.metric("Mutant Length", f"{align_res['mut_length']} aa")
        
        st.markdown("#### Detected Mutations Summary")
        muts = align_res.get("mutations_detected", [])
        if muts:
            df_muts = pd.DataFrame(muts)
            st.dataframe(df_muts, use_container_width=True, hide_index=True)
            
            # Download CSV
            csv_data = align_res.get("csv_summary", "")
            st.download_button(
                "📥 Download Mutation Summary (CSV)",
                data=csv_data,
                file_name=f"{st.session_state.selected_protein}_mutations.csv",
                mime="text/csv"
            )
        else:
            st.info("No sequence mutations detected (sequences are identical).")
            
        st.markdown("#### Pairwise Alignment Viewer")
        st.code(align_res.get("alignment_repr", ""), language="text")

    render_step_footer(2)

# ==========================================
# PAGE 4: TERTIARY STRUCTURE PREDICTION
# ==========================================
elif page.startswith("4."):
    render_pipeline_stepper(3)
    st.markdown("### Tertiary Structure Prediction (ESMFold Hosted Service)")
    
    p_name = st.session_state.selected_protein
    struct_dir = os.path.join(APP_DIR, "data", "structures", p_name)
    
    st.info(
        "Tertiary structures are predicted using the hosted ESMFold API "
        "(`https://api.esmatlas.com/foldSequence/v1/pdb/`). Structures are validated and stored locally."
    )
    
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Reference Structure")
        ref_pdb = os.path.join(struct_dir, "reference.pdb")
        if os.path.exists(ref_pdb):
            valid, msg, ca_cnt = StructureService.validate_pdb(ref_pdb)
            st.success(f"✅ Ready: {os.path.basename(ref_pdb)} ({ca_cnt} C-alpha residues)")
            st.caption(f"Path: `{ref_pdb}`")
            with open(ref_pdb, "r") as f:
                st.download_button("📥 Download Reference PDB", f.read(), file_name=f"{p_name}_reference.pdb")
        else:
            st.warning("Structure not yet predicted.")
            if st.button("Predict Reference Structure via ESMFold"):
                with st.spinner("Submitting to ESMFold..."):
                    try:
                        seq = st.session_state.ref_seq
                        if len(seq) > 400:
                            seq = seq[:400]
                        StructureService.predict_with_esmfold(seq, ref_pdb)
                        st.success("Reference structure successfully predicted and saved!")
                        st.rerun()
                    except Exception as e:
                        st.error(f"Prediction failed: {e}")
                        
    with col2:
        st.subheader("Mutant Structure")
        ref_pdb, mut_pdb, _ = get_mutant_pdb_paths(struct_dir, st.session_state.selected_mutant_file)
                
        if mut_pdb and os.path.exists(mut_pdb):
            valid, msg, ca_cnt = StructureService.validate_pdb(mut_pdb)
            st.success(f"✅ Ready: {os.path.basename(mut_pdb)} ({ca_cnt} C-alpha residues)")
            st.caption(f"Path: `{mut_pdb}`")
            with open(mut_pdb, "r") as f:
                st.download_button("📥 Download Mutant PDB", f.read(), file_name=os.path.basename(mut_pdb))
        else:
            st.warning("Structure not yet predicted.")
            if st.button("Predict Mutant Structure via ESMFold"):
                with st.spinner("Submitting to ESMFold..."):
                    try:
                        seq = st.session_state.mut_seq
                        if len(seq) > 400:
                            seq = seq[:400]
                        mut_tag = "mutant"
                        if st.session_state.selected_mutant_file:
                            parts = os.path.basename(st.session_state.selected_mutant_file).replace(".fasta", "").split("_")
                            if len(parts) >= 2:
                                mut_tag = parts[1]
                        out_p = os.path.join(struct_dir, f"{mut_tag}_mutant.pdb")
                        StructureService.predict_with_esmfold(seq, out_p)
                        st.success("Mutant structure successfully predicted and saved!")
                        st.rerun()
                    except Exception as e:
                        st.error(f"Prediction failed: {e}")
                        
    # pLDDT Confidence Section
    st.markdown("### Prediction Confidence (pLDDT Analysis)")
    if ref_pdb and os.path.exists(ref_pdb):
        plddt_ref = AnalysisService.extract_plddt_scores(ref_pdb)
        if plddt_ref.get("mean_plddt"):
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Mean pLDDT", f"{plddt_ref['mean_plddt']:.1f}")
            c2.metric("Median pLDDT", f"{plddt_ref['median_plddt']:.1f}")
            c3.metric("Min pLDDT", f"{plddt_ref['min_plddt']:.1f}")
            c4.metric("Max pLDDT", f"{plddt_ref['max_plddt']:.1f}")
            
            # Confidence distribution chart
            bands = plddt_ref["confidence_bands"]
            fig = px.bar(
                x=["Very High (>90)", "Confident (70-90)", "Low (50-70)", "Very Low (<50)"],
                y=[bands["very_high_plddt_gt90"], bands["confident_plddt_70_90"], bands["low_confidence_plddt_50_70"], bands["very_low_confidence_plddt_lt50"]],
                labels={"x": "pLDDT Confidence Band", "y": "Percentage of Residues (%)"},
                title="Residue pLDDT Confidence Distribution",
                color=["#00539C", "#65AFFF", "#FFDF00", "#FF4500"],
                color_discrete_map="identity"
            )
            st.plotly_chart(fig, use_container_width=True)

    render_step_footer(3)

# ==========================================
# PAGE 5 & 6: SUPERIMPOSITION & RMSD
# ==========================================
elif page.startswith("5."):
    render_pipeline_stepper(4)
    st.markdown("### Structural Superimposition & RMSD Quantification")
    st.markdown("""
    Structural superimposition is executed automatically using Python via **Bio.PDB.Superimposer** and NumPy 
    implementing the **Kabsch algorithm**. Mismatched residues and gaps are excluded by sequence-based alignment.
    """)
    
    p_name = st.session_state.selected_protein
    struct_dir = os.path.join(APP_DIR, "data", "structures", p_name)
    ref_pdb, mut_pdb, existing_sup = get_mutant_pdb_paths(struct_dir, st.session_state.selected_mutant_file)
                
    if not os.path.exists(ref_pdb) or not mut_pdb or not os.path.exists(mut_pdb):
        st.warning("Please ensure structures are generated first (use '⚡ Run Full Pipeline Now' in sidebar).")
    else:
        mut_base = os.path.splitext(os.path.basename(mut_pdb))[0]
        sup_pdb = os.path.join(struct_dir, f"{mut_base}_superimposed.pdb")
        
        # Determine mutation reference position from current analysis or alignment
        local_pos = None
        if st.session_state.current_analysis:
            cur_muts = st.session_state.current_analysis.get("alignment", {}).get("mutations", [])
            if cur_muts:
                local_pos = cur_muts[0].get("ref_pos")
        elif st.session_state.ref_seq and st.session_state.mut_seq:
            try:
                aln_quick = SequenceService.align_sequences(st.session_state.ref_seq, st.session_state.mut_seq, p_name)
                if aln_quick.get("mutations"):
                    local_pos = aln_quick["mutations"][0].get("ref_pos")
            except Exception:
                pass

        # Trigger superimposition
        sup_res = SuperimpositionService.superimpose_structures(
            ref_pdb_path=ref_pdb,
            mut_pdb_path=mut_pdb,
            output_superimposed_path=sup_pdb,
            mutation_ref_pos=local_pos
        )
        st.session_state.ref_pdb_path = ref_pdb
        st.session_state.mut_pdb_path = mut_pdb
        st.session_state.sup_pdb_path = sup_pdb
        
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Global C-alpha RMSD", f"{sup_res['global_rmsd']:.4f} Å")
        c2.metric("Local RMSD", f"{sup_res['local_rmsd']:.4f} Å" if sup_res.get('local_rmsd') is not None else "N/A")
        c3.metric("Alignment Coverage", f"{sup_res['alignment_coverage']:.1f}%")
        c4.metric("TM-Score (est.)", f"{sup_res['tm_score']:.4f}" if sup_res.get('tm_score') else "N/A")
        
        st.markdown("#### Mathematical Formulation")
        st.latex(r"\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^N \|X_i - Y_i\|^2}")
        st.caption("Where $N$ is the number of matched C-alpha atoms, and $X_i, Y_i$ represent the coordinates after optimal Kabsch rigid-body transformation.")
        
        # Details Table
        df_details = pd.DataFrame([{
            "Reference Structure": os.path.basename(ref_pdb),
            "Mutant Structure": os.path.basename(mut_pdb),
            "Matched CA Atoms": sup_res["matched_ca_atoms_count"],
            "Global RMSD (Å)": sup_res["global_rmsd"],
            "Local RMSD (Å)": sup_res.get("local_rmsd", "N/A"),
            "Alignment Coverage": f"{sup_res['alignment_coverage']}%",
            "Superimposed PDB File": os.path.basename(sup_pdb)
        }])
        st.dataframe(df_details, use_container_width=True, hide_index=True)
        
        with open(sup_pdb, "r") as f:
            st.download_button(
                "📥 Download Superimposed Mutant PDB",
                f.read(),
                file_name=os.path.basename(sup_pdb),
                mime="chemical/x-pdb"
            )

    render_step_footer(4)

# ==========================================
# PAGE 7: INTERACTIVE 3D VISUALIZATION
# ==========================================
elif page.startswith("6."):
    render_pipeline_stepper(5)
    st.markdown("### Interactive 3D Protein Visualization")
    
    p_name = st.session_state.selected_protein
    struct_dir = os.path.join(APP_DIR, "data", "structures", p_name)
    ref_pdb, mut_pdb, sup_pdb = get_mutant_pdb_paths(struct_dir, st.session_state.selected_mutant_file)
                
    if not os.path.exists(ref_pdb) or not mut_pdb or not os.path.exists(mut_pdb):
        st.warning("Structures not found. Please click '⚡ Run Full Pipeline Now' in the sidebar.")
    else:
        # Viewer Controls
        ctrl_c1, ctrl_c2, ctrl_c3, ctrl_c4 = st.columns(4)
        with ctrl_c1:
            view_mode = st.selectbox(
                "Display Mode:",
                ["superimposed", "reference", "mutant"],
                format_func=lambda x: "Superimposed Overlay (Ref + Mut)" if x == "superimposed" else f"{x.capitalize()} Only"
            )
        with ctrl_c2:
            style_type = st.selectbox("Representation:", ["cartoon", "stick", "sphere"])
        with ctrl_c3:
            color_scheme = st.selectbox("Color Scheme:", ["protein", "plddt"], format_func=lambda x: "By Chain / Protein" if x == "protein" else "By pLDDT Confidence")
        with ctrl_c4:
            highlight_pos = st.number_input("Highlight Residue #:", min_value=1, max_value=2000, value=7 if p_name == "HBB" else 1)
            
        ref_text = open(ref_pdb).read() if os.path.exists(ref_pdb) else ""
        mut_text = open(mut_pdb).read() if os.path.exists(mut_pdb) else ""
        sup_text = open(sup_pdb).read() if sup_pdb and os.path.exists(sup_pdb) else mut_text
        
        # Render 3D viewer
        render_3d_viewer(
            ref_pdb_content=ref_text,
            mut_pdb_content=mut_text,
            superimposed_mut_content=sup_text,
            view_mode=view_mode,
            style=style_type,
            color_scheme=color_scheme,
            mutation_res_num=highlight_pos,
            height=580
        )
        
        # Downloads row
        dl_c1, dl_c2, dl_c3 = st.columns(3)
        with dl_c1:
            st.download_button("📥 Download Original Reference PDB", ref_text, file_name=f"{p_name}_reference.pdb")
        with dl_c2:
            st.download_button("📥 Download Original Mutant PDB", mut_text, file_name=os.path.basename(mut_pdb))
        with dl_c3:
            st.download_button("📥 Download Superimposed Mutant PDB", sup_text, file_name=os.path.basename(sup_pdb or mut_pdb))

    render_step_footer(5)

# ==========================================
# PAGE 8: RESULTS & BIOPHYSICAL DASHBOARD
# ==========================================
elif page.startswith("7."):
    render_pipeline_stepper(6)
    st.markdown("### Structural Effect Analysis & Biophysical Dashboard")
    
    if not st.session_state.current_analysis:
        st.info("No active analysis in memory. Please click '⚡ Run Full Pipeline Now' in the sidebar.")
    else:
        ana = st.session_state.current_analysis.get("analysis", {}).get("analysis", {})
        p_name = st.session_state.selected_protein
        
        # Biophysical Changes Table
        st.markdown("#### Biophysical Property Alterations")
        bio_changes = ana.get("biophysical_changes", [])
        if bio_changes:
            df_bio = pd.DataFrame(bio_changes)
            st.dataframe(df_bio, use_container_width=True, hide_index=True)
            
            # Biophysical delta bar chart
            deltas = []
            for b in bio_changes:
                if b.get("delta_hydropathy") is not None:
                    deltas.append({"Property": "ΔHydropathy (Kyte-Doolittle)", "Value": b["delta_hydropathy"]})
                    deltas.append({"Property": "ΔCharge", "Value": b["delta_charge"]})
                    deltas.append({"Property": "ΔVolume (Å³)", "Value": b["delta_volume"]})
                    deltas.append({"Property": "ΔMolecular Weight (Da)", "Value": b["delta_molecular_weight"]})
            if deltas:
                fig_bio = px.bar(
                    pd.DataFrame(deltas),
                    x="Property",
                    y="Value",
                    color="Property",
                    title="Biophysical Delta (Mutant - Reference)"
                )
                st.plotly_chart(fig_bio, use_container_width=True)
                
        # Scientific Interpretation
        st.markdown("#### Scientific Interpretation & Discussion")
        disc = ana.get("discussion", "")
        st.markdown(disc)
        
        # Mandatory Limitations Box
        st.markdown(f'<div class="disclaimer-box">{ana.get("disclaimer", "")}</div>', unsafe_allow_html=True)

    render_step_footer(6)

# ==========================================
# PAGE 9: PDF REPORT GENERATION
# ==========================================
elif page.startswith("8."):
    render_pipeline_stepper(7)
    st.markdown("### Publication-Grade PDF Report Generation")
    
    reports_dir = os.path.join(APP_DIR, "data", "reports")
    pdf_files = [f for f in os.listdir(reports_dir) if f.endswith(".pdf")] if os.path.exists(reports_dir) else []
    
    if pdf_files:
        st.success(f"Generated Reports Available: {len(pdf_files)}")
        for pdf_f in reversed(pdf_files[-5:]):
            p_path = os.path.join(reports_dir, pdf_f)
            c1, c2 = st.columns([3, 1])
            c1.markdown(f"📄 **{pdf_f}** ({os.path.getsize(p_path) / 1024:.1f} KB)")
            with open(p_path, "rb") as f:
                c2.download_button(
                    label="📥 Download PDF",
                    data=f.read(),
                    file_name=pdf_f,
                    mime="application/pdf",
                    key=f"dl_{pdf_f}"
                )
    else:
        st.info("No PDF reports generated yet. Click '⚡ Run Full Pipeline Now' in the sidebar to generate a report.")
    render_step_footer(7)
