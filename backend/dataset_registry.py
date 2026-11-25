"""
Centralized Dataset Registry for Protein Mutation Analysis Suite.
Manages all benchmark and extended protein datasets, UniProt accessions,
mature-protein chain offsets, mutant FASTA filenames, and clinical annotations.
"""

from typing import Dict, Any, List, Optional

# Full registry of all 28 reference proteins and their mutants
DATASET_REGISTRY: Dict[str, Dict[str, Any]] = {
    # =========================================================================
    # 1. ORIGINAL BENCHMARK PROTEINS (Tiers 1-4)
    # =========================================================================
    "HBB": {
        "name": "HBB",
        "uniprot_id": "P68871",
        "description": "Hemoglobin subunit beta (Homo sapiens)",
        "category": "Original Benchmark",
        "mature_offset": -1,  # Precursor Met1 cleaved off; mature Val1 ... Glu6 is precursor E7
        "ref_file": "HBB_reference.fasta",
        "ref_length": 147,
        "mutants": [
            {
                "mutation": "E6V",
                "label": "E6V (Sickle Cell Anemia - Original)",
                "filename": "HBB_E6V_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 7,
                "mature_pos": 6,
                "clinical_note": "Sickle cell disease (mature numbering)",
                "is_original": True,
                "available": True
            },
            {
                "mutation": "E6K",
                "label": "E6K (Hemoglobin C Disease)",
                "filename": "HBB_E6K_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 7,
                "mature_pos": 6,
                "clinical_note": "HbC (mature numbering)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "E26K",
                "label": "E26K (Hemoglobin E Disease)",
                "filename": "HBB_E26K_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 27,
                "mature_pos": 26,
                "clinical_note": "HbE (mature numbering)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "E121Q",
                "label": "E121Q (Hb D-Punjab Variant)",
                "filename": "HBB_E121Q_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 122,
                "mature_pos": 121,
                "clinical_note": "Hb D-Punjab (mature numbering)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "CFTR": {
        "name": "CFTR",
        "uniprot_id": "P13569",
        "description": "Cystic fibrosis transmembrane conductance regulator (Homo sapiens)",
        "category": "Original Benchmark",
        "mature_offset": 0,
        "ref_file": "CFTR_reference.fasta",
        "ref_length": 1480,
        "mutants": [
            {
                "mutation": "F508del",
                "label": "F508del (In-frame Deletion - Original)",
                "filename": "CFTR_F508del_mutant.fasta",
                "type": "DELETION",
                "uniprot_pos": 508,
                "mature_pos": 508,
                "clinical_note": "Primary cause of Cystic Fibrosis (~70% alleles, ERAD degradation)",
                "is_original": True,
                "available": True
            },
            {
                "mutation": "G551D",
                "label": "G551D (Class III Gating - Original)",
                "filename": "CFTR_G551D_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 551,
                "mature_pos": 551,
                "clinical_note": "Cystic fibrosis ABC signature gating defect (ivacaftor-responsive)",
                "is_original": True,
                "available": True
            },
            {
                "mutation": "R117H",
                "label": "R117H (CFTR-related disorder / Mild CF)",
                "filename": "CFTR_R117H_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 117,
                "mature_pos": 117,
                "clinical_note": "Cystic fibrosis / CFTR-related disorder (pore gating attenuation)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "EGFR": {
        "name": "EGFR",
        "uniprot_id": "P00533",
        "description": "Epidermal growth factor receptor (Homo sapiens)",
        "category": "Original Benchmark",
        "mature_offset": 0,
        "ref_file": "EGFR_reference.fasta",
        "ref_length": 1210,
        "mutants": [
            {
                "mutation": "L858R",
                "label": "L858R (Kinase Activation Driver - Original)",
                "filename": "EGFR_L858R_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 858,
                "mature_pos": 858,
                "clinical_note": "Non-small-cell lung cancer (activation loop driver)",
                "is_original": True,
                "available": True
            },
            {
                "mutation": "T790M",
                "label": "T790M (Gatekeeper Resistance Mutation)",
                "filename": "EGFR_T790M_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 790,
                "mature_pos": 790,
                "clinical_note": "Acquired resistance to 1st/2nd gen TKIs in lung cancer",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "G719S",
                "label": "G719S (P-Loop Activation Variant)",
                "filename": "EGFR_G719S_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 719,
                "mature_pos": 719,
                "clinical_note": "Non-small-cell lung cancer (exon 18 P-loop variant)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "TP53": {
        "name": "TP53",
        "uniprot_id": "P04637",
        "description": "Cellular tumor antigen p53 (Homo sapiens)",
        "category": "Original Benchmark",
        "mature_offset": 0,
        "ref_file": "TP53_reference.fasta",
        "ref_length": 393,
        "mutants": [
            {
                "mutation": "R175H",
                "label": "R175H (Conformational Zinc Hotspot - Original)",
                "filename": "TP53_R175H_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 175,
                "mature_pos": 175,
                "clinical_note": "Li-Fraumeni / somatic cancers (zinc coordination disruption)",
                "is_original": True,
                "available": True
            },
            {
                "mutation": "R273H",
                "label": "R273H (Direct DNA Contact Hotspot - Original)",
                "filename": "TP53_R273H_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 273,
                "mature_pos": 273,
                "clinical_note": "Multiple cancers (direct DNA phosphate backbone contact mutation)",
                "is_original": True,
                "available": True
            },
            {
                "mutation": "R248Q",
                "label": "R248Q (DNA Minor-Groove Contact Hotspot)",
                "filename": "TP53_R248Q_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 248,
                "mature_pos": 248,
                "clinical_note": "Multiple cancers (L3 loop minor groove contact)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "R249S",
                "label": "R249S (Aflatoxin B1 / Hepatocellular Hotspot)",
                "filename": "TP53_R249S_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 249,
                "mature_pos": 249,
                "clinical_note": "Hepatocellular carcinoma associated with aflatoxin exposure",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "Y220C",
                "label": "Y220C (Surface Cavity Destabilization Hotspot)",
                "filename": "TP53_Y220C_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 220,
                "mature_pos": 220,
                "clinical_note": "Multiple cancers (creates targetable surface crevice)",
                "is_original": False,
                "available": True
            }
        ]
    },

    # =========================================================================
    # 2. EXTENDED / NEWLY ADDED PROTEIN DATASETS
    # =========================================================================
    "KRAS": {
        "name": "KRAS",
        "uniprot_id": "P01116",
        "description": "GTPase KRas (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "KRAS_reference.fasta",
        "ref_length": 189,
        "mutants": [
            {
                "mutation": "G12D",
                "label": "G12D (Pancreatic / Colorectal Driver)",
                "filename": "KRAS_G12D_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 12,
                "mature_pos": 12,
                "clinical_note": "Pancreatic / colorectal cancers (impairs intrinsic GTP hydrolysis)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "G12V",
                "label": "G12V (GTPase-Impairment Driver)",
                "filename": "KRAS_G12V_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 12,
                "mature_pos": 12,
                "clinical_note": "Several cancers (colorectal, lung, pancreatic)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "G12C",
                "label": "G12C (Sotorasib / Adagrasib Target)",
                "filename": "KRAS_G12C_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 12,
                "mature_pos": 12,
                "clinical_note": "Non-small cell lung cancer (covalent switch-II pocket target)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "G13D",
                "label": "G13D (Colorectal Adenocarcinoma Driver)",
                "filename": "KRAS_G13D_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 13,
                "mature_pos": 13,
                "clinical_note": "Colorectal cancer",
                "is_original": False,
                "available": True
            }
        ]
    },
    "BRAF": {
        "name": "BRAF",
        "uniprot_id": "P15056",
        "description": "Serine/threonine-protein kinase B-raf (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "BRAF_reference.fasta",
        "ref_length": 766,
        "mutants": [
            {
                "mutation": "V600E",
                "label": "V600E (Melanoma Constitutive Kinase Driver)",
                "filename": "BRAF_V600E_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 600,
                "mature_pos": 600,
                "clinical_note": "Melanoma, colorectal, thyroid cancers (Vemurafenib-sensitive)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "V600K",
                "label": "V600K (Melanoma Activating Mutation)",
                "filename": "BRAF_V600K_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 600,
                "mature_pos": 600,
                "clinical_note": "Melanoma (constitutive monomer activation)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "PIK3CA": {
        "name": "PIK3CA",
        "uniprot_id": "P42336",
        "description": "Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "PIK3CA_reference.fasta",
        "ref_length": 1068,
        "mutants": [
            {
                "mutation": "E545K",
                "label": "E545K (Helical Domain Oncogenic Driver)",
                "filename": "PIK3CA_E545K_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 545,
                "mature_pos": 545,
                "clinical_note": "Breast, colorectal, and endometrial cancers (alpelisib-sensitive)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "H1047R",
                "label": "H1047R (Kinase Domain Activation Driver)",
                "filename": "PIK3CA_H1047R_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 1047,
                "mature_pos": 1047,
                "clinical_note": "Breast cancer activation loop hotspot",
                "is_original": False,
                "available": True
            }
        ]
    },
    "KIT": {
        "name": "KIT",
        "uniprot_id": "P10721",
        "description": "Mast/stem cell growth factor receptor Kit (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "KIT_reference.fasta",
        "ref_length": 976,
        "mutants": [
            {
                "mutation": "D816V",
                "label": "D816V (Systemic Mastocytosis / GIST Driver)",
                "filename": "KIT_D816V_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 816,
                "mature_pos": 816,
                "clinical_note": "Systemic mastocytosis, GIST (imatinib-resistant, avapritinib-sensitive)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "JAK2": {
        "name": "JAK2",
        "uniprot_id": "O60674",
        "description": "Tyrosine-protein kinase JAK2 (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "JAK2_reference.fasta",
        "ref_length": 1132,
        "mutants": [
            {
                "mutation": "V617F",
                "label": "V617F (Polycythemia Vera / MPN Driver)",
                "filename": "JAK2_V617F_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 617,
                "mature_pos": 617,
                "clinical_note": "Myeloproliferative neoplasms (disrupts autoinhibitory JH2 domain)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "IDH1": {
        "name": "IDH1",
        "uniprot_id": "O75874",
        "description": "Isocitrate dehydrogenase [NADP] cytoplasmic (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "IDH1_reference.fasta",
        "ref_length": 414,
        "mutants": [
            {
                "mutation": "R132H",
                "label": "R132H (Glioma / AML Oncometabolite Producer)",
                "filename": "IDH1_R132H_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 132,
                "mature_pos": 132,
                "clinical_note": "Glioma, AML, chondrosarcoma (neomorphic synthesis of D-2-HG)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "IDH2": {
        "name": "IDH2",
        "uniprot_id": "P48735",
        "description": "Isocitrate dehydrogenase [NADP], mitochondrial (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "IDH2_reference.fasta",
        "ref_length": 452,
        "mutants": [
            {
                "mutation": "R140Q",
                "label": "R140Q (AML / Glioma Oncometabolite Producer)",
                "filename": "IDH2_R140Q_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 140,
                "mature_pos": 140,
                "clinical_note": "Acute myeloid leukemia (enasidenib-responsive)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "FGFR3": {
        "name": "FGFR3",
        "uniprot_id": "P22607",
        "description": "Fibroblast growth factor receptor 3 (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "FGFR3_reference.fasta",
        "ref_length": 806,
        "mutants": [
            {
                "mutation": "G380R",
                "label": "G380R (Achondroplasia Genetic Driver)",
                "filename": "FGFR3_G380R_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 380,
                "mature_pos": 380,
                "clinical_note": "Achondroplasia (stabilizes transmembrane helix dimer)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "RET": {
        "name": "RET",
        "uniprot_id": "P07949",
        "description": "Proto-oncogene tyrosine-protein kinase receptor Ret (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "RET_reference.fasta",
        "ref_length": 1114,
        "mutants": [
            {
                "mutation": "M918T",
                "label": "M918T (MEN2B / Medullary Thyroid Carcinoma)",
                "filename": "RET_M918T_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 918,
                "mature_pos": 918,
                "clinical_note": "Multiple endocrine neoplasia type 2B (alters substrate specificity)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "SOD1": {
        "name": "SOD1",
        "uniprot_id": "P00441",
        "description": "Superoxide dismutase [Cu-Zn] (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": -1,  # Cleaved Met1; precursor A5 is mature A4
        "ref_file": "SOD1_reference.fasta",
        "ref_length": 154,
        "mutants": [
            {
                "mutation": "A4V",
                "label": "A4V (Familial ALS Severe Rapid Progression)",
                "filename": "SOD1_A4V_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 5,
                "mature_pos": 4,
                "clinical_note": "Familial amyotrophic lateral sclerosis (dimer destabilization)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "G93A",
                "label": "G93A (Classical ALS Benchmark Variant)",
                "filename": "SOD1_G93A_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 94,
                "mature_pos": 93,
                "clinical_note": "Familial ALS model mutation (accelerates neurotoxic misfolding)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "MAPT": {
        "name": "MAPT",
        "uniprot_id": "P10636-8",
        "description": "Microtubule-associated protein tau (Isoform 441-aa, Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "MAPT_reference.fasta",
        "ref_length": 441,
        "mutants": [
            {
                "mutation": "P301L",
                "label": "P301L (Frontotemporal Dementia / Tauopathy)",
                "filename": "MAPT_P301L_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 301,
                "mature_pos": 301,
                "clinical_note": "Frontotemporal dementia with parkinsonism (accelerates filament assembly)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "SNCA": {
        "name": "SNCA",
        "uniprot_id": "P37840",
        "description": "Alpha-synuclein (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "SNCA_reference.fasta",
        "ref_length": 140,
        "mutants": [
            {
                "mutation": "A53T",
                "label": "A53T (Familial Early-Onset Parkinson Disease)",
                "filename": "SNCA_A53T_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 53,
                "mature_pos": 53,
                "clinical_note": "Autosomal dominant Parkinson disease (accelerates protofibril kinetics)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "A30P",
                "label": "A30P (Familial Parkinson Disease Variant)",
                "filename": "SNCA_A30P_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 30,
                "mature_pos": 30,
                "clinical_note": "Familial Parkinson disease (disrupts membrane helical affinity)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "APP": {
        "name": "APP",
        "uniprot_id": "P05067",
        "description": "Amyloid-beta precursor protein (Isoform APP770, Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "APP_reference.fasta",
        "ref_length": 770,
        "mutants": [
            {
                "mutation": "V717I",
                "label": "V717I (London Mutation / Familial Alzheimer)",
                "filename": "APP_V717I_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 717,
                "mature_pos": 717,
                "clinical_note": "Familial Alzheimer disease (shifts gamma-secretase cleavage to toxic Abeta42)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "PSEN1": {
        "name": "PSEN1",
        "uniprot_id": "P49768",
        "description": "Presenilin-1 (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "PSEN1_reference.fasta",
        "ref_length": 467,
        "mutants": [
            {
                "mutation": "M146L",
                "label": "M146L (Early-Onset Alzheimer Disease)",
                "filename": "PSEN1_M146L_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 146,
                "mature_pos": 146,
                "clinical_note": "Autosomal dominant early-onset Alzheimer disease",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "E280A",
                "label": "E280A (Paisa Mutation / Familial Alzheimer)",
                "filename": "PSEN1_E280A_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 280,
                "mature_pos": 280,
                "clinical_note": "Familial Alzheimer disease (Antioquia/Paisa kindred hotspot)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "PAH": {
        "name": "PAH",
        "uniprot_id": "P00439",
        "description": "Phenylalanine-4-hydroxylase (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "PAH_reference.fasta",
        "ref_length": 452,
        "mutants": [
            {
                "mutation": "R408W",
                "label": "R408W (Classic Severe Phenylketonuria)",
                "filename": "PAH_R408W_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 408,
                "mature_pos": 408,
                "clinical_note": "Classic phenylketonuria (catalytic domain aggregation)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "Y414C",
                "label": "Y414C (BH4-Responsive Phenylketonuria)",
                "filename": "PAH_Y414C_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 414,
                "mature_pos": 414,
                "clinical_note": "Mild PKU / BH4-responsive hyperphenylalaninemia",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "R261Q",
                "label": "R261Q (Intermediate Phenylketonuria)",
                "filename": "PAH_R261Q_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 261,
                "mature_pos": 261,
                "clinical_note": "Phenylketonuria (unstable catalytic core fold)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "HFE": {
        "name": "HFE",
        "uniprot_id": "Q30201",
        "description": "Hereditary hemochromatosis protein (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "HFE_reference.fasta",
        "ref_length": 348,
        "mutants": [
            {
                "mutation": "C282Y",
                "label": "C282Y (Primary Hereditary Hemochromatosis)",
                "filename": "HFE_C282Y_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 282,
                "mature_pos": 282,
                "clinical_note": "Hereditary hemochromatosis (disrupts beta-2-microglobulin disulfide bond)",
                "is_original": False,
                "available": True
            },
            {
                "mutation": "H63D",
                "label": "H63D (Low-Penetrance Iron Overload Variant)",
                "filename": "HFE_H63D_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 63,
                "mature_pos": 63,
                "clinical_note": "Iron overload risk cofactor",
                "is_original": False,
                "available": True
            }
        ]
    },
    "FBN1": {
        "name": "FBN1",
        "uniprot_id": "P35555",
        "description": "Fibrillin-1 (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "FBN1_reference.fasta",
        "ref_length": 2871,
        "mutants": [
            {
                "mutation": "C1039G",
                "label": "C1039G (Classic Marfan Syndrome Driver)",
                "filename": "FBN1_C1039G_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 1039,
                "mature_pos": 1039,
                "clinical_note": "Marfan syndrome (disrupts calcium-binding cbEGF domain disulfide)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "F9": {
        "name": "F9",
        "uniprot_id": "P00740",
        "description": "Coagulation factor IX (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": -46,  # 46-aa signal peptide and propeptide cleaved in mature plasma factor IX
        "ref_file": "F9_reference.fasta",
        "ref_length": 461,
        "mutants": [
            {
                "mutation": "R338L",
                "label": "R338L (Factor IX Padua / Thrombophilia)",
                "filename": "F9_R338L_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 384,
                "mature_pos": 338,
                "clinical_note": "Factor IX Padua (hyperactive catalytic domain with 8-fold clotting activity)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "GBA1": {
        "name": "GBA1",
        "uniprot_id": "P04062",
        "description": "Glucosylceramidase (Beta-glucocerebrosidase, Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": -39,  # 39-aa signal peptide cleaved; precursor N409 corresponds to legacy N370S
        "ref_file": "GBA1_reference.fasta",
        "ref_length": 536,
        "mutants": [
            {
                "mutation": "N409S",
                "label": "N409S (Gaucher Disease Type 1 / Legacy N370S)",
                "filename": "GBA1_N409S_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 409,
                "mature_pos": 370,
                "clinical_note": "Gaucher disease type 1 & Parkinson risk (legacy nomenclature N370S)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "HEXA": {
        "name": "HEXA",
        "uniprot_id": "P06865",
        "description": "Beta-hexosaminidase subunit alpha (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "HEXA_reference.fasta",
        "ref_length": 529,
        "mutants": [
            {
                "mutation": "R178H",
                "label": "R178H (B1 Variant Tay-Sachs Disease)",
                "filename": "HEXA_R178H_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 178,
                "mature_pos": 178,
                "clinical_note": "Tay-Sachs disease B1 variant (active-site catalytic charge disruption)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "ATP7B": {
        "name": "ATP7B",
        "uniprot_id": "P35670",
        "description": "Copper-transporting ATPase 2 (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "ATP7B_reference.fasta",
        "ref_length": 1465,
        "mutants": [
            {
                "mutation": "H1069Q",
                "label": "H1069Q (Classic Wilson Disease Hotspot)",
                "filename": "ATP7B_H1069Q_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 1069,
                "mature_pos": 1069,
                "clinical_note": "Wilson disease (ATP-binding SEHPL domain misfolding and retention)",
                "is_original": False,
                "available": True
            }
        ]
    },
    "OTC": {
        "name": "OTC",
        "uniprot_id": "P00480",
        "description": "Ornithine carbamoyltransferase, mitochondrial (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "OTC_reference.fasta",
        "ref_length": 354,
        "mutants": [
            {
                "mutation": "R40H",
                "label": "R40H (Ornithine Transcarbamylase Deficiency)",
                "filename": "OTC_R40H_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 40,
                "mature_pos": 40,
                "clinical_note": "Hyperammonemia due to urea cycle disorder",
                "is_original": False,
                "available": True
            }
        ]
    },
    "G6PD": {
        "name": "G6PD",
        "uniprot_id": "P11413",
        "description": "Glucose-6-phosphate 1-dehydrogenase (Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "G6PD_reference.fasta",
        "ref_length": 515,
        "mutants": [
            {
                "mutation": "N126D",
                "label": "N126D (G6PD A- Haplotype Component)",
                "filename": "G6PD_N126D_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": 126,
                "mature_pos": 126,
                "clinical_note": "G6PD A+ / A- African polymorphic haplotype component",
                "is_original": False,
                "available": True
            }
        ]
    },
    "RYR1": {
        "name": "RYR1",
        "uniprot_id": "P21817",
        "description": "Ryanodine receptor 1 (Skeletal muscle, Homo sapiens)",
        "category": "Extended Dataset",
        "mature_offset": 0,
        "ref_file": "RYR1_reference.fasta",
        "ref_length": 5038,
        "mutants": [
            {
                "mutation": "R615C",
                "label": "R615C (Malignant Hyperthermia - Sequence Discrepancy)",
                "filename": "RYR1_R615C_mutant.fasta",
                "type": "SUBSTITUTION",
                "uniprot_pos": None,
                "mature_pos": None,
                "clinical_note": (
                    "Malignant hyperthermia susceptibility. Note: In the canonical UniProt sequence P21817, "
                    "residue 615 is Serine (S) rather than Arginine (R) (neighbours: A612-V613-R614-S615-N616-Q617-D618). "
                    "Mutant sequence not supplied in dataset archive."
                ),
                "is_original": False,
                "available": False
            }
        ]
    }
}


def get_all_proteins() -> List[Dict[str, Any]]:
    """Returns a list of all registered proteins with their metadata."""
    return list(DATASET_REGISTRY.values())


def get_protein(gene_name: str) -> Optional[Dict[str, Any]]:
    """Returns metadata for a specific protein by gene name."""
    return DATASET_REGISTRY.get(gene_name.upper())


def get_mutants_for_protein(gene_name: str) -> List[Dict[str, Any]]:
    """Returns available mutant datasets for a given protein."""
    prot = get_protein(gene_name)
    if prot:
        return [m for m in prot.get("mutants", []) if m.get("available", True)]
    return []


def get_mature_offsets() -> Dict[str, int]:
    """Returns a mapping of gene symbols to mature-protein numbering offsets."""
    return {k: v.get("mature_offset", 0) for k, v in DATASET_REGISTRY.items()}


def get_dataset_map() -> Dict[str, Any]:
    """
    Returns the dataset mapping dictionary formatted for the Streamlit sidebar selection.
    """
    mapping = {}
    for gene, info in DATASET_REGISTRY.items():
        mut_tuples = [
            (m["label"], m["filename"])
            for m in info.get("mutants", [])
            if m.get("available", True)
        ]
        mapping[gene] = {
            "uniprot": info["uniprot_id"],
            "description": info["description"],
            "category": info["category"],
            "ref_file": info["ref_file"],
            "mutants": mut_tuples,
            "has_mutants": len(mut_tuples) > 0
        }
    return mapping


# Comprehensive biological literature context for automated scientific discussion
KNOWN_BIOLOGICAL_CONTEXT: Dict[str, str] = {
    # Original 6 benchmarks
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
    ),

    # Additional HBB mutants
    "HBB_E6K": (
        "Hemoglobin C (HbC) variant: Glu6Lys replaces a negatively charged glutamate with a positively charged lysine at mature position 6, "
        "reducing solubility in erythrocytes and causing crystallization of oxy-HbC with mild hemolytic anemia."
    ),
    "HBB_E26K": (
        "Hemoglobin E (HbE) variant: Glu26Lys mutation activates an alternate cryptic splice site in exon 1 of beta-globin mRNA, "
        "producing mild beta-thalassemic phenotype and microcytic anemia common in Southeast Asia."
    ),
    "HBB_E121Q": (
        "Hemoglobin D-Punjab (Hb D-Los Angeles): Glu121Gln substitution in the alpha-1/beta-1 interface alters electrostatic surface potential; "
        "co-inheritance with HbS causes clinically significant sickle cell disease."
    ),

    # Additional CFTR mutant
    "CFTR_R117H": (
        "Cystic Fibrosis conductance/gating mutation: Arg117His in the first transmembrane segment attenuates chloride channel open "
        "conductance and gating, often resulting in congenital bilateral absence of the vas deferens (CBAVD) or mild atypical CF."
    ),

    # Additional EGFR mutants
    "EGFR_T790M": (
        "EGFR gatekeeper resistance mutation: Thr790Met replaces the steric gatekeeper threonine with bulky methionine, sterically hindering "
        "first-generation TKIs (gefitinib/erlotinib) while increasing ATP binding affinity; targeted by 3rd-generation osimertinib."
    ),
    "EGFR_G719S": (
        "EGFR exon 18 P-loop mutation: Gly719Ser destabilizes the autoinhibitory alphaC-helix conformation in the N-terminal lobe, "
        "rendering tumors sensitive to afatinib and osimertinib."
    ),

    # Additional TP53 mutants
    "TP53_R248Q": (
        "p53 DNA-contact hotspot: Arg248 in the L3 loop projects directly into the minor groove of response elements; substitution to Gln "
        "abrogates DNA binding with dominant-negative gain-of-function oncogenic properties."
    ),
    "TP53_R249S": (
        "p53 aflatoxin B1 fingerprint mutation: Arg249Ser hotspot in hepatocellular carcinoma associated with dietary aflatoxin B1 exposure "
        "and hepatitis B virus co-infection."
    ),
    "TP53_Y220C": (
        "p53 surface cavity mutation: Tyr220Cys destabilizes the core domain by creating a 120-cubic-Angstrom hydrophobic crevice on the "
        "beta-sandwich surface, targeted by small-molecule conformational stabilizers (e.g. PC14586)."
    ),

    # KRAS
    "KRAS_G12D": (
        "KRAS G12D driver: Gly12Asp in the P-loop sterically hinders GAP-mediated GTP hydrolysis, freezing KRAS in the active GTP-bound state "
        "and driving constitutive RAF-MEK-ERK signaling in pancreatic adenocarcinoma and colorectal cancer."
    ),
    "KRAS_G12V": (
        "KRAS G12V driver: Valine substitution at Gly12 introduces steric bulk into the catalytic pocket, locking the switch-I and switch-II "
        "regions into active effector-binding conformations."
    ),
    "KRAS_G12C": (
        "KRAS G12C driver: Cysteine substitution at Gly12 provides a nucleophilic sulfhydryl group targeted by covalent allosteric inhibitors "
        "(sotorasib, adagrasib) that lock the protein into the inactive GDP-bound state in NSCLC."
    ),
    "KRAS_G13D": (
        "KRAS G13D driver: Gly13Asp impairs GAP-stimulated GTP hydrolysis, frequently observed in colorectal cancer with distinct sensitivity "
        "profiles compared to codon 12 alterations."
    ),

    # BRAF
    "BRAF_V600E": (
        "BRAF V600E driver: Val600Glu in the activation loop introduces negative charge mimicking activating phosphorylation, causing a 500-fold "
        "increase in kinase activity as an active monomer independent of RAS-dimerization; highly responsive to vemurafenib/dabrafenib."
    ),
    "BRAF_V600K": (
        "BRAF V600K driver: Lysine substitution at Val600 destabilizes inactive autoinhibitory contacts, driving potent constitutive monomeric "
        "signaling in metastatic melanoma."
    ),

    # PIK3CA
    "PIK3CA_E545K": (
        "PIK3CA helical domain hotspot: Glu545Lys disrupts inhibitory interaction with the nSH2 domain of regulatory subunit p85, promoting "
        "constitutive PI3K lipid kinase signaling in breast cancer."
    ),
    "PIK3CA_H1047R": (
        "PIK3CA kinase domain hotspot: His1047Arg in the activation loop increases affinity for lipid membranes, hyperactivating catalytic "
        "phosphorylation of PIP2 to PIP3 in hormone receptor-positive breast cancers."
    ),

    # KIT & JAK2
    "KIT_D816V": (
        "KIT D816V driver: Asp816Val in the activation loop stabilizes the active kinase state independent of stem cell factor (SCF), "
        "driving systemic mastocytosis and conferring resistance to imatinib."
    ),
    "JAK2_V617F": (
        "JAK2 V617F driver: Val617Phe in the pseudokinase (JH2) domain abrogates its physiological autoinhibition on the adjacent kinase (JH1) "
        "domain, driving autonomous erythropoiesis in polycythemia vera."
    ),

    # IDH1 & IDH2
    "IDH1_R132H": (
        "IDH1 oncometabolite driver: Arg132His in the catalytic active site confers neomorphic enzymatic activity converting alpha-ketoglutarate "
        "to the oncometabolite D-2-hydroxyglutarate (D-2-HG), causing DNA/histone hypermethylation in low-grade gliomas and AML."
    ),
    "IDH2_R140Q": (
        "IDH2 oncometabolite driver: Arg140Gln in the mitochondrial isocitrate dehydrogenase active site catalyzes production of D-2-HG, "
        "inhibiting TET methylcytosine dioxygenases and blocking hematopoietic differentiation in AML."
    ),

    # FGFR3 & RET
    "FGFR3_G380R": (
        "FGFR3 achondroplasia driver: Gly380Arg in the transmembrane domain stabilizes ligand-independent receptor homodimers, sending "
        "constitutive STAT1/MAPK inhibitory signals to growth plate chondrocytes."
    ),
    "RET_M918T": (
        "RET M918T MEN2B driver: Met918Thr in the kinase activation loop alters substrate recognition kinetics and drives ligand-independent "
        "oncogenic signaling in aggressive medullary thyroid carcinoma."
    ),

    # Neurodegenerative proteins
    "SOD1_A4V": (
        "SOD1 A4V ALS driver: Ala4Val (precursor Ala5Val) at the dimer interface dramatically destabilizes the Cu/Zn superoxide dismutase "
        "homodimer, promoting rapid toxic oligomerization and motor neuron apoptosis."
    ),
    "SOD1_G93A": (
        "SOD1 G93A ALS benchmark: Gly93Ala accelerates beta-barrel misfolding and aggregation in familial amyotrophic lateral sclerosis."
    ),
    "MAPT_P301L": (
        "Tau P301L frontotemporal dementia driver: Pro301Leu in the microtubule-binding repeat region diminishes microtubule assembly and "
        "drastically accelerates assembly into paired helical filaments and neurofibrillary tangles."
    ),
    "SNCA_A53T": (
        "Alpha-synuclein A53T Parkinson driver: Ala53Thr in the amphipathic lipid-binding domain accelerates nucleation and formation of "
        "insoluble amyloid protofibrils and Lewy bodies in familial Parkinson disease."
    ),
    "SNCA_A30P": (
        "Alpha-synuclein A30P variant: Ala30Pro disrupts membrane-associated alpha-helical conformation, shifting the equilibrium toward "
        "cytosolic oligomers."
    ),
    "APP_V717I": (
        "APP London mutation: Val717Ile near the gamma-secretase transmembrane cleavage site shifts the cleavage ratio toward the highly "
        "fibrillogenic and neurotoxic amyloid-beta 42 peptide in early-onset familial Alzheimer disease."
    ),
    "PSEN1_M146L": (
        "Presenilin-1 M146L driver: Met146Leu in transmembrane domain 2 of the gamma-secretase catalytic core alters proteolytic precision, "
        "elevating the Abeta42/Abeta40 ratio in familial Alzheimer disease."
    ),
    "PSEN1_E280A": (
        "Presenilin-1 E280A Paisa mutation: Glu280Ala hotspot in the transmembrane domain alters catalytic positioning, causing early-onset "
        "familial Alzheimer disease."
    ),

    # Metabolic diseases
    "PAH_R408W": (
        "PAH R408W PKU driver: Arg408Trp in the catalytic domain causes rapid misfolding and clearance, resulting in severe classic "
        "phenylketonuria and toxic phenylalanine accumulation."
    ),
    "PAH_Y414C": (
        "PAH Y414C variant: Tyr414Cys destabilizes the catalytic core but retains residual activity and responds clinically to sapropterin (BH4)."
    ),
    "PAH_R261Q": (
        "PAH R261Q variant: Arg261Gln destabilizes the dimer-tetramer interface in intermediate phenylketonuria."
    ),
    "HFE_C282Y": (
        "HFE C282Y hemochromatosis driver: Cys282Tyr breaks a critical disulfide bridge required for association with beta-2-microglobulin, "
        "preventing surface expression and dysregulating hepcidin-mediated iron homeostasis."
    ),
    "HFE_H63D": (
        "HFE H63D cofactor: His63Asp in the alpha-1 domain moderately alters transferrin receptor affinity, acting as a low-penetrance "
        "risk factor for systemic iron overload."
    ),
    "FBN1_C1039G": (
        "Fibrillin-1 C1039G Marfan driver: Cys1039Gly eliminates a conserved disulfide bond in a calcium-binding EGF-like domain, "
        "causing microfibril fragmentation and aortic aneurysm susceptibility."
    ),
    "F9_R338L": (
        "Factor IX Padua thrombophilia: Arg338Leu (precursor Arg384Leu) enhances catalytic activity 8-fold by optimizing interaction with "
        "factor VIIIa, causing familial thrombophilia."
    ),
    "GBA1_N409S": (
        "GBA1 N409S (legacy N370S): Asn409Ser in beta-glucocerebrosidase impairs catalytic turnover and lysosomal stability, causing "
        "glucosylceramide storage in macrophages in Gaucher disease type 1."
    ),
    "HEXA_R178H": (
        "HEXA R178H Tay-Sachs B1 variant: Arg178His alters the electrostatic charge in the catalytic pocket of beta-hexosaminidase A, "
        "preventing GM2 ganglioside degradation."
    ),
    "ATP7B_H1069Q": (
        "ATP7B H1069Q Wilson disease driver: His1069Gln in the cytoplasmic nucleotide-binding domain causes thermal instability and ER "
        "retention, blocking biliary copper excretion and causing hepatic and neurological toxicity."
    ),
    "OTC_R40H": (
        "OTC R40H deficiency: Arg40His in mitochondrial ornithine transcarbamylase disrupts the substrate-binding cleft, causing severe "
        "hyperammonemic urea cycle crisis."
    ),
    "G6PD_N126D": (
        "G6PD N126D variant: Asn126Asp is the characteristic polymorphic marker of the G6PD A variant common in African populations."
    ),
    "RYR1_R615C": (
        "RYR1 R615C Malignant Hyperthermia: Note: in canonical UniProt sequence P21817, residue 615 is Serine (S) rather than Arginine (R). "
        "This dataset exhibits an isoform/numbering discrepancy with the reference sequence."
    )
}
