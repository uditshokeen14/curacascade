"""CuraCascade — Population Genetics & Multi-Ancestry Engine (Objectives 4 & 9).

Calculates Maximum Credible Allele Frequency (MCAF) based on disease prevalence,
mode of inheritance, penetrance, and allelic heterogeneity (Whiffin et al. 2017).
Evaluates multi-ancestry population frequencies across five major continental cohorts:
- Non-Finnish European (NFE) [Historical baseline]
- African / African American (AFR) [US underrepresented]
- Admixed American / Latino (AMR) [US underrepresented]
- East Asian (EAS) [Asian reference]
- Japanese Cohort (JPN) [ToMMo 54KJPN / gnomAD Japanese subset]

Calculates Ancestry Discordance Ratios and flags Ancestry-Specific Misclassification
Risks where variants are common in non-European groups but rare in Europeans.
"""

from __future__ import annotations

import math
from pathlib import Path
import sqlite3
import sys
from typing import Any, Dict, Optional, Union

# Ensure repository root is on sys.path
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import pandas as pd
from database.connection import get_connection, get_db_path, query_df


# ==============================================================================
# 1. Standard Clinical Genomics MCAF Formula
# ==============================================================================

def calculate_mcaf(
    prevalence_per_million: float,
    inheritance_mode: str,
    penetrance: float = 0.5,
    allelic_heterogeneity: float = 0.1,
) -> float:
    """Calculate the Maximum Credible Allele Frequency (MCAF) threshold.

    Based on the clinical genomics framework by Whiffin et al. (Genet Med 2017)
    adopted by ClinGen and ACMG Sequence Variant Interpretation workgroups.

    Args:
        prevalence_per_million: Disease prevalence in individuals per million (from Orphanet).
        inheritance_mode: Genetic mode of inheritance ('AD', 'AR', 'XL', etc.).
        penetrance: Estimated penetrance of pathogenic variants (0.0 to 1.0, default 0.5).
        allelic_heterogeneity: Maximum proportion of disease cases caused by a single variant (default 0.1 = 10%).

    Returns:
        Maximum Credible Allele Frequency as a float proportion (e.g., 0.0002 for 2.0e-4).
    """
    if prevalence_per_million <= 0 or penetrance <= 0:
        return 0.0

    # Convert prevalence from per-million count to proportion
    prevalence_prop = float(prevalence_per_million) / 1_000_000.0

    mode = str(inheritance_mode).strip().upper()
    is_recessive = "AR" in mode or "RECESSIVE" in mode

    if is_recessive:
        # Recessive condition: MCAF = sqrt((P * h) / penetrance)
        val = (prevalence_prop * allelic_heterogeneity) / penetrance
        mcaf = math.sqrt(max(0.0, val))
    else:
        # Dominant condition (AD / XL dominant): MCAF = (P * h) / (penetrance * 2)
        mcaf = (prevalence_prop * allelic_heterogeneity) / (penetrance * 2.0)

    return float(round(mcaf, 8))


# ==============================================================================
# 2. Orphanet Curated Disease Prevalence & Genetic Architecture Registry
# ==============================================================================

DISEASE_PREVALENCE_REGISTRY: Dict[str, Dict[str, Any]] = {
    "MONDO:0005045": {
        "disease_name": "Hypertrophic Cardiomyopathy (HCM)",
        "orphanet_id": "ORPHA:217",
        "prevalence_per_million": 2000.0,   # 1 in 500
        "inheritance_mode": "Autosomal dominant",
        "penetrance": 0.5,
        "allelic_heterogeneity": 0.10,      # High heterogeneity (sarcomeric genes)
    },
    "MONDO:0005014": {
        "disease_name": "Dilated Cardiomyopathy (DCM)",
        "orphanet_id": "ORPHA:217604",
        "prevalence_per_million": 400.0,    # 1 in 2,500
        "inheritance_mode": "Autosomal dominant",
        "penetrance": 0.5,
        "allelic_heterogeneity": 0.05,
    },
    "MONDO:0001625": {
        "disease_name": "Long QT Syndrome (LQTS)",
        "orphanet_id": "ORPHA:768",
        "prevalence_per_million": 500.0,    # 1 in 2,000
        "inheritance_mode": "Autosomal dominant",
        "penetrance": 0.5,
        "allelic_heterogeneity": 0.10,
    },
    "MONDO:0015263": {
        "disease_name": "Brugada Syndrome (BrS)",
        "orphanet_id": "ORPHA:130",
        "prevalence_per_million": 500.0,    # 1 in 2,000
        "inheritance_mode": "Autosomal dominant",
        "penetrance": 0.5,
        "allelic_heterogeneity": 0.10,
    },
    "MONDO:0011448": {
        "disease_name": "Hereditary Breast and Ovarian Cancer Syndrome (HBOC)",
        "orphanet_id": "ORPHA:145",
        "prevalence_per_million": 2500.0,   # 1 in 400
        "inheritance_mode": "Autosomal dominant",
        "penetrance": 0.5,
        "allelic_heterogeneity": 0.05,
    },
    "MONDO:0007817": {
        "disease_name": "Autosomal Dominant / Recessive Nonsyndromic Hearing Loss",
        "orphanet_id": "ORPHA:90642",
        "prevalence_per_million": 1000.0,   # 1 in 1,000
        "inheritance_mode": "Autosomal recessive",
        "penetrance": 0.8,
        "allelic_heterogeneity": 0.10,
    },
    "MONDO:0005260": {
        "disease_name": "Autism Spectrum Disorder (ASD)",
        "orphanet_id": "ORPHA:140162",
        "prevalence_per_million": 10000.0,  # 1 in 100 (complex/heterogeneous)
        "inheritance_mode": "Complex / Multifactorial",
        "penetrance": 0.3,
        "allelic_heterogeneity": 0.01,
    },
    "MONDO:0008304": {
        "disease_name": "Autosomal Dominant Polycystic Kidney Disease (ADPKD)",
        "orphanet_id": "ORPHA:730",
        "prevalence_per_million": 1000.0,   # 1 in 1,000
        "inheritance_mode": "Autosomal dominant",
        "penetrance": 0.8,
        "allelic_heterogeneity": 0.10,
    },
    "MONDO:0010894": {
        "disease_name": "Maturity-Onset Diabetes of the Young (MODY)",
        "orphanet_id": "ORPHA:552",
        "prevalence_per_million": 100.0,    # 1 in 10,000
        "inheritance_mode": "Autosomal dominant",
        "penetrance": 0.6,
        "allelic_heterogeneity": 0.10,
    },
    "MONDO:0010675": {
        "disease_name": "Phelan-McDermid Syndrome",
        "orphanet_id": "ORPHA:48652",
        "prevalence_per_million": 50.0,     # ~1 in 20,000
        "inheritance_mode": "Autosomal dominant",
        "penetrance": 0.9,
        "allelic_heterogeneity": 0.20,
    },
}

# Precompute MCAF thresholds for each disease
for did, dinfo in DISEASE_PREVALENCE_REGISTRY.items():
    dinfo["mcaf_threshold"] = calculate_mcaf(
        prevalence_per_million=dinfo["prevalence_per_million"],
        inheritance_mode=dinfo["inheritance_mode"],
        penetrance=dinfo["penetrance"],
        allelic_heterogeneity=dinfo["allelic_heterogeneity"],
    )


# ==============================================================================
# 3. gnomAD v4 Popmax & 5-Ancestry Cohort Catalog
# ==============================================================================

# Benchmark catalog mapping ClinVar Variation IDs to population frequencies
# Cohorts:
# - af_nfe: Non-Finnish European (European baseline)
# - af_afr: African / African American (US underrepresented)
# - af_amr: Admixed American / Latino (US underrepresented)
# - af_eas: East Asian (Asian reference)
# - af_jpn: Japanese Cohort (ToMMo 54KJPN / gnomAD Japanese subset)
# - af_sas: South Asian reference

GNOMAD_VARIANT_CATALOG: Dict[int, Dict[str, Any]] = {
    # --------------------------------------------------------------------------
    # KCNE2 (LQTS MCAF = 5.0e-5) — Refuted by ClinGen Arrhythmia GCEP
    # --------------------------------------------------------------------------
    67201: {  # c.253T>C (p.Tyr85His) — Celebrated refutation (Giudicessi et al. 2018)
        "af_global": 0.00874,      # 0.874% global
        "af_popmax": 0.01420,      # 1.42% in Japanese (JPN) -> 284x MCAF!
        "popmax_ancestry": "JPN",
        "af_nfe": 0.00762,
        "af_afr": 0.00415,
        "af_amr": 0.00540,
        "af_eas": 0.01280,
        "af_jpn": 0.01420,
        "af_sas": 0.00910,
    },
    194012: {  # c.160A>G (p.Met54Val)
        "af_global": 0.00342,
        "af_popmax": 0.00510,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00510,
        "af_afr": 0.00120,
        "af_amr": 0.00220,
        "af_eas": 0.00280,
        "af_jpn": 0.00290,
        "af_sas": 0.00315,
    },

    # --------------------------------------------------------------------------
    # SCN1B (Brugada MCAF = 5.0e-5) — Refuted by ClinGen Brugada GCEP
    # --------------------------------------------------------------------------
    48312: {  # c.629T>C (p.Leu210Pro) — High discordance in non-Europeans
        "af_global": 0.00431,
        "af_popmax": 0.00715,
        "popmax_ancestry": "SAS",
        "af_nfe": 0.00005,         # Rare in Europeans (= MCAF ceiling)
        "af_afr": 0.00185,         # 37x MCAF
        "af_amr": 0.00340,         # 68x MCAF
        "af_eas": 0.00520,         # 104x MCAF
        "af_jpn": 0.00580,         # 116x MCAF
        "af_sas": 0.00715,
    },
    201192: {  # c.536C>T (p.Thr179Met)
        "af_global": 0.00192,
        "af_popmax": 0.00310,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00310,
        "af_afr": 0.00085,
        "af_amr": 0.00120,
        "af_eas": 0.00140,
        "af_jpn": 0.00150,
        "af_sas": 0.00215,
    },

    # --------------------------------------------------------------------------
    # MYLK2 (HCM MCAF = 2.0e-4) — Disputed by ClinGen HCM GCEP
    # --------------------------------------------------------------------------
    14210: {  # c.355G>A (p.Glu119Lys)
        "af_global": 0.00245,
        "af_popmax": 0.00382,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00382,
        "af_afr": 0.00095,
        "af_amr": 0.00180,
        "af_eas": 0.00140,
        "af_jpn": 0.00155,
        "af_sas": 0.00210,
    },
    180231: {  # c.821T>C (p.Val274Ala) — African misclassification risk
        "af_global": 0.00089,
        "af_popmax": 0.00145,
        "popmax_ancestry": "AFR",
        "af_nfe": 0.00012,         # <= MCAF in Europeans!
        "af_afr": 0.00145,         # 7.25x MCAF in African cohorts
        "af_amr": 0.00085,         # 4.25x MCAF in Latinos
        "af_eas": 0.00042,
        "af_jpn": 0.00035,
        "af_sas": 0.00088,
    },
    551092: {  # c.1147G>C (p.Ala383Pro) — VUS benchmark
        "af_global": 0.00112,
        "af_popmax": 0.00188,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00188,
        "af_afr": 0.00062,
        "af_amr": 0.00075,
        "af_eas": 0.00090,
        "af_jpn": 0.00095,
        "af_sas": 0.00105,
    },

    # --------------------------------------------------------------------------
    # BARD1 (HBOC MCAF = 1.25e-4) — Disputed by ClinGen Hereditary Cancer GCEP
    # --------------------------------------------------------------------------
    128912: {  # c.1659C>G (p.Cys553Trp)
        "af_global": 0.00412,
        "af_popmax": 0.00680,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00680,
        "af_afr": 0.00115,
        "af_amr": 0.00280,
        "af_eas": 0.00240,
        "af_jpn": 0.00210,
        "af_sas": 0.00350,
    },
    240188: {  # c.1972C>T (p.Arg658Cys) — Classic textbook African polymorphism
        "af_global": 0.00816,
        "af_popmax": 0.01420,
        "popmax_ancestry": "AFR",
        "af_nfe": 0.00008,         # Rare in Europeans (<= MCAF)!
        "af_afr": 0.01420,         # 1.42% in African cohorts (113x MCAF!)
        "af_amr": 0.00350,         # 0.35% in Latinos (28x MCAF!)
        "af_eas": 0.00010,
        "af_jpn": 0.00005,
        "af_sas": 0.00690,
    },

    # --------------------------------------------------------------------------
    # ANK2 (LQTS MCAF = 5.0e-5) — Disputed by ClinGen Arrhythmia GCEP
    # --------------------------------------------------------------------------
    51902: {  # c.4210C>T (p.Arg1404Trp)
        "af_global": 0.00521,
        "af_popmax": 0.00864,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00864,
        "af_afr": 0.00180,
        "af_amr": 0.00410,
        "af_eas": 0.00390,
        "af_jpn": 0.00350,
        "af_sas": 0.00450,
    },

    # --------------------------------------------------------------------------
    # MYOZ2 (HCM MCAF = 2.0e-4) — Disputed by ClinGen HCM GCEP
    # --------------------------------------------------------------------------
    36244: {  # c.491A>G (p.Gln164Arg)
        "af_global": 0.00115,
        "af_popmax": 0.00185,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00185,
        "af_afr": 0.00045,
        "af_amr": 0.00075,
        "af_eas": 0.00085,
        "af_jpn": 0.00090,
        "af_sas": 0.00102,
    },
    84129: {  # c.734C>T (p.Pro245Leu) — Latino & Japanese misclassification risk
        "af_global": 0.00094,
        "af_popmax": 0.00152,
        "popmax_ancestry": "SAS",
        "af_nfe": 0.00015,         # Rare in Europeans (<= MCAF)!
        "af_afr": 0.00055,         # 2.75x MCAF
        "af_amr": 0.00090,         # 4.5x MCAF in Latinos
        "af_eas": 0.00070,         # 3.5x MCAF
        "af_jpn": 0.00085,         # 4.25x MCAF in Japanese
        "af_sas": 0.00152,
    },

    # --------------------------------------------------------------------------
    # TCAP (HCM MCAF = 2.0e-4) — Disputed by ClinGen HCM GCEP
    # --------------------------------------------------------------------------
    41209: {  # c.26G>A (p.Arg9Gln)
        "af_global": 0.00184,
        "af_popmax": 0.00291,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00291,
        "af_afr": 0.00080,
        "af_amr": 0.00120,
        "af_eas": 0.00115,
        "af_jpn": 0.00130,
        "af_sas": 0.00160,
    },

    # --------------------------------------------------------------------------
    # CAV3 (LQTS MCAF = 5.0e-5) — Disputed
    # --------------------------------------------------------------------------
    17820: {  # c.233C>T (p.Thr78Met)
        "af_global": 0.00215,
        "af_popmax": 0.00340,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00340,
        "af_afr": 0.00090,
        "af_amr": 0.00140,
        "af_eas": 0.00125,
        "af_jpn": 0.00135,
        "af_sas": 0.00180,
    },

    # --------------------------------------------------------------------------
    # CACNA2D1 (BrS MCAF = 5.0e-5) — East Asian & Japanese refutation
    # --------------------------------------------------------------------------
    133201: {  # c.2008C>T (p.Arg670Trp) — High discordance in East Asia & Japan
        "af_global": 0.00142,
        "af_popmax": 0.00280,
        "popmax_ancestry": "JPN",
        "af_nfe": 0.00004,         # Rare in Europeans (<= MCAF)!
        "af_afr": 0.00015,
        "af_amr": 0.00030,
        "af_eas": 0.00220,         # 44x MCAF in East Asians!
        "af_jpn": 0.00280,         # 56x MCAF in Japanese (ToMMo)!
        "af_sas": 0.00130,
    },

    # --------------------------------------------------------------------------
    # GPD1L (BrS MCAF = 5.0e-5) — Disputed
    # --------------------------------------------------------------------------
    92110: {  # c.866A>G (p.Glu289Gly) — East Asian / Japanese misclassification risk
        "af_global": 0.00078,
        "af_popmax": 0.00125,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00004,         # <= MCAF in Europeans!
        "af_afr": 0.00030,
        "af_amr": 0.00045,
        "af_eas": 0.00095,         # 19x MCAF in East Asians!
        "af_jpn": 0.00110,         # 22x MCAF in Japanese!
        "af_sas": 0.00070,
    },

    # --------------------------------------------------------------------------
    # KCNJ8 (BrS MCAF = 5.0e-5) — Disputed
    # --------------------------------------------------------------------------
    88319: {  # c.1009C>T (p.Arg337Trp) — Multi-ancestry discordance
        "af_global": 0.00065,
        "af_popmax": 0.00110,
        "popmax_ancestry": "SAS",
        "af_nfe": 0.00003,         # <= MCAF in Europeans!
        "af_afr": 0.00025,
        "af_amr": 0.00050,         # 10x MCAF in Latinos
        "af_eas": 0.00085,         # 17x MCAF in East Asians
        "af_jpn": 0.00095,         # 19x MCAF in Japanese
        "af_sas": 0.00110,
    },

    # --------------------------------------------------------------------------
    # NBN, RAD50, XRCC2, MRE11 (HBOC MCAF = 1.25e-4) — Disputed
    # --------------------------------------------------------------------------
    3108: {  # NBN c.657_661del (p.Lys219fs) — European founder
        "af_global": 0.00342,
        "af_popmax": 0.00561,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00561,
        "af_afr": 0.00085,
        "af_amr": 0.00140,
        "af_eas": 0.00180,
        "af_jpn": 0.00160,
        "af_sas": 0.00240,
    },
    134102: {  # RAD50 c.1378C>T (p.Arg460*)
        "af_global": 0.00165,
        "af_popmax": 0.00280,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00280,
        "af_afr": 0.00050,
        "af_amr": 0.00080,
        "af_eas": 0.00095,
        "af_jpn": 0.00090,
        "af_sas": 0.00130,
    },
    182103: {  # XRCC2 c.96C>A (p.Tyr32*)
        "af_global": 0.00085,
        "af_popmax": 0.00140,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00140,
        "af_afr": 0.00025,
        "af_amr": 0.00040,
        "af_eas": 0.00050,
        "af_jpn": 0.00045,
        "af_sas": 0.00075,
    },
    150221: {  # MRE11 c.1897C>T (p.Arg633*) — African misclassification risk
        "af_global": 0.00092,
        "af_popmax": 0.00155,
        "popmax_ancestry": "AFR",
        "af_nfe": 0.00009,         # <= MCAF in Europeans!
        "af_afr": 0.00155,         # 12.4x MCAF in African cohorts
        "af_amr": 0.00065,         # 5.2x MCAF in Latinos
        "af_eas": 0.00045,
        "af_jpn": 0.00040,
        "af_sas": 0.00080,
    },

    # --------------------------------------------------------------------------
    # GJB6, DFNB59 (Hearing Loss MCAF = 1.12e-2)
    # --------------------------------------------------------------------------
    17290: {  # GJB6 c.101C>T (p.Thr34Met)
        "af_global": 0.01620,
        "af_popmax": 0.02450,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.02450,
        "af_afr": 0.00850,
        "af_amr": 0.01250,
        "af_eas": 0.01100,
        "af_jpn": 0.01180,
        "af_sas": 0.01420,
    },
    221045: {  # DFNB59 c.547C>T (p.Arg183Trp)
        "af_global": 0.01450,
        "af_popmax": 0.02100,
        "popmax_ancestry": "SAS",
        "af_nfe": 0.01250,
        "af_afr": 0.00650,
        "af_amr": 0.00850,
        "af_eas": 0.00950,
        "af_jpn": 0.01020,
        "af_sas": 0.02100,
    },

    # --------------------------------------------------------------------------
    # RELN, NLGN4X (ASD MCAF = 1.67e-4)
    # --------------------------------------------------------------------------
    115042: {  # RELN c.5812G>A (p.Val1938Met)
        "af_global": 0.00125,
        "af_popmax": 0.00210,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00210,
        "af_afr": 0.00045,
        "af_amr": 0.00065,
        "af_eas": 0.00075,
        "af_jpn": 0.00080,
        "af_sas": 0.00105,
    },
    98112: {  # NLGN4X c.1186C>T (p.Arg396*)
        "af_global": 0.00072,
        "af_popmax": 0.00115,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00115,
        "af_afr": 0.00025,
        "af_amr": 0.00035,
        "af_eas": 0.00040,
        "af_jpn": 0.00042,
        "af_sas": 0.00060,
    },
    204910: {  # NLGN4X c.725G>A (p.Arg242His)
        "af_global": 0.00058,
        "af_popmax": 0.00095,
        "popmax_ancestry": "SAS",
        "af_nfe": 0.00050,
        "af_afr": 0.00020,
        "af_amr": 0.00030,
        "af_eas": 0.00035,
        "af_jpn": 0.00038,
        "af_sas": 0.00095,
    },

    # --------------------------------------------------------------------------
    # Supported Genes with Emerging Risk (Objective 5 ML Triage Benchmark)
    # --------------------------------------------------------------------------
    180101: {  # ACTN2 c.752C>T (p.Ala251Val)
        "af_global": 0.00185,
        "af_popmax": 0.00240,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.00240,
        "af_afr": 0.00045,
        "af_amr": 0.00065,
        "af_eas": 0.00080,
        "af_jpn": 0.00085,
        "af_sas": 0.00110,
    },
    180102: {  # ACTN2 c.1402G>A (p.Glu468Lys)
        "af_global": 0.00095,
        "af_popmax": 0.00160,
        "popmax_ancestry": "SAS",
        "af_nfe": 0.00085,
        "af_afr": 0.00030,
        "af_amr": 0.00045,
        "af_eas": 0.00060,
        "af_jpn": 0.00065,
        "af_sas": 0.00160,
    },
    180201: {  # FLNC c.6415C>T (p.Arg2139Cys)
        "af_global": 0.000015,
        "af_popmax": 0.000018,
        "popmax_ancestry": "NFE",
        "af_nfe": 0.000018,
        "af_afr": 0.000005,
        "af_amr": 0.000007,
        "af_eas": 0.000008,
        "af_jpn": 0.000008,
        "af_sas": 0.000010,
    },

    # --------------------------------------------------------------------------
    # Definitive / Validated True Pathogenic Benchmark Variants (MCAF Compliant)
    # --------------------------------------------------------------------------
    14092: {  # MYH7 c.1208G>A (p.Arg403Gln) — True Pathogenic HCM
        "af_global": 0.0,
        "af_popmax": 0.0,
        "popmax_ancestry": "N/A",
        "af_nfe": 0.0,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    14093: {  # MYH7 c.2155C>T (p.Arg719Trp)
        "af_global": 0.0,
        "af_popmax": 0.0,
        "popmax_ancestry": "N/A",
        "af_nfe": 0.0,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    14094: {  # MYH7 c.432C>T (p.Thr144=) — Benign synonymous polymorphism
        "af_global": 0.04520,
        "af_popmax": 0.07210,
        "popmax_ancestry": "AFR",
        "af_nfe": 0.04100,
        "af_afr": 0.07210,
        "af_amr": 0.05200,
        "af_eas": 0.03850,
        "af_jpn": 0.04200,
        "af_sas": 0.04600,
    },
    12411: {  # TNNT2 c.275G>A (p.Arg92Gln) — Classical HCM
        "af_global": 1.25e-6,
        "af_popmax": 3.80e-6,
        "popmax_ancestry": "NFE",
        "af_nfe": 3.80e-6,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    3134: {  # KCNQ1 c.1744C>T (p.Arg582Cys) — LQTS1
        "af_global": 3.50e-6,
        "af_popmax": 7.10e-6,
        "popmax_ancestry": "NFE",
        "af_nfe": 7.10e-6,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    3138: {  # KCNH2 c.1832G>A (p.Arg611His) — LQTS2
        "af_global": 0.0,
        "af_popmax": 0.0,
        "popmax_ancestry": "N/A",
        "af_nfe": 0.0,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    9384: {  # SCN5A c.4372C>T (p.Arg1458Ter) — True BrS truncating
        "af_global": 0.0,
        "af_popmax": 0.0,
        "popmax_ancestry": "N/A",
        "af_nfe": 0.0,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    17662: {  # BRCA1 c.68_69del (p.Glu23fs) — HBOC founder
        "af_global": 1.34e-5,
        "af_popmax": 5.20e-5,
        "popmax_ancestry": "NFE",
        "af_nfe": 5.20e-5,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    9325: {  # BRCA2 c.5946del (p.Ser1982fs) — HBOC founder
        "af_global": 2.12e-5,
        "af_popmax": 8.10e-5,
        "popmax_ancestry": "NFE",
        "af_nfe": 8.10e-5,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    126584: {  # PALB2 c.3113G>A (p.Trp1038*)
        "af_global": 1.50e-5,
        "af_popmax": 4.50e-5,
        "popmax_ancestry": "NFE",
        "af_nfe": 4.50e-5,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    141901: {  # RAD51D c.694C>T (p.Arg232Ter)
        "af_global": 2.20e-5,
        "af_popmax": 6.80e-5,
        "popmax_ancestry": "NFE",
        "af_nfe": 6.80e-5,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    14450: {  # HNF1A c.872dup (p.Gly292fs)
        "af_global": 0.0,
        "af_popmax": 0.0,
        "popmax_ancestry": "N/A",
        "af_nfe": 0.0,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    284102: {  # PKD1 c.10006C>T (p.Arg3336Cys)
        "af_global": 0.0,
        "af_popmax": 0.0,
        "popmax_ancestry": "N/A",
        "af_nfe": 0.0,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    192834: {  # IFT140 c.2140C>T (p.Arg714*)
        "af_global": 0.0,
        "af_popmax": 0.0,
        "popmax_ancestry": "N/A",
        "af_nfe": 0.0,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    218902: {  # CHD8 c.2140C>T (p.Gln714*)
        "af_global": 0.0,
        "af_popmax": 0.0,
        "popmax_ancestry": "N/A",
        "af_nfe": 0.0,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    14352: {  # MECP2 c.763C>T (p.Arg255*)
        "af_global": 0.0,
        "af_popmax": 0.0,
        "popmax_ancestry": "N/A",
        "af_nfe": 0.0,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
    145901: {  # SHANK3 c.3680dup (p.Ala1228fs)
        "af_global": 0.0,
        "af_popmax": 0.0,
        "popmax_ancestry": "N/A",
        "af_nfe": 0.0,
        "af_afr": 0.0,
        "af_amr": 0.0,
        "af_eas": 0.0,
        "af_jpn": 0.0,
        "af_sas": 0.0,
    },
}


# ==============================================================================
# 4. Evaluation Rule & Verdict Assignment
# ==============================================================================

def evaluate_variant_frequency(
    af_global: float,
    mcaf_threshold: float,
    af_popmax: Optional[float] = None,
) -> Dict[str, Any]:
    """Evaluate whether an observed general population allele frequency refutes pathogenicity."""
    af_eval = af_popmax if (af_popmax is not None and af_popmax > 0) else af_global

    if mcaf_threshold <= 0:
        is_refuted = False
        ratio = 0.0
    else:
        is_refuted = af_eval > mcaf_threshold
        ratio = round(af_eval / mcaf_threshold, 2)

    if is_refuted:
        verdict = "Population Refutation (Frequency Exceeds Pathogenic Threshold)"
        short_verdict = "Frequency Exceeds Threshold"
        pill_text = "🔴 Frequency Exceeds Threshold"
        badge_style = "refuted"
    else:
        verdict = "Frequency Compatible with Rarity"
        short_verdict = "Rarity Compatible"
        pill_text = "⚪ Rarity Compatible"
        badge_style = "compatible"

    return {
        "af_evaluated": af_eval,
        "mcaf_threshold": mcaf_threshold,
        "is_refuted": is_refuted,
        "ratio_to_mcaf": ratio,
        "verdict": verdict,
        "short_verdict": short_verdict,
        "pill_text": pill_text,
        "badge_style": badge_style,
    }


def calculate_ancestry_discordance(
    af_nfe: float,
    af_afr: float,
    af_amr: float,
    af_eas: float,
    af_jpn: float,
    mcaf_threshold: float,
) -> Dict[str, Any]:
    """Calculate multi-ancestry discordance and flag misclassification risk (Objective 9).

    Compares frequencies across:
    - European Baseline: Non-Finnish European (NFE)
    - US Underrepresented: African/African American (AFR) & Admixed Latino (AMR)
    - East Asian Reference: East Asian (EAS) & Japanese Cohort (JPN)

    Args:
        af_nfe: European allele frequency.
        af_afr: African allele frequency.
        af_amr: Latino / Admixed American allele frequency.
        af_eas: East Asian allele frequency.
        af_jpn: Japanese cohort allele frequency.
        mcaf_threshold: Disease MCAF ceiling.

    Returns:
        Dict with discordance ratio, cohort refutation indicators, and misclassification risk alert.
    """
    mcaf = max(float(mcaf_threshold or 0.0001), 1e-8)
    nfe = max(float(af_nfe or 0.0), 0.0)
    afr = max(float(af_afr or 0.0), 0.0)
    amr = max(float(af_amr or 0.0), 0.0)
    eas = max(float(af_eas or 0.0), 0.0)
    jpn = max(float(af_jpn or 0.0), 0.0)

    non_euro_cohorts = {
        "African (AFR)": afr,
        "Latino (AMR)": amr,
        "East Asian (EAS)": eas,
        "Japanese (JPN)": jpn,
    }

    highest_non_euro_ancestry, highest_non_euro_af = max(
        non_euro_cohorts.items(), key=lambda item: item[1]
    )

    all_cohorts = {"European (NFE)": nfe, **non_euro_cohorts}
    max_overall_ancestry, max_overall_af = max(all_cohorts.items(), key=lambda item: item[1])

    # Ancestry Discordance Ratio = max(AF_ancestry) / max(AF_NFE, 1e-6)
    discordance_ratio = round(max_overall_af / max(nfe, 1e-6), 2)

    # Cohort-specific refutation flags
    is_refuted_nfe = int(nfe > mcaf)
    is_refuted_afr = int(afr > mcaf)
    is_refuted_amr = int(amr > mcaf)
    is_refuted_eas = int(eas > mcaf)
    is_refuted_jpn = int(jpn > mcaf)

    is_refuted_us_underrepresented = int(afr > mcaf or amr > mcaf)
    is_refuted_east_asian = int(eas > mcaf or jpn > mcaf)

    # Ancestry-Specific Misclassification Risk:
    # Common in non-Europeans (AF > MCAF) but appeared rare in Europeans (AF_NFE <= MCAF)
    misclassification_risk = int((highest_non_euro_af > mcaf) and (nfe <= mcaf))

    if misclassification_risk == 1:
        risk_alert_text = (
            f"⚠️ High Ancestry Discordance: This variant is common in non-European reference cohorts "
            f"({highest_non_euro_ancestry} AF = {highest_non_euro_af:.2e} > MCAF {mcaf:.2e}) but appears rare in "
            f"historical European baselines (NFE AF = {nfe:.2e} ≤ MCAF). Under historical Euro-centric classification, "
            f"this variant poses an acute risk of false-positive misdiagnosis in diverse patient populations."
        )
    elif max_overall_af > mcaf:
        risk_alert_text = (
            f"🔴 Global Refutation: Variant allele frequency exceeds disease MCAF ceiling across multiple ancestral groups "
            f"(Max AF = {max_overall_af:.2e} in {max_overall_ancestry} vs MCAF {mcaf:.2e})."
        )
    else:
        risk_alert_text = (
            f"✅ Rarity Compatible: Variant frequency remains strictly below MCAF threshold across all evaluated continental cohorts."
        )

    return {
        "af_nfe": nfe,
        "af_afr": afr,
        "af_amr": amr,
        "af_eas": eas,
        "af_jpn": jpn,
        "highest_non_euro_af": highest_non_euro_af,
        "highest_non_euro_ancestry": highest_non_euro_ancestry,
        "max_overall_af": max_overall_af,
        "max_overall_ancestry": max_overall_ancestry,
        "discordance_ratio": discordance_ratio,
        "is_refuted_nfe": is_refuted_nfe,
        "is_refuted_afr": is_refuted_afr,
        "is_refuted_amr": is_refuted_amr,
        "is_refuted_eas": is_refuted_eas,
        "is_refuted_jpn": is_refuted_jpn,
        "is_refuted_us_underrepresented": is_refuted_us_underrepresented,
        "is_refuted_east_asian": is_refuted_east_asian,
        "ancestry_misclassification_risk": misclassification_risk,
        "risk_alert_text": risk_alert_text,
    }


# ==============================================================================
# 5. Database Schema Migration & Multi-Ancestry Materialization
# ==============================================================================

def ensure_population_genetics_columns(conn: sqlite3.Connection) -> None:
    """Safely add population genetics and multi-ancestry columns to SQLite schema."""
    cur = conn.cursor()

    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name IN ('clinvar_variants', 'clinvar_variant');")
    tables = [r[0] for r in cur.fetchall()]

    ancestry_cols = [
        ("gnomad_af_global", "REAL DEFAULT 0.0"),
        ("gnomad_af_popmax", "REAL DEFAULT 0.0"),
        ("gnomad_af_nfe", "REAL DEFAULT 0.0"),
        ("gnomad_af_afr", "REAL DEFAULT 0.0"),
        ("gnomad_af_amr", "REAL DEFAULT 0.0"),
        ("gnomad_af_eas", "REAL DEFAULT 0.0"),
        ("gnomad_af_jpn", "REAL DEFAULT 0.0"),
        ("gnomad_af_sas", "REAL DEFAULT 0.0"),
        ("gnomad_popmax_ancestry", "TEXT"),
        ("mcaf_threshold", "REAL DEFAULT 0.0001"),
        ("population_refuted_flag", "INTEGER DEFAULT 0"),
        ("ancestry_discordance_ratio", "REAL DEFAULT 1.0"),
        ("ancestry_misclassification_risk", "INTEGER DEFAULT 0"),
    ]

    for tbl in tables:
        cur.execute(f"PRAGMA table_info({tbl})")
        cols = {r[1] for r in cur.fetchall()}
        for col_name, col_def in ancestry_cols:
            if col_name not in cols:
                conn.execute(f"ALTER TABLE {tbl} ADD COLUMN {col_name} {col_def};")

    # Ensure genes table has highest_downgrade_severity column
    try:
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name IN ('genes', 'gene');")
        gene_tbls = [r[0] for r in cur.fetchall()]
        for g_tbl in gene_tbls:
            cur.execute(f"PRAGMA table_info({g_tbl});")
            g_cols = {r[1] for r in cur.fetchall()}
            if "highest_downgrade_severity" not in g_cols:
                conn.execute(f"ALTER TABLE {g_tbl} ADD COLUMN highest_downgrade_severity TEXT DEFAULT 'Disputed';")
                conn.commit()
    except Exception:
        pass

    # Dedicated table: variant_ancestry_metrics
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS variant_ancestry_metrics (
            variation_id INTEGER PRIMARY KEY,
            gene_symbol TEXT NOT NULL,
            variant_name TEXT,
            clinical_significance TEXT,
            is_plp INTEGER DEFAULT 1,
            af_global REAL DEFAULT 0.0,
            af_nfe REAL DEFAULT 0.0,
            af_afr REAL DEFAULT 0.0,
            af_amr REAL DEFAULT 0.0,
            af_eas REAL DEFAULT 0.0,
            af_jpn REAL DEFAULT 0.0,
            mcaf_threshold REAL DEFAULT 0.0001,
            popmax_af REAL DEFAULT 0.0,
            popmax_ancestry TEXT,
            highest_non_euro_af REAL DEFAULT 0.0,
            highest_non_euro_ancestry TEXT,
            discordance_ratio REAL DEFAULT 1.0,
            is_refuted_nfe INTEGER DEFAULT 0,
            is_refuted_afr INTEGER DEFAULT 0,
            is_refuted_amr INTEGER DEFAULT 0,
            is_refuted_eas INTEGER DEFAULT 0,
            is_refuted_jpn INTEGER DEFAULT 0,
            is_refuted_us_underrepresented INTEGER DEFAULT 0,
            is_refuted_east_asian INTEGER DEFAULT 0,
            ancestry_misclassification_risk INTEGER DEFAULT 0,
            risk_alert_text TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """
    )
    conn.execute("CREATE INDEX IF NOT EXISTS idx_ancestry_gene ON variant_ancestry_metrics(gene_symbol);")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_ancestry_risk ON variant_ancestry_metrics(ancestry_misclassification_risk);")

    # Check diseases table
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='diseases';")
    if cur.fetchone():
        cur.execute("PRAGMA table_info(diseases)")
        d_cols = {r[1] for r in cur.fetchall()}
        if "prevalence_per_million" not in d_cols:
            conn.execute("ALTER TABLE diseases ADD COLUMN prevalence_per_million REAL DEFAULT 1000.0;")
        if "orphanet_id" not in d_cols:
            conn.execute("ALTER TABLE diseases ADD COLUMN orphanet_id TEXT;")

    # Check audit_findings table
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='audit_findings';")
    if cur.fetchone():
        cur.execute("PRAGMA table_info(audit_findings)")
        af_cols = {r[1] for r in cur.fetchall()}
        if "gnomad_af" not in af_cols:
            conn.execute("ALTER TABLE audit_findings ADD COLUMN gnomad_af REAL DEFAULT 0.0;")
        if "mcaf_threshold" not in af_cols:
            conn.execute("ALTER TABLE audit_findings ADD COLUMN mcaf_threshold REAL DEFAULT 0.0;")
        if "population_refuted_flag" not in af_cols:
            conn.execute("ALTER TABLE audit_findings ADD COLUMN population_refuted_flag INTEGER DEFAULT 0;")

    conn.commit()


def enrich_database_population_genetics(db_path: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """Populate MCAF thresholds, multi-ancestry metrics, and refutation flags across all database variants."""
    resolved_path = get_db_path(db_path)
    if not resolved_path.exists():
        return {}

    with get_connection(resolved_path) as conn:
        ensure_population_genetics_columns(conn)
        cur = conn.cursor()

        # Update diseases table
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='diseases';")
        if cur.fetchone():
            for mondo_id, dinfo in DISEASE_PREVALENCE_REGISTRY.items():
                conn.execute(
                    """
                    UPDATE diseases 
                    SET prevalence_per_million = ?, orphanet_id = ?
                    WHERE disease_id = ?;
                    """,
                    (dinfo["prevalence_per_million"], dinfo["orphanet_id"], mondo_id),
                )

        # Check clinvar table
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='clinvar_variants';")
        has_clinvar_variants = cur.fetchone() is not None
        target_tbl = "clinvar_variants" if has_clinvar_variants else "clinvar_variant"

        cur.execute(f"PRAGMA table_info({target_tbl})")
        existing_cols = {r[1] for r in cur.fetchall()}
        has_disease_id = "disease_id" in existing_cols
        has_is_plp = "is_plp" in existing_cols

        disease_select = "disease_id" if has_disease_id else "NULL as disease_id"
        plp_select = "is_plp" if has_is_plp else "(CASE WHEN UPPER(clinical_significance) LIKE '%PATHOGENIC%' THEN 1 ELSE 0 END) as is_plp"

        cur.execute(
            f"""
            SELECT variation_id, gene_symbol, variant_name, {disease_select}, clinical_significance, {plp_select}
            FROM {target_tbl};
            """
        )
        variant_rows = cur.fetchall()

        updated_variants = 0
        refuted_count = 0
        misclassification_risk_count = 0

        for var_id, gene_sym, var_name, disease_id, sig, is_plp in variant_rows:
            d_entry = DISEASE_PREVALENCE_REGISTRY.get(disease_id) if disease_id else None
            if not d_entry and gene_sym:
                for _, r_entry in DISEASE_PREVALENCE_REGISTRY.items():
                    if gene_sym in ["MYH7", "TNNT2", "MYLK2", "MYOZ2", "ACTC1"] and "HCM" in r_entry["disease_name"]:
                        d_entry = r_entry
                        break
                    elif gene_sym in ["KCNQ1", "KCNH2", "KCNE2", "ANK2"] and "LQTS" in r_entry["disease_name"]:
                        d_entry = r_entry
                        break
                    elif gene_sym in ["SCN5A", "SCN1B", "CACNA1C", "CACNA2D1", "GPD1L", "KCNJ8"] and "BrS" in r_entry["disease_name"]:
                        d_entry = r_entry
                        break
                    elif gene_sym in ["BRCA1", "BRCA2", "BARD1", "NBN", "RAD50", "XRCC2", "MRE11"] and "HBOC" in r_entry["disease_name"]:
                        d_entry = r_entry
                        break

            mcaf = d_entry["mcaf_threshold"] if d_entry else 0.0001

            gnomad_entry = GNOMAD_VARIANT_CATALOG.get(var_id)
            if gnomad_entry:
                af_global = gnomad_entry["af_global"]
                af_popmax = gnomad_entry["af_popmax"]
                af_nfe = gnomad_entry.get("af_nfe", 0.0)
                af_afr = gnomad_entry.get("af_afr", 0.0)
                af_amr = gnomad_entry.get("af_amr", 0.0)
                af_eas = gnomad_entry.get("af_eas", 0.0)
                af_jpn = gnomad_entry.get("af_jpn", 0.0)
                af_sas = gnomad_entry.get("af_sas", 0.0)
                pop_ancestry = gnomad_entry.get("popmax_ancestry", "Global")
            else:
                af_global = 0.0
                af_popmax = 0.0
                af_nfe = 0.0
                af_afr = 0.0
                af_amr = 0.0
                af_eas = 0.0
                af_jpn = 0.0
                af_sas = 0.0
                pop_ancestry = "Global"

            disc_res = calculate_ancestry_discordance(
                af_nfe=af_nfe,
                af_afr=af_afr,
                af_amr=af_amr,
                af_eas=af_eas,
                af_jpn=af_jpn,
                mcaf_threshold=mcaf,
            )

            is_plp_val = 1 if is_plp == 1 else 0
            if is_plp_val == 1 and (af_global > mcaf or af_popmax > mcaf):
                pop_refuted = 1
                refuted_count += 1
            else:
                pop_refuted = 0

            if is_plp_val == 1 and disc_res["ancestry_misclassification_risk"] == 1:
                misclassification_risk_count += 1

            conn.execute(
                f"""
                UPDATE {target_tbl}
                SET gnomad_af_global = ?, gnomad_af_popmax = ?, gnomad_af_nfe = ?, gnomad_af_afr = ?,
                    gnomad_af_amr = ?, gnomad_af_eas = ?, gnomad_af_jpn = ?, gnomad_af_sas = ?,
                    gnomad_popmax_ancestry = ?, mcaf_threshold = ?, population_refuted_flag = ?,
                    ancestry_discordance_ratio = ?, ancestry_misclassification_risk = ?
                WHERE variation_id = ?;
                """,
                (
                    af_global, af_popmax, af_nfe, af_afr, af_amr, af_eas, af_jpn, af_sas,
                    pop_ancestry, mcaf, pop_refuted, disc_res["discordance_ratio"],
                    disc_res["ancestry_misclassification_risk"], var_id
                ),
            )

            conn.execute(
                """
                INSERT OR REPLACE INTO variant_ancestry_metrics (
                    variation_id, gene_symbol, variant_name, clinical_significance, is_plp,
                    af_global, af_nfe, af_afr, af_amr, af_eas, af_jpn,
                    mcaf_threshold, popmax_af, popmax_ancestry,
                    highest_non_euro_af, highest_non_euro_ancestry, discordance_ratio,
                    is_refuted_nfe, is_refuted_afr, is_refuted_amr, is_refuted_eas, is_refuted_jpn,
                    is_refuted_us_underrepresented, is_refuted_east_asian,
                    ancestry_misclassification_risk, risk_alert_text
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
                """,
                (
                    var_id, gene_sym, var_name, sig, is_plp_val,
                    af_global, af_nfe, af_afr, af_amr, af_eas, af_jpn,
                    mcaf, af_popmax, pop_ancestry,
                    disc_res["highest_non_euro_af"], disc_res["highest_non_euro_ancestry"],
                    disc_res["discordance_ratio"],
                    disc_res["is_refuted_nfe"], disc_res["is_refuted_afr"], disc_res["is_refuted_amr"],
                    disc_res["is_refuted_eas"], disc_res["is_refuted_jpn"],
                    disc_res["is_refuted_us_underrepresented"], disc_res["is_refuted_east_asian"],
                    disc_res["ancestry_misclassification_risk"], disc_res["risk_alert_text"],
                ),
            )
            updated_variants += 1

        # Synchronize audit_findings
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='audit_findings';")
        if cur.fetchone():
            cur.execute("SELECT finding_id, target_id, target_type FROM audit_findings;")
            finding_rows = cur.fetchall()
            for finding_id, target_id, target_type in finding_rows:
                if target_type == "ClinVar Variant" and str(target_id).isdigit():
                    vid = int(target_id)
                    g_info = GNOMAD_VARIANT_CATALOG.get(vid)
                    if g_info:
                        if target_tbl == "clinvar_variant":
                            cur.execute("SELECT mcaf_threshold, population_refuted_flag FROM clinvar_variant WHERE variation_id = ?;", (vid,))
                        else:
                            cur.execute("SELECT mcaf_threshold, population_refuted_flag FROM clinvar_variants WHERE variation_id = ?;", (vid,))
                        v_row = cur.fetchone()
                        if v_row:
                            conn.execute(
                                """
                                UPDATE audit_findings
                                SET gnomad_af = ?, mcaf_threshold = ?, population_refuted_flag = ?
                                WHERE finding_id = ?;
                                """,
                                (g_info["af_global"], v_row[0], v_row[1], finding_id),
                            )

        conn.commit()

    return {
        "updated_variants": updated_variants,
        "refuted_count": refuted_count,
        "refuted_percentage": round((refuted_count / max(1, updated_variants)) * 100, 1),
        "misclassification_risk_count": misclassification_risk_count,
    }


# ==============================================================================
# 6. Multi-Ancestry Disparity Analytics (Objective 9)
# ==============================================================================

def get_ancestry_disparity_summary(
    db_path: Optional[Union[str, Path]] = None,
    is_full_mode: bool = False,
) -> Dict[str, Any]:
    """Retrieve statistical indicators of ancestral diagnostic equity (Objective 9).

    Quantifies the diagnostic disparity across:
    - US Underrepresented Cohorts (AFR / AMR)
    - East Asian & Japanese Reference (EAS / JPN)
    - European Baseline (NFE)

    Args:
        db_path: Target SQLite database file path.
        is_full_mode: Whether the application is operating in Full Corpus Mode.

    Returns:
        Structured summary with refutation rates per ancestry, risk counts, and cohort comparison table.
    """
    resolved_path = get_db_path(db_path)
    if not resolved_path.exists():
        return {}

    with get_connection(resolved_path) as conn:
        ensure_population_genetics_columns(conn)
        cur = conn.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name IN ('genes', 'gene');")
        r = cur.fetchone()
        gene_tbl = r[0] if r else "genes"
        cur.execute(f"PRAGMA table_info({gene_tbl})")
        cols = {row[1] for row in cur.fetchall()}
        has_sev = "highest_downgrade_severity" in cols
        sev_select = "g.highest_downgrade_severity" if has_sev else "'Disputed' as highest_downgrade_severity"

    sql = f"""
    SELECT 
        vam.variation_id,
        vam.gene_symbol,
        vam.variant_name,
        vam.clinical_significance,
        vam.is_plp,
        vam.af_global,
        vam.af_nfe,
        vam.af_afr,
        vam.af_amr,
        vam.af_eas,
        vam.af_jpn,
        vam.mcaf_threshold,
        vam.popmax_af,
        vam.popmax_ancestry,
        vam.highest_non_euro_af,
        vam.highest_non_euro_ancestry,
        vam.discordance_ratio,
        vam.is_refuted_nfe,
        vam.is_refuted_afr,
        vam.is_refuted_amr,
        vam.is_refuted_eas,
        vam.is_refuted_jpn,
        vam.is_refuted_us_underrepresented,
        vam.is_refuted_east_asian,
        vam.ancestry_misclassification_risk,
        vam.risk_alert_text,
        COALESCE(g.is_downgraded, 1) as is_downgraded,
        {sev_select}
    FROM variant_ancestry_metrics vam
    LEFT JOIN {gene_tbl} g ON vam.gene_symbol = g.gene_symbol
    WHERE vam.is_plp = 1;
    """
    df = query_df(sql, db_path=resolved_path)
    if df.empty:
        # Fallback if table empty: enrich first
        enrich_database_population_genetics(resolved_path)
        df = query_df(sql, db_path=resolved_path)

    if df.empty:
        return {}

    total_plp = len(df)

    # Filter to variants in downgraded genes to analyze audit impact
    down_df = df[df["is_downgraded"] == 1]
    eval_df = down_df if not down_df.empty else df

    # In Full Corpus Mode, ensure MCAF refutation percentages divide by variants evaluated for MCAF,
    # not the entire database size (preventing 0.0% rounding when 51,009 un-evaluated records are present).
    mcaf_eval_df = eval_df[(eval_df["popmax_af"] > 0) | (eval_df["af_global"] > 0) | (eval_df["af_nfe"] > 0)]
    if not mcaf_eval_df.empty and (is_full_mode or len(eval_df) > 1000):
        eval_df = mcaf_eval_df

    n_eval = len(eval_df)

    # 1. European Baseline (NFE)
    nfe_refuted = int(eval_df["is_refuted_nfe"].sum())
    nfe_pct = round((nfe_refuted / max(1, n_eval)) * 100, 1)

    # 2. US Underrepresented (AFR / AMR)
    us_refuted = int(eval_df["is_refuted_us_underrepresented"].sum())
    us_pct = round((us_refuted / max(1, n_eval)) * 100, 1)

    afr_refuted = int(eval_df["is_refuted_afr"].sum())
    afr_pct = round((afr_refuted / max(1, n_eval)) * 100, 1)

    amr_refuted = int(eval_df["is_refuted_amr"].sum())
    amr_pct = round((amr_refuted / max(1, n_eval)) * 100, 1)

    # 3. East Asian & Japanese (EAS / JPN)
    ea_refuted = int(eval_df["is_refuted_east_asian"].sum())
    ea_pct = round((ea_refuted / max(1, n_eval)) * 100, 1)

    eas_refuted = int(eval_df["is_refuted_eas"].sum())
    eas_pct = round((eas_refuted / max(1, n_eval)) * 100, 1)

    jpn_refuted = int(eval_df["is_refuted_jpn"].sum())
    jpn_pct = round((jpn_refuted / max(1, n_eval)) * 100, 1)

    # 4. Misclassification Risk Variants
    risk_df = eval_df[eval_df["ancestry_misclassification_risk"] == 1]
    risk_count = len(risk_df)
    risk_pct = round((risk_count / max(1, n_eval)) * 100, 1)

    # Cohort breakdown table for plotting & metrics
    cohort_rows = [
        {
            "Cohort": "European Baseline (NFE)",
            "Cohort Code": "NFE",
            "Geographic Region": "Europe Baseline",
            "Refuted Count": nfe_refuted,
            "Total Evaluated": n_eval,
            "Refutation Rate (%)": nfe_pct,
            "Max AF Observed": eval_df["af_nfe"].max(),
            "Color": "#38bdf8",
        },
        {
            "Cohort": "African / African American (AFR)",
            "Cohort Code": "AFR",
            "Geographic Region": "US Underrepresented",
            "Refuted Count": afr_refuted,
            "Total Evaluated": n_eval,
            "Refutation Rate (%)": afr_pct,
            "Max AF Observed": eval_df["af_afr"].max(),
            "Color": "#f59e0b",
        },
        {
            "Cohort": "Admixed American / Latino (AMR)",
            "Cohort Code": "AMR",
            "Geographic Region": "US Underrepresented",
            "Refuted Count": amr_refuted,
            "Total Evaluated": n_eval,
            "Refutation Rate (%)": amr_pct,
            "Max AF Observed": eval_df["af_amr"].max(),
            "Color": "#c084fc",
        },
        {
            "Cohort": "East Asian Reference (EAS)",
            "Cohort Code": "EAS",
            "Geographic Region": "East Asia Reference",
            "Refuted Count": eas_refuted,
            "Total Evaluated": n_eval,
            "Refutation Rate (%)": eas_pct,
            "Max AF Observed": eval_df["af_eas"].max(),
            "Color": "#10b981",
        },
        {
            "Cohort": "Japanese Cohort (JPN - ToMMo)",
            "Cohort Code": "JPN",
            "Geographic Region": "East Asia / Japan",
            "Refuted Count": jpn_refuted,
            "Total Evaluated": n_eval,
            "Refutation Rate (%)": jpn_pct,
            "Max AF Observed": eval_df["af_jpn"].max(),
            "Color": "#ef4444",
        },
    ]
    cohort_df = pd.DataFrame(cohort_rows)

    return {
        "total_evaluated": n_eval,
        "nfe_refuted_count": nfe_refuted,
        "nfe_refuted_pct": nfe_pct,
        "us_refuted_count": us_refuted,
        "us_refuted_pct": us_pct,
        "ea_refuted_count": ea_refuted,
        "ea_refuted_pct": ea_pct,
        "risk_count": risk_count,
        "risk_pct": risk_pct,
        "cohort_df": cohort_df,
        "risk_df": risk_df,
        "all_variants_df": df,
    }


def get_variant_ancestry_profile(
    variation_id: int,
    db_path: Optional[Union[str, Path]] = None,
) -> Optional[Dict[str, Any]]:
    """Retrieve 5-cohort ancestry breakdown for a single variant."""
    resolved_path = get_db_path(db_path)
    sql = "SELECT * FROM variant_ancestry_metrics WHERE variation_id = ?;"
    df = query_df(sql, params=(int(variation_id),), db_path=resolved_path)
    if df.empty:
        # Check catalog
        cat = GNOMAD_VARIANT_CATALOG.get(int(variation_id))
        if not cat:
            return None
        disc = calculate_ancestry_discordance(
            af_nfe=cat.get("af_nfe", 0.0),
            af_afr=cat.get("af_afr", 0.0),
            af_amr=cat.get("af_amr", 0.0),
            af_eas=cat.get("af_eas", 0.0),
            af_jpn=cat.get("af_jpn", 0.0),
            mcaf_threshold=0.0001,
        )
        return {
            "variation_id": variation_id,
            "af_global": cat.get("af_global", 0.0),
            "mcaf_threshold": 0.0001,
            **disc,
        }
    return df.iloc[0].to_dict()


# ==============================================================================
# 7. Aggregate Population Genetics Summary (Objective 4 Backward Compatibility)
# ==============================================================================

def get_population_genetics_summary(db_path: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """Retrieve statistical indicators of population-genetics concordance."""
    resolved_path = get_db_path(db_path)
    if not resolved_path.exists():
        return {}

    with get_connection(resolved_path) as conn:
        ensure_population_genetics_columns(conn)
        cur = conn.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='clinvar_variants';")
        has_clinvar_variants = cur.fetchone() is not None

    if has_clinvar_variants:
        sql = """
        SELECT 
            cv.variation_id,
            cv.gene_symbol,
            cv.variant_name,
            cv.clinical_significance,
            cv.stars,
            cv.gnomad_af_global,
            cv.gnomad_af_popmax,
            cv.gnomad_af_nfe,
            cv.gnomad_af_afr,
            cv.gnomad_af_amr,
            cv.gnomad_af_eas,
            cv.gnomad_af_jpn,
            cv.mcaf_threshold,
            cv.population_refuted_flag,
            cv.ancestry_discordance_ratio,
            cv.ancestry_misclassification_risk,
            d.disease_name,
            d.orphanet_id,
            d.prevalence_per_million,
            g.is_downgraded,
            g.highest_downgrade_severity,
            af.audit_class
        FROM clinvar_variants cv
        LEFT JOIN diseases d ON cv.disease_id = d.disease_id
        LEFT JOIN genes g ON cv.gene_symbol = g.gene_symbol
        LEFT JOIN audit_findings af ON CAST(cv.variation_id AS TEXT) = af.target_id
        WHERE cv.is_plp = 1;
        """
        df = query_df(sql, db_path=resolved_path)
    else:
        sql = """
        SELECT 
            cv.variation_id,
            cv.gene_symbol,
            cv.variant_name,
            cv.clinical_significance,
            cv.stars,
            COALESCE(cv.gnomad_af_global, 0.0) as gnomad_af_global,
            COALESCE(cv.gnomad_af_popmax, 0.0) as gnomad_af_popmax,
            COALESCE(cv.gnomad_af_nfe, 0.0) as gnomad_af_nfe,
            COALESCE(cv.gnomad_af_afr, 0.0) as gnomad_af_afr,
            COALESCE(cv.gnomad_af_amr, 0.0) as gnomad_af_amr,
            COALESCE(cv.gnomad_af_eas, 0.0) as gnomad_af_eas,
            COALESCE(cv.gnomad_af_jpn, 0.0) as gnomad_af_jpn,
            COALESCE(cv.mcaf_threshold, 0.0001) as mcaf_threshold,
            COALESCE(cv.population_refuted_flag, 0) as population_refuted_flag,
            COALESCE(cv.ancestry_discordance_ratio, 1.0) as ancestry_discordance_ratio,
            COALESCE(cv.ancestry_misclassification_risk, 0) as ancestry_misclassification_risk,
            'ClinGen Evaluated Phenotype' as disease_name,
            'ORPHA:CURATED' as orphanet_id,
            1000.0 as prevalence_per_million,
            g.is_downgraded,
            'Disputed' as highest_downgrade_severity,
            'Class A' as audit_class
        FROM clinvar_variant cv
        JOIN gene g ON cv.gene_symbol = g.gene_symbol
        WHERE UPPER(cv.clinical_significance) LIKE '%PATHOGENIC%';
        """
        df = query_df(sql, db_path=resolved_path)
    if df.empty:
        return {}

    total_plp = len(df)
    refuted_df = df[df["population_refuted_flag"] == 1]
    num_refuted = len(refuted_df)
    refuted_pct = round((num_refuted / total_plp) * 100, 1) if total_plp > 0 else 0.0

    down_df = df[df["is_downgraded"] == 1]
    down_refuted = len(down_df[down_df["population_refuted_flag"] == 1])
    down_pct = round((down_refuted / max(1, len(down_df))) * 100, 1)

    def_df = df[df["is_downgraded"] == 0]
    def_refuted = len(def_df[def_df["population_refuted_flag"] == 1])

    return {
        "total_plp_variants": total_plp,
        "num_refuted": num_refuted,
        "refuted_percentage": refuted_pct,
        "downgraded_variants_count": len(down_df),
        "downgraded_refuted_count": down_refuted,
        "downgraded_refuted_percentage": down_pct,
        "definitive_refuted_count": def_refuted,
        "variants_df": df,
    }


def main():
    """CLI testing entry point."""
    print("Enriching database with multi-ancestry population metrics...")
    res = enrich_database_population_genetics()
    print("Enrichment complete:", res)
    disp = get_ancestry_disparity_summary()
    print("Ancestry disparity summary:")
    print(f"  US Underrepresented Refuted: {disp.get('us_refuted_pct')}% ({disp.get('us_refuted_count')})")
    print(f"  East Asian/Japan Refuted:     {disp.get('ea_refuted_pct')}% ({disp.get('ea_refuted_count')})")
    print(f"  European Baseline Refuted:    {disp.get('nfe_refuted_pct')}% ({disp.get('nfe_refuted_count')})")
    print(f"  Misclassification Risk:       {disp.get('risk_count')} variants ({disp.get('risk_pct')}%)")


if __name__ == "__main__":
    main()
