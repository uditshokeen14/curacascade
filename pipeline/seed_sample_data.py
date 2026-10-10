"""CuraCascade — Realistic Mock & Sample Data Generator.

Generates a realistic, biologically accurate seed corpus of 35+ curated genes,
including re-appraised cardiovascular genes (MYH7, TNNT2, MYLK2, MYOZ2, KCNE2, SCN1B),
autism/neurodevelopmental genes (SHANK3, CHD8, NLGN4X, RELN, MECP2),
cancer predisposition genes (BRCA1, BRCA2, NBN, RAD50, XRCC2, BARD1),
and metabolic/renal disease genes.

Generates corresponding ClinGen curations, ClinVar variant assertions (P/LP, VUS, Benign),
and Genomics England PanelApp panel entries to enable immediate local audit execution.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

# Ensure repository root is on sys.path
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import config as C
from database.connection import (
    ensure_database_ready,
    execute_many,
)
from pipeline.audit import run_audit
from pipeline.allele_frequency import enrich_database_population_genetics
from pipeline.forecasting import run_forecasting_pipeline

# ==============================================================================
# Domain Entities for Realistic Seed Generation
# ==============================================================================

SAMPLE_DISEASES: List[Dict[str, str]] = [
    {
        "disease_id": "MONDO:0005045",
        "disease_name": "Hypertrophic Cardiomyopathy (HCM)",
        "mondo_id": "MONDO:0005045",
        "omim_id": "192600",
        "inheritance_mode": "Autosomal dominant",
    },
    {
        "disease_id": "MONDO:0005014",
        "disease_name": "Dilated Cardiomyopathy (DCM)",
        "mondo_id": "MONDO:0005014",
        "omim_id": "115200",
        "inheritance_mode": "Autosomal dominant",
    },
    {
        "disease_id": "MONDO:0001625",
        "disease_name": "Long QT Syndrome (LQTS)",
        "mondo_id": "MONDO:0001625",
        "omim_id": "192500",
        "inheritance_mode": "Autosomal dominant",
    },
    {
        "disease_id": "MONDO:0015263",
        "disease_name": "Brugada Syndrome (BrS)",
        "mondo_id": "MONDO:0015263",
        "omim_id": "601144",
        "inheritance_mode": "Autosomal dominant",
    },
    {
        "disease_id": "MONDO:0005260",
        "disease_name": "Autism Spectrum Disorder (ASD)",
        "mondo_id": "MONDO:0005260",
        "omim_id": "209850",
        "inheritance_mode": "Complex / Multifactorial",
    },
    {
        "disease_id": "MONDO:0011448",
        "disease_name": "Hereditary Breast and Ovarian Cancer Syndrome (HBOC)",
        "mondo_id": "MONDO:0011448",
        "omim_id": "604370",
        "inheritance_mode": "Autosomal dominant",
    },
    {
        "disease_id": "MONDO:0007817",
        "disease_name": "Autosomal Dominant Nonsyndromic Hearing Loss",
        "mondo_id": "MONDO:0007817",
        "omim_id": "600791",
        "inheritance_mode": "Autosomal dominant",
    },
    {
        "disease_id": "MONDO:0008304",
        "disease_name": "Autosomal Dominant Polycystic Kidney Disease (ADPKD)",
        "mondo_id": "MONDO:0008304",
        "omim_id": "173900",
        "inheritance_mode": "Autosomal dominant",
    },
    {
        "disease_id": "MONDO:0010894",
        "disease_name": "Maturity-Onset Diabetes of the Young (MODY)",
        "mondo_id": "MONDO:0010894",
        "omim_id": "600496",
        "inheritance_mode": "Autosomal dominant",
    },
    {
        "disease_id": "MONDO:0010675",
        "disease_name": "Phelan-McDermid Syndrome",
        "mondo_id": "MONDO:0010675",
        "omim_id": "606232",
        "inheritance_mode": "Autosomal dominant",
    },
]

# Master Gene Curation Records
GENE_CATALOG: List[Dict[str, Any]] = [
    # --------------------------------------------------------------------------
    # Group 1: Re-Appraised & Disputed/Refuted Cardiovascular Genes (Class A / C candidates)
    # --------------------------------------------------------------------------
    {
        "gene_symbol": "MYLK2",
        "hgnc_id": "HGNC:7591",
        "entrez_id": "85366",
        "gene_name": "Myosin light chain kinase 2",
        "chromosome": "20q11.21",
        "curations": [{
            "disease_id": "MONDO:0005045",
            "classification": "Disputed",
            "classification_date": "2019-11-15",
            "gcep": "ClinGen Hypertrophic Cardiomyopathy GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_7591_1",
            "pmids": "31784485, 27532257",
        }],
        "clinvar_variants": [
            {
                "variation_id": 14210,
                "variant_name": "NM_005597.4(MYLK2):c.355G>A (p.Glu119Lys)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 4,
                "submitter_names": "Invitae; GeneDx; Ambry Genetics; Laboratory for Molecular Medicine",
                "last_evaluated": "2018-05-12",
                "phenotype_list": "Hypertrophic cardiomyopathy; Familial hypertrophic cardiomyopathy 1",
            },
            {
                "variation_id": 180231,
                "variant_name": "NM_005597.4(MYLK2):c.821T>C (p.Val274Ala)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Blueprint Genetics",
                "last_evaluated": "2017-09-20",
                "phenotype_list": "Cardiomyopathy",
            },
            {
                "variation_id": 551092,
                "variant_name": "NM_005597.4(MYLK2):c.1147G>C (p.Ala383Pro)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Uncertain significance",
                "is_plp": 0,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 3,
                "submitter_names": "GeneDx; Illumina Clinical Services",
                "last_evaluated": "2021-02-14",
                "phenotype_list": "Hypertrophic cardiomyopathy",
            },
        ],
        "panelapp_entries": [
            {
                "panel_id": 134,
                "panel_name": "Hypertrophic cardiomyopathy - adult",
                "panel_version": "2.4",
                "confidence": "green",  # Inconsistent! Class C
                "mode_of_inheritance": "MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted",
                "phenotypes": "HCM; sudden cardiac death",
            },
            {
                "panel_id": 218,
                "panel_name": "Sudden unexplained cardiac death",
                "panel_version": "1.8",
                "confidence": "amber",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Arrhythmia; SCD",
            },
        ],
    },
    {
        "gene_symbol": "MYOZ2",
        "hgnc_id": "HGNC:15951",
        "entrez_id": "51741",
        "gene_name": "Myozenin 2",
        "chromosome": "4q26",
        "curations": [{
            "disease_id": "MONDO:0005045",
            "classification": "Disputed",
            "classification_date": "2020-04-20",
            "gcep": "ClinGen Hypertrophic Cardiomyopathy GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_15951_1",
            "pmids": "31784485, 11062473",
        }],
        "clinvar_variants": [
            {
                "variation_id": 36244,
                "variant_name": "NM_016599.5(MYOZ2):c.491A>G (p.Gln164Arg)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 3,
                "submitter_names": "LMM; Invitae; GeneDx",
                "last_evaluated": "2019-03-10",
                "phenotype_list": "Familial hypertrophic cardiomyopathy 16",
            },
            {
                "variation_id": 84129,
                "variant_name": "NM_016599.5(MYOZ2):c.734C>T (p.Pro245Leu)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Illumina Clinical Services Laboratory",
                "last_evaluated": "2018-11-05",
                "phenotype_list": "Hypertrophic cardiomyopathy",
            },
        ],
        "panelapp_entries": [
            {
                "panel_id": 134,
                "panel_name": "Hypertrophic cardiomyopathy - adult",
                "panel_version": "2.1",
                "confidence": "red",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "HCM",
            }
        ],
    },
    {
        "gene_symbol": "KCNE2",
        "hgnc_id": "HGNC:6242",
        "entrez_id": "9992",
        "gene_name": "Potassium voltage-gated channel subfamily E regulatory subunit 2",
        "chromosome": "21q22.11",
        "curations": [{
            "disease_id": "MONDO:0001625",
            "classification": "Refuted",
            "classification_date": "2020-09-12",
            "gcep": "ClinGen Arrhythmia GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_6242_1",
            "pmids": "32675988, 10471489",
        }],
        "clinvar_variants": [
            {
                "variation_id": 67201,
                "variant_name": "NM_172257.2(KCNE2):c.253T>C (p.Tyr85His)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 5,
                "submitter_names": "GeneDx; Invitae; Ambry Genetics; EGL Genetic Diagnostics",
                "last_evaluated": "2019-06-18",
                "phenotype_list": "Long QT syndrome type 6; Arrhythmia",
            },
            {
                "variation_id": 194012,
                "variant_name": "NM_172257.2(KCNE2):c.160A>G (p.Met54Val)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Laboratory for Molecular Medicine",
                "last_evaluated": "2018-02-11",
                "phenotype_list": "Long QT syndrome",
            },
        ],
        "panelapp_entries": [
            {
                "panel_id": 137,
                "panel_name": "Inherited Long QT syndrome",
                "panel_version": "3.0",
                "confidence": "green",  # Inconsistent! Class C
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Syncope; ventricular fibrillation; prolonged QTc",
            }
        ],
    },
    {
        "gene_symbol": "SCN1B",
        "hgnc_id": "HGNC:10586",
        "entrez_id": "6324",
        "gene_name": "Sodium voltage-gated channel beta subunit 1",
        "chromosome": "19q13.11",
        "curations": [{
            "disease_id": "MONDO:0015263",
            "classification": "Refuted",
            "classification_date": "2021-02-18",
            "gcep": "ClinGen Brugada Syndrome GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_10586_1",
            "pmids": "33499824, 18562304",
        }],
        "clinvar_variants": [
            {
                "variation_id": 48312,
                "variant_name": "NM_001037.5(SCN1B):c.629T>C (p.Leu210Pro)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 2,
                "submitter_names": "Invitae; GeneDx",
                "last_evaluated": "2020-01-14",
                "phenotype_list": "Brugada syndrome 5",
            },
            {
                "variation_id": 201192,
                "variant_name": "NM_001037.5(SCN1B):c.536C>T (p.Thr179Met)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Ambry Genetics",
                "last_evaluated": "2019-10-30",
                "phenotype_list": "Cardiac arrhythmia",
            },
        ],
        "panelapp_entries": [
            {
                "panel_id": 143,
                "panel_name": "Brugada syndrome",
                "panel_version": "2.2",
                "confidence": "green",  # Inconsistent! Class C
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Brugada ECG pattern; ventricular tachycardia",
            }
        ],
    },
    {
        "gene_symbol": "GPD1L",
        "hgnc_id": "HGNC:28475",
        "entrez_id": "23171",
        "gene_name": "Glycerol-3-phosphate dehydrogenase 1 like",
        "chromosome": "3p22.3",
        "curations": [{
            "disease_id": "MONDO:0015263",
            "classification": "Refuted",
            "classification_date": "2020-06-10",
            "gcep": "ClinGen Brugada Syndrome GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_28475_1",
            "pmids": "33499824, 17565368",
        }],
        "clinvar_variants": [
            {
                "variation_id": 92110,
                "variant_name": "NM_015141.4(GPD1L):c.866A>G (p.Glu289Gly)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "OMIM; Baylor Genetics",
                "last_evaluated": "2017-04-19",
                "phenotype_list": "Brugada syndrome 2",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 143,
                "panel_name": "Brugada syndrome",
                "panel_version": "2.1",
                "confidence": "red",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Brugada syndrome",
            }
        ],
    },
    {
        "gene_symbol": "TCAP",
        "hgnc_id": "HGNC:11611",
        "entrez_id": "8557",
        "gene_name": "Titin-cap",
        "chromosome": "17q12",
        "curations": [{
            "disease_id": "MONDO:0005045",
            "classification": "Disputed",
            "classification_date": "2019-12-04",
            "gcep": "ClinGen Hypertrophic Cardiomyopathy GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_11611_1",
            "pmids": "31784485, 11479158",
        }],
        "clinvar_variants": [
            {
                "variation_id": 41209,
                "variant_name": "NM_003673.4(TCAP):c.26G>A (p.Arg9Gln)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "GeneDx",
                "last_evaluated": "2018-08-15",
                "phenotype_list": "Familial hypertrophic cardiomyopathy 10",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 134,
                "panel_name": "Hypertrophic cardiomyopathy - adult",
                "panel_version": "2.3",
                "confidence": "amber",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "HCM",
            }
        ],
    },
    {
        "gene_symbol": "CACNA2D1",
        "hgnc_id": "HGNC:1399",
        "entrez_id": "781",
        "gene_name": "Calcium voltage-gated channel auxiliary subunit alpha2delta 1",
        "chromosome": "7q21.11",
        "curations": [{
            "disease_id": "MONDO:0015263",
            "classification": "Refuted",
            "classification_date": "2020-03-14",
            "gcep": "ClinGen Brugada Syndrome GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_1399_1",
            "pmids": "33499824, 20664032",
        }],
        "clinvar_variants": [
            {
                "variation_id": 133201,
                "variant_name": "NM_000722.4(CACNA2D1):c.2008C>T (p.Arg670Trp)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Laboratory for Molecular Medicine",
                "last_evaluated": "2016-12-05",
                "phenotype_list": "Brugada syndrome 3",
            }
        ],
        "panelapp_entries": [],
    },
    {
        "gene_symbol": "ANK2",
        "hgnc_id": "HGNC:493",
        "entrez_id": "287",
        "gene_name": "Ankyrin 2",
        "chromosome": "4q22.1",
        "curations": [{
            "disease_id": "MONDO:0001625",
            "classification": "Disputed",
            "classification_date": "2020-05-18",
            "gcep": "ClinGen Arrhythmia GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_493_1",
            "pmids": "32675988, 12567173",
        }],
        "clinvar_variants": [
            {
                "variation_id": 51902,
                "variant_name": "NM_001148.5(ANK2):c.4210C>T (p.Arg1404Trp)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 3,
                "submitter_names": "Invitae; GeneDx; Ambry Genetics",
                "last_evaluated": "2019-01-22",
                "phenotype_list": "Long QT syndrome 4; Ankyrin-B syndrome",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 137,
                "panel_name": "Inherited Long QT syndrome",
                "panel_version": "2.8",
                "confidence": "amber",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Sinus bradycardia; LQTS",
            }
        ],
    },
    {
        "gene_symbol": "CAV3",
        "hgnc_id": "HGNC:1529",
        "entrez_id": "859",
        "gene_name": "Caveolin 3",
        "chromosome": "3p25.3",
        "curations": [{
            "disease_id": "MONDO:0001625",
            "classification": "Disputed",
            "classification_date": "2020-08-22",
            "gcep": "ClinGen Arrhythmia GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_1529_1",
            "pmids": "32675988, 16935999",
        }],
        "clinvar_variants": [
            {
                "variation_id": 17820,
                "variant_name": "NM_033337.3(CAV3):c.233C>T (p.Thr78Met)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "GeneDx",
                "last_evaluated": "2017-03-12",
                "phenotype_list": "Long QT syndrome 9",
            }
        ],
        "panelapp_entries": [],
    },
    {
        "gene_symbol": "KCNJ8",
        "hgnc_id": "HGNC:6265",
        "entrez_id": "3764",
        "gene_name": "Potassium inwardly rectifying channel subfamily J member 8",
        "chromosome": "12p12.1",
        "curations": [{
            "disease_id": "MONDO:0015263",
            "classification": "Disputed",
            "classification_date": "2021-01-15",
            "gcep": "ClinGen Brugada Syndrome GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_6265_1",
            "pmids": "33499824, 21876532",
        }],
        "clinvar_variants": [
            {
                "variation_id": 88319,
                "variant_name": "NM_004982.4(KCNJ8):c.1009C>T (p.Arg337Trp)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Invitae",
                "last_evaluated": "2018-09-08",
                "phenotype_list": "Brugada syndrome",
            }
        ],
        "panelapp_entries": [],
    },

    # --------------------------------------------------------------------------
    # Group 2: Limited-Evidence Genes (Class B candidates)
    # --------------------------------------------------------------------------
    {
        "gene_symbol": "NLGN4X",
        "hgnc_id": "HGNC:14287",
        "entrez_id": "57599",
        "gene_name": "Neuroligin 4 X-linked",
        "chromosome": "Xp22.32-p22.31",
        "curations": [{
            "disease_id": "MONDO:0005260",
            "classification": "Limited",
            "classification_date": "2021-07-14",
            "gcep": "ClinGen Neurodevelopmental Disorders GCEP",
            "moi": "XL",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_14287_1",
            "pmids": "34212345, 12669065",
        }],
        "clinvar_variants": [
            {
                "variation_id": 98112,
                "variant_name": "NM_020742.4(NLGN4X):c.1186C>T (p.Arg396*)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "GeneDx",
                "last_evaluated": "2019-04-10",
                "phenotype_list": "Autism spectrum disorder, X-linked",
            },
            {
                "variation_id": 204910,
                "variant_name": "NM_020742.4(NLGN4X):c.725G>A (p.Arg242His)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 2,
                "submitter_names": "Invitae; Ambry Genetics",
                "last_evaluated": "2020-11-19",
                "phenotype_list": "Neurodevelopmental disorder",
            },
        ],
        "panelapp_entries": [
            {
                "panel_id": 234,
                "panel_name": "Intellectual disability and autism",
                "panel_version": "3.1",
                "confidence": "amber",
                "mode_of_inheritance": "X-LINKED: hemizygous mutation in males",
                "phenotypes": "Autism; ID",
            }
        ],
    },
    {
        "gene_symbol": "RELN",
        "hgnc_id": "HGNC:9957",
        "entrez_id": "5649",
        "gene_name": "Reelin",
        "chromosome": "7q22.1",
        "curations": [{
            "disease_id": "MONDO:0005260",
            "classification": "Limited",
            "classification_date": "2021-11-05",
            "gcep": "ClinGen Neurodevelopmental Disorders GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_9957_1",
            "pmids": "34212345, 15880352",
        }],
        "clinvar_variants": [
            {
                "variation_id": 115042,
                "variant_name": "NM_005045.4(RELN):c.5812G>A (p.Val1938Met)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Baylor Genetics",
                "last_evaluated": "2018-06-25",
                "phenotype_list": "Autism spectrum disorder",
            }
        ],
        "panelapp_entries": [],
    },
    {
        "gene_symbol": "BARD1",
        "hgnc_id": "HGNC:952",
        "entrez_id": "580",
        "gene_name": "BRCA1 associated RING domain 1",
        "chromosome": "2q35",
        "curations": [{
            "disease_id": "MONDO:0011448",
            "classification": "Limited",
            "classification_date": "2021-09-17",
            "gcep": "ClinGen Hereditary Breast Cancer GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_952_1",
            "pmids": "33947868, 26976419",
        }],
        "clinvar_variants": [
            {
                "variation_id": 128912,
                "variant_name": "NM_000465.4(BARD1):c.1659C>G (p.Cys553Trp)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 4,
                "submitter_names": "Myriad Genetic Laboratories; Ambry; GeneDx; Invitae",
                "last_evaluated": "2020-03-15",
                "phenotype_list": "Hereditary breast and ovarian cancer",
            },
            {
                "variation_id": 240188,
                "variant_name": "NM_000465.4(BARD1):c.1972C>T (p.Arg658Cys)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Color Genomics",
                "last_evaluated": "2019-12-08",
                "phenotype_list": "Breast cancer susceptibility",
            },
        ],
        "panelapp_entries": [
            {
                "panel_id": 55,
                "panel_name": "Inherited breast and ovarian cancer",
                "panel_version": "2.9",
                "confidence": "amber",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Breast cancer",
            }
        ],
    },
    {
        "gene_symbol": "RAD51D",
        "hgnc_id": "HGNC:9823",
        "entrez_id": "5892",
        "gene_name": "RAD51 paralog D",
        "chromosome": "17q12",
        "curations": [{
            "disease_id": "MONDO:0011448",
            "classification": "Limited",
            "classification_date": "2021-08-11",
            "gcep": "ClinGen Hereditary Breast Cancer GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_9823_1",
            "pmids": "33947868, 21822267",
        }],
        "clinvar_variants": [
            {
                "variation_id": 141901,
                "variant_name": "NM_002878.3(RAD51D):c.694C>T (p.Arg232Ter)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 3,
                "submitter_names": "Ambry Genetics; Invitae; GeneDx",
                "last_evaluated": "2019-07-21",
                "phenotype_list": "Familial breast cancer",
            }
        ],
        "panelapp_entries": [],
    },
    {
        "gene_symbol": "MRE11",
        "hgnc_id": "HGNC:7230",
        "entrez_id": "4361",
        "gene_name": "MRE11 homolog, double strand break repair nuclease",
        "chromosome": "11q21",
        "curations": [{
            "disease_id": "MONDO:0011448",
            "classification": "Limited",
            "classification_date": "2021-06-25",
            "gcep": "ClinGen Hereditary Breast Cancer GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_7230_1",
            "pmids": "33947868, 16905642",
        }],
        "clinvar_variants": [
            {
                "variation_id": 150221,
                "variant_name": "NM_005591.4(MRE11):c.1897C>T (p.Arg633*)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Invitae",
                "last_evaluated": "2020-05-14",
                "phenotype_list": "Hereditary cancer-predisposing syndrome",
            }
        ],
        "panelapp_entries": [],
    },
    {
        "gene_symbol": "IFT140",
        "hgnc_id": "HGNC:29077",
        "entrez_id": "9742",
        "gene_name": "Intraflagellar transport 140",
        "chromosome": "16p13.3",
        "curations": [{
            "disease_id": "MONDO:0008304",
            "classification": "Limited",
            "classification_date": "2022-01-18",
            "gcep": "ClinGen Kidney Disease GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_29077_1",
            "pmids": "35030248",
        }],
        "clinvar_variants": [
            {
                "variation_id": 192834,
                "variant_name": "NM_014714.4(IFT140):c.2140C>T (p.Arg714*)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Natera",
                "last_evaluated": "2020-08-20",
                "phenotype_list": "Polycystic kidney disease",
            }
        ],
        "panelapp_entries": [],
    },

    # --------------------------------------------------------------------------
    # Group 3: Disputed/Refuted Cancer & Sensory Genes (Class A / C candidates)
    # --------------------------------------------------------------------------
    {
        "gene_symbol": "NBN",
        "hgnc_id": "HGNC:7652",
        "entrez_id": "4683",
        "gene_name": "Nibrin",
        "chromosome": "8q21.3",
        "curations": [{
            "disease_id": "MONDO:0011448",
            "classification": "Disputed",
            "classification_date": "2021-04-12",
            "gcep": "ClinGen Hereditary Breast Cancer GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_7652_1",
            "pmids": "33947868, 9590180",
        }],
        "clinvar_variants": [
            {
                "variation_id": 3108,
                "variant_name": "NM_002485.5(NBN):c.657_661del (p.Lys219fs)",
                "variant_type": "deletion",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 6,
                "submitter_names": "Invitae; Myriad Genetics; Ambry; GeneDx; Illumina",
                "last_evaluated": "2020-02-17",
                "phenotype_list": "Hereditary breast cancer; Nijmegen breakage syndrome",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 55,
                "panel_name": "Inherited breast and ovarian cancer",
                "panel_version": "2.7",
                "confidence": "green",  # Inconsistent! Class C
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Breast cancer susceptibility",
            }
        ],
    },
    {
        "gene_symbol": "RAD50",
        "hgnc_id": "HGNC:9816",
        "entrez_id": "10111",
        "gene_name": "RAD50 double strand break repair protein",
        "chromosome": "5q31.1",
        "curations": [{
            "disease_id": "MONDO:0011448",
            "classification": "Disputed",
            "classification_date": "2021-05-19",
            "gcep": "ClinGen Hereditary Breast Cancer GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_9816_1",
            "pmids": "33947868, 19688582",
        }],
        "clinvar_variants": [
            {
                "variation_id": 134102,
                "variant_name": "NM_005732.4(RAD50):c.1378C>T (p.Arg460*)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "GeneDx",
                "last_evaluated": "2019-09-11",
                "phenotype_list": "Hereditary cancer predisposition",
            }
        ],
        "panelapp_entries": [],
    },
    {
        "gene_symbol": "XRCC2",
        "hgnc_id": "HGNC:12826",
        "entrez_id": "7516",
        "gene_name": "X-ray repair cross complementing 2",
        "chromosome": "7q36.1",
        "curations": [{
            "disease_id": "MONDO:0011448",
            "classification": "Refuted",
            "classification_date": "2020-12-08",
            "gcep": "ClinGen Hereditary Breast Cancer GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_12826_1",
            "pmids": "33947868, 22464252",
        }],
        "clinvar_variants": [
            {
                "variation_id": 182103,
                "variant_name": "NM_005431.2(XRCC2):c.96C>A (p.Tyr32*)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Ambry Genetics",
                "last_evaluated": "2017-10-15",
                "phenotype_list": "Breast cancer susceptibility",
            }
        ],
        "panelapp_entries": [],
    },
    {
        "gene_symbol": "GJB6",
        "hgnc_id": "HGNC:4289",
        "entrez_id": "10804",
        "gene_name": "Gap junction protein beta 6",
        "chromosome": "13q12.11",
        "curations": [{
            "disease_id": "MONDO:0007817",
            "classification": "Disputed",
            "classification_date": "2019-08-11",
            "gcep": "ClinGen Hearing Loss GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_4289_1",
            "pmids": "31461234, 10508513",
        }],
        "clinvar_variants": [
            {
                "variation_id": 17290,
                "variant_name": "NM_006783.5(GJB6):c.101C>T (p.Thr34Met)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 3,
                "submitter_names": "Invitae; GeneDx; PreventionGenetics",
                "last_evaluated": "2018-04-18",
                "phenotype_list": "Autosomal dominant nonsyndromic sensorineural deafness 3B",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 189,
                "panel_name": "Hearing impairment",
                "panel_version": "2.0",
                "confidence": "green",  # Inconsistent! Class C
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Non-syndromic deafness",
            }
        ],
    },
    {
        "gene_symbol": "DFNB59",
        "hgnc_id": "HGNC:29399",
        "entrez_id": "494513",
        "gene_name": "Pejvakin (PJVK)",
        "chromosome": "2q31.2",
        "curations": [{
            "disease_id": "MONDO:0007817",
            "classification": "Disputed",
            "classification_date": "2020-01-20",
            "gcep": "ClinGen Hearing Loss GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_29399_1",
            "pmids": "31461234, 16826534",
        }],
        "clinvar_variants": [
            {
                "variation_id": 221045,
                "variant_name": "NM_001042702.4(DFNB59):c.547C>T (p.Arg183Trp)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "Blueprint Genetics",
                "last_evaluated": "2018-12-04",
                "phenotype_list": "Autosomal dominant auditory neuropathy",
            }
        ],
        "panelapp_entries": [],
    },
    {
        "gene_symbol": "MECP2",
        "hgnc_id": "HGNC:6980",
        "entrez_id": "4204",
        "gene_name": "Methyl-CpG binding protein 2",
        "chromosome": "Xq28",
        "curations": [
            {
                "disease_id": "MONDO:0005260",  # Disputed for non-syndromic autism
                "classification": "Disputed",
                "classification_date": "2020-10-15",
                "gcep": "ClinGen Neurodevelopmental Disorders GCEP",
                "moi": "XL",
                "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_6980_1",
                "pmids": "33076543, 11528392",
            }
        ],
        "clinvar_variants": [
            {
                "variation_id": 14352,
                "variant_name": "NM_004992.4(MECP2):c.763C>T (p.Arg255*)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "reviewed by expert panel",
                "stars": 3,
                "num_submitters": 7,
                "submitter_names": "RettBASE; GeneDx; Invitae; Ambry Genetics",
                "last_evaluated": "2019-11-01",
                "phenotype_list": "Autism spectrum disorder, non-syndromic; Rett syndrome",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 234,
                "panel_name": "Intellectual disability and autism",
                "panel_version": "4.0",
                "confidence": "green",  # Class C for autism indication
                "mode_of_inheritance": "X-LINKED",
                "phenotypes": "Microcephaly; ASD; Rett syndrome",
            }
        ],
    },

    # --------------------------------------------------------------------------
    # Group 4: Definitive Clinical Benchmarks (Controls — NOT Downgraded)
    # --------------------------------------------------------------------------
    {
        "gene_symbol": "MYH7",
        "hgnc_id": "HGNC:7577",
        "entrez_id": "4625",
        "gene_name": "Myosin heavy chain 7",
        "chromosome": "14q11.2",
        "curations": [{
            "disease_id": "MONDO:0005045",
            "classification": "Definitive",
            "classification_date": "2017-06-01",
            "gcep": "ClinGen Hypertrophic Cardiomyopathy GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_7577_1",
            "pmids": "28574440, 2185209",
        }],
        "clinvar_variants": [
            {
                "variation_id": 14092,
                "variant_name": "NM_000257.4(MYH7):c.1208G>A (p.Arg403Gln)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "reviewed by expert panel",
                "stars": 3,
                "num_submitters": 12,
                "submitter_names": "ClinGen Cardiovascular CDWG; Invitae; GeneDx; LMM; Ambry",
                "last_evaluated": "2021-08-10",
                "phenotype_list": "Hypertrophic cardiomyopathy 1",
            },
            {
                "variation_id": 14093,
                "variant_name": "NM_000257.4(MYH7):c.2155C>T (p.Arg719Trp)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 8,
                "submitter_names": "Invitae; GeneDx; Ambry Genetics",
                "last_evaluated": "2022-03-14",
                "phenotype_list": "Familial hypertrophic cardiomyopathy",
            },
            {
                "variation_id": 14094,
                "variant_name": "NM_000257.4(MYH7):c.432C>T (p.Thr144=)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Benign",
                "is_plp": 0,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 5,
                "submitter_names": "Invitae; GeneDx",
                "last_evaluated": "2020-05-19",
                "phenotype_list": "Cardiomyopathy",
            },
        ],
        "panelapp_entries": [
            {
                "panel_id": 134,
                "panel_name": "Hypertrophic cardiomyopathy - adult",
                "panel_version": "3.1",
                "confidence": "green",  # Correct! Not inconsistent
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "HCM",
            }
        ],
    },
    {
        "gene_symbol": "TNNT2",
        "hgnc_id": "HGNC:11949",
        "entrez_id": "7139",
        "gene_name": "Troponin T2, cardiac type",
        "chromosome": "1q32.1",
        "curations": [{
            "disease_id": "MONDO:0005045",
            "classification": "Definitive",
            "classification_date": "2017-06-01",
            "gcep": "ClinGen Hypertrophic Cardiomyopathy GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_11949_1",
            "pmids": "28574440, 7522223",
        }],
        "clinvar_variants": [
            {
                "variation_id": 12411,
                "variant_name": "NM_000364.3(TNNT2):c.275G>A (p.Arg92Gln)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "reviewed by expert panel",
                "stars": 3,
                "num_submitters": 9,
                "submitter_names": "ClinGen Cardiovascular CDWG; Invitae; GeneDx",
                "last_evaluated": "2021-09-02",
                "phenotype_list": "Hypertrophic cardiomyopathy 2",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 134,
                "panel_name": "Hypertrophic cardiomyopathy - adult",
                "panel_version": "3.1",
                "confidence": "green",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "HCM",
            }
        ],
    },
    {
        "gene_symbol": "BRCA1",
        "hgnc_id": "HGNC:1100",
        "entrez_id": "672",
        "gene_name": "BRCA1 DNA repair associated",
        "chromosome": "17q21.31",
        "curations": [{
            "disease_id": "MONDO:0011448",
            "classification": "Definitive",
            "classification_date": "2018-03-01",
            "gcep": "ClinGen Hereditary Breast Cancer GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_1100_1",
            "pmids": "29543210, 7545954",
        }],
        "clinvar_variants": [
            {
                "variation_id": 17662,
                "variant_name": "NM_007294.4(BRCA1):c.68_69del (p.Glu23fs)",
                "variant_type": "deletion",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "reviewed by expert panel",
                "stars": 3,
                "num_submitters": 15,
                "submitter_names": "ENIGMA; Ambry Genetics; Myriad; Invitae",
                "last_evaluated": "2023-01-10",
                "phenotype_list": "Hereditary breast and ovarian cancer syndrome",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 55,
                "panel_name": "Inherited breast and ovarian cancer",
                "panel_version": "3.2",
                "confidence": "green",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "HBOC",
            }
        ],
    },
    {
        "gene_symbol": "BRCA2",
        "hgnc_id": "HGNC:1101",
        "entrez_id": "675",
        "gene_name": "BRCA2 DNA repair associated",
        "chromosome": "13q13.1",
        "curations": [{
            "disease_id": "MONDO:0011448",
            "classification": "Definitive",
            "classification_date": "2018-03-01",
            "gcep": "ClinGen Hereditary Breast Cancer GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_1101_1",
            "pmids": "29543210, 8524414",
        }],
        "clinvar_variants": [
            {
                "variation_id": 9325,
                "variant_name": "NM_000059.4(BRCA2):c.5946del (p.Ser1982fs)",
                "variant_type": "deletion",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "reviewed by expert panel",
                "stars": 3,
                "num_submitters": 14,
                "submitter_names": "ENIGMA; Invitae; Myriad Genetics",
                "last_evaluated": "2023-02-15",
                "phenotype_list": "Hereditary breast and ovarian cancer syndrome",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 55,
                "panel_name": "Inherited breast and ovarian cancer",
                "panel_version": "3.2",
                "confidence": "green",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "HBOC",
            }
        ],
    },
    {
        "gene_symbol": "SHANK3",
        "hgnc_id": "HGNC:14294",
        "entrez_id": "85358",
        "gene_name": "SH3 and multiple ankyrin repeat domains 3",
        "chromosome": "22q13.33",
        "curations": [{
            "disease_id": "MONDO:0010675",
            "classification": "Definitive",
            "classification_date": "2018-07-15",
            "gcep": "ClinGen Neurodevelopmental Disorders GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_14294_1",
            "pmids": "30012345, 17179978",
        }],
        "clinvar_variants": [
            {
                "variation_id": 145901,
                "variant_name": "NM_001372044.2(SHANK3):c.3680dup (p.Ala1228fs)",
                "variant_type": "duplication",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 4,
                "submitter_names": "GeneDx; Ambry Genetics; Invitae",
                "last_evaluated": "2022-05-18",
                "phenotype_list": "Phelan-McDermid syndrome",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 234,
                "panel_name": "Intellectual disability and autism",
                "panel_version": "3.8",
                "confidence": "green",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Global developmental delay; ASD",
            }
        ],
    },
    {
        "gene_symbol": "CHD8",
        "hgnc_id": "HGNC:20153",
        "entrez_id": "57680",
        "gene_name": "Chromodomain helicase DNA binding protein 8",
        "chromosome": "14q11.2",
        "curations": [{
            "disease_id": "MONDO:0005260",
            "classification": "Definitive",
            "classification_date": "2019-03-20",
            "gcep": "ClinGen Neurodevelopmental Disorders GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_20153_1",
            "pmids": "30898765, 25042858",
        }],
        "clinvar_variants": [
            {
                "variation_id": 218902,
                "variant_name": "NM_001170629.2(CHD8):c.2140C>T (p.Gln714*)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 3,
                "submitter_names": "GeneDx; Ambry Genetics",
                "last_evaluated": "2021-10-12",
                "phenotype_list": "Autism spectrum disorder with macrocephaly",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 234,
                "panel_name": "Intellectual disability and autism",
                "panel_version": "3.8",
                "confidence": "green",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Macrocephaly; ASD",
            }
        ],
    },
    {
        "gene_symbol": "KCNQ1",
        "hgnc_id": "HGNC:6289",
        "entrez_id": "3784",
        "gene_name": "Potassium voltage-gated channel subfamily Q member 1",
        "chromosome": "11p15.5-p15.4",
        "curations": [{
            "disease_id": "MONDO:0001625",
            "classification": "Definitive",
            "classification_date": "2017-04-10",
            "gcep": "ClinGen Arrhythmia GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_6289_1",
            "pmids": "28412345, 8528249",
        }],
        "clinvar_variants": [
            {
                "variation_id": 3134,
                "variant_name": "NM_000218.3(KCNQ1):c.1744C>T (p.Arg582Cys)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "reviewed by expert panel",
                "stars": 3,
                "num_submitters": 8,
                "submitter_names": "ClinGen Cardiovascular CDWG; Invitae; GeneDx",
                "last_evaluated": "2022-04-01",
                "phenotype_list": "Long QT syndrome 1",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 137,
                "panel_name": "Inherited Long QT syndrome",
                "panel_version": "3.1",
                "confidence": "green",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "LQTS1",
            }
        ],
    },
    {
        "gene_symbol": "KCNH2",
        "hgnc_id": "HGNC:6251",
        "entrez_id": "3757",
        "gene_name": "Potassium voltage-gated channel subfamily H member 2",
        "chromosome": "7q36.1",
        "curations": [{
            "disease_id": "MONDO:0001625",
            "classification": "Definitive",
            "classification_date": "2017-04-10",
            "gcep": "ClinGen Arrhythmia GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_6251_1",
            "pmids": "28412345, 7724170",
        }],
        "clinvar_variants": [
            {
                "variation_id": 3138,
                "variant_name": "NM_000238.4(KCNH2):c.1832G>A (p.Arg611His)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "reviewed by expert panel",
                "stars": 3,
                "num_submitters": 10,
                "submitter_names": "ClinGen Cardiovascular CDWG; Invitae; GeneDx",
                "last_evaluated": "2022-06-20",
                "phenotype_list": "Long QT syndrome 2",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 137,
                "panel_name": "Inherited Long QT syndrome",
                "panel_version": "3.1",
                "confidence": "green",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "LQTS2",
            }
        ],
    },
    {
        "gene_symbol": "SCN5A",
        "hgnc_id": "HGNC:10593",
        "entrez_id": "6331",
        "gene_name": "Sodium voltage-gated channel alpha subunit 5",
        "chromosome": "3p22.2",
        "curations": [{
            "disease_id": "MONDO:0015263",
            "classification": "Definitive",
            "classification_date": "2018-05-15",
            "gcep": "ClinGen Brugada Syndrome GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_10593_1",
            "pmids": "29876543, 9497256",
        }],
        "clinvar_variants": [
            {
                "variation_id": 9384,
                "variant_name": "NM_198056.3(SCN5A):c.4372C>T (p.Arg1458Ter)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "reviewed by expert panel",
                "stars": 3,
                "num_submitters": 8,
                "submitter_names": "ClinGen Cardiovascular CDWG; Invitae",
                "last_evaluated": "2023-01-20",
                "phenotype_list": "Brugada syndrome 1",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 143,
                "panel_name": "Brugada syndrome",
                "panel_version": "3.0",
                "confidence": "green",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "Brugada syndrome",
            }
        ],
    },
    {
        "gene_symbol": "HNF1A",
        "hgnc_id": "HGNC:11621",
        "entrez_id": "6927",
        "gene_name": "HNF1 homeobox A",
        "chromosome": "12q24.31",
        "curations": [{
            "disease_id": "MONDO:0010894",
            "classification": "Definitive",
            "classification_date": "2019-09-10",
            "gcep": "ClinGen Monogenic Diabetes GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_11621_1",
            "pmids": "31543210, 8825583",
        }],
        "clinvar_variants": [
            {
                "variation_id": 14450,
                "variant_name": "NM_000545.8(HNF1A):c.872dup (p.Gly292fs)",
                "variant_type": "duplication",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "reviewed by expert panel",
                "stars": 3,
                "num_submitters": 5,
                "submitter_names": "Monogenic Diabetes Variant Curation Expert Panel",
                "last_evaluated": "2021-12-05",
                "phenotype_list": "MODY3; diabetes mellitus",
            }
        ],
        "panelapp_entries": [],
    },
    {
        "gene_symbol": "PKD1",
        "hgnc_id": "HGNC:9008",
        "entrez_id": "5310",
        "gene_name": "Polycystin 1, transient receptor potential channel interacting",
        "chromosome": "16p13.3",
        "curations": [{
            "disease_id": "MONDO:0008304",
            "classification": "Definitive",
            "classification_date": "2019-04-12",
            "gcep": "ClinGen Kidney Disease GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_9008_1",
            "pmids": "31043210, 7780076",
        }],
        "clinvar_variants": [
            {
                "variation_id": 284102,
                "variant_name": "NM_001009944.3(PKD1):c.10006C>T (p.Arg3336*)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 4,
                "submitter_names": "GeneDx; Invitae; Natera",
                "last_evaluated": "2022-07-15",
                "phenotype_list": "Autosomal dominant polycystic kidney disease",
            }
        ],
        "panelapp_entries": [],
    },
    {
        "gene_symbol": "PALB2",
        "hgnc_id": "HGNC:26144",
        "entrez_id": "79728",
        "gene_name": "Partner and localizer of BRCA2",
        "chromosome": "16p12.2",
        "curations": [{
            "disease_id": "MONDO:0011448",
            "classification": "Definitive",
            "classification_date": "2019-10-10",
            "gcep": "ClinGen Hereditary Breast Cancer GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_26144_1",
            "pmids": "31654321, 17200668",
        }],
        "clinvar_variants": [
            {
                "variation_id": 126584,
                "variant_name": "NM_024675.4(PALB2):c.3113G>A (p.Trp1038*)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, multiple submitters, no conflicts",
                "stars": 2,
                "num_submitters": 7,
                "submitter_names": "Ambry; Invitae; GeneDx; Myriad Genetics",
                "last_evaluated": "2023-03-10",
                "phenotype_list": "Hereditary breast cancer",
            }
        ],
        "panelapp_entries": [
            {
                "panel_id": 55,
                "panel_name": "Inherited breast and ovarian cancer",
                "panel_version": "3.2",
                "confidence": "green",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "HBOC",
            }
        ],
    },
    # --------------------------------------------------------------------------
    # Group 4: Currently Supported Genes with Emerging Discordance (Objective 5 ML Triage)
    # --------------------------------------------------------------------------
    {
        "gene_symbol": "ACTN2",
        "hgnc_id": "HGNC:164",
        "entrez_id": "88",
        "gene_name": "Actinin alpha 2",
        "chromosome": "1q43",
        "curations": [{
            "disease_id": "MONDO:0005045",
            "classification": "Moderate",
            "classification_date": "2018-04-12",
            "gcep": "ClinGen Hypertrophic Cardiomyopathy GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_164_1",
            "pmids": "29653963, 10766170",
        }],
        "clinvar_variants": [
            {
                "variation_id": 180101,
                "variant_name": "NM_001103.3(ACTN2):c.752C>T (p.Ala251Val)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 1,
                "submitter_names": "CardioVascular Diagnostic Lab",
                "last_evaluated": "2018-02-14",
                "phenotype_list": "Hypertrophic cardiomyopathy",
            },
            {
                "variation_id": 180102,
                "variant_name": "NM_001103.3(ACTN2):c.1402G>A (p.Glu468Lys)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "no assertion criteria provided",
                "stars": 0,
                "num_submitters": 1,
                "submitter_names": "Laboratory Corp",
                "last_evaluated": "2017-09-20",
                "phenotype_list": "Cardiomyopathy",
            },
        ],
        "panelapp_entries": [
            {
                "panel_id": 134,
                "panel_name": "Hypertrophic cardiomyopathy - adult",
                "panel_version": "2.0",
                "confidence": "amber",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "HCM",
            }
        ],
    },
    {
        "gene_symbol": "FLNC",
        "hgnc_id": "HGNC:3756",
        "entrez_id": "2318",
        "gene_name": "Filamin C",
        "chromosome": "7q32.1",
        "curations": [{
            "disease_id": "MONDO:0005014",
            "classification": "Moderate",
            "classification_date": "2019-02-15",
            "gcep": "ClinGen Dilated Cardiomyopathy GCEP",
            "moi": "AD",
            "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_3756_1",
            "pmids": "30833202, 27814427",
        }],
        "clinvar_variants": [
            {
                "variation_id": 180201,
                "variant_name": "NM_001458.5(FLNC):c.6415C>T (p.Arg2139Cys)",
                "variant_type": "single nucleotide variant",
                "clinical_significance": "Likely pathogenic",
                "is_plp": 1,
                "review_status": "criteria provided, single submitter",
                "stars": 1,
                "num_submitters": 2,
                "submitter_names": "GeneDx; Blueprint Genetics",
                "last_evaluated": "2019-08-11",
                "phenotype_list": "Dilated cardiomyopathy",
            },
        ],
        "panelapp_entries": [
            {
                "panel_id": 135,
                "panel_name": "Dilated cardiomyopathy - adult",
                "panel_version": "3.0",
                "confidence": "green",
                "mode_of_inheritance": "MONOALLELIC, autosomal",
                "phenotypes": "DCM",
            }
        ],
    },
]


def seed_database(
    db_path: Optional[Union[str, Path]] = None,
    drop_existing: bool = True,
    run_audit_after_seed: bool = True,
) -> Dict[str, Any]:
    """Populate the database with the realistic seed dataset and compute audits.

    Args:
        db_path: Target SQLite database path.
        drop_existing: If True, recreates schema from scratch.
        run_audit_after_seed: If True, executes run_audit() immediately.

    Returns:
        Audit report metrics dictionary.
    """
    ensure_database_ready(db_path=db_path, drop_existing=drop_existing)

    # 1. Insert Diseases
    disease_rows = [
        (
            d["disease_id"],
            d["disease_name"],
            d.get("mondo_id"),
            d.get("omim_id"),
            d.get("inheritance_mode"),
        )
        for d in SAMPLE_DISEASES
    ]
    execute_many(
        """
        INSERT INTO diseases (disease_id, disease_name, mondo_id, omim_id, inheritance_mode)
        VALUES (?, ?, ?, ?, ?);
        """,
        disease_rows,
        db_path=db_path,
    )

    # 2. Prepare Genes, Curations, Variants, and Panels
    gene_rows = []
    curation_rows = []
    variant_rows = []
    panel_rows = []

    for item in GENE_CATALOG:
        symbol = item["gene_symbol"]
        curations = item.get("curations", [])

        # Determine downgrade status
        classifications = [c["classification"] for c in curations]
        is_downgraded = 1 if any(c in C.DOWNGRADED_CLASSIFICATIONS for c in classifications) else 0

        # Highest severity calculation
        highest_severity = "None"
        for sev in ("Refuted", "Disputed", "Limited", "Moderate", "Strong", "Definitive"):
            if sev in classifications:
                highest_severity = sev
                break

        gene_rows.append((
            symbol,
            item["hgnc_id"],
            item.get("entrez_id"),
            item.get("gene_name"),
            item.get("chromosome"),
            is_downgraded,
            highest_severity,
        ))

        # Curations
        for c in curations:
            curation_rows.append((
                symbol,
                c["disease_id"],
                c.get("moi"),
                c["classification"],
                c.get("classification_date"),
                c.get("gcep"),
                c.get("sop", "SOP v8"),
                c.get("report_url"),
                c.get("pmids"),
            ))

        # ClinVar Variants
        for v in item.get("clinvar_variants", []):
            variant_rows.append((
                v["variation_id"],
                symbol,
                curations[0]["disease_id"] if curations else None,
                v["variant_name"],
                v.get("variant_type", "single nucleotide variant"),
                v["clinical_significance"],
                v.get("is_plp", 1),
                v.get("review_status", "criteria provided, single submitter"),
                v.get("stars", 1),
                v.get("num_submitters", 1),
                v.get("submitter_names"),
                v.get("last_evaluated"),
                v.get("rcv_accessions"),
                v.get("phenotype_list"),
            ))

        # PanelApp Entries
        for p in item.get("panelapp_entries", []):
            panel_rows.append((
                symbol,
                curations[0]["disease_id"] if curations else None,
                p["panel_id"],
                p["panel_name"],
                p.get("panel_version", "1.0"),
                p.get("confidence", "amber").lower(),
                p.get("mode_of_inheritance"),
                p.get("phenotypes"),
            ))

    # Insert into database
    execute_many(
        """
        INSERT INTO genes (gene_symbol, hgnc_id, entrez_id, gene_name, chromosome, is_downgraded, highest_downgrade_severity)
        VALUES (?, ?, ?, ?, ?, ?, ?);
        """,
        gene_rows,
        db_path=db_path,
    )

    execute_many(
        """
        INSERT INTO clingen_curations (
            gene_symbol, disease_id, moi, classification, classification_date, gcep, sop, report_url, pmids
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);
        """,
        curation_rows,
        db_path=db_path,
    )

    execute_many(
        """
        INSERT INTO clinvar_variants (
            variation_id, gene_symbol, disease_id, variant_name, variant_type, clinical_significance,
            is_plp, review_status, stars, num_submitters, submitter_names, last_evaluated, rcv_accessions, phenotype_list
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
        """,
        variant_rows,
        db_path=db_path,
    )

    execute_many(
        """
        INSERT INTO panelapp_entries (
            gene_symbol, disease_id, panel_id, panel_name, panel_version, confidence, mode_of_inheritance, phenotypes
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?);
        """,
        panel_rows,
        db_path=db_path,
    )

    print(f"Successfully seeded {len(gene_rows)} genes, {len(curation_rows)} curations, "
          f"{len(variant_rows)} variants, and {len(panel_rows)} panel entries.")

    audit_result = {}
    if run_audit_after_seed:
        audit_result = run_audit(db_path=db_path)

    # Enrich database with gnomAD population frequencies and Orphanet MCAF thresholds (Objective 4)
    enrich_database_population_genetics(db_path=db_path)

    # Train and materialize downgrade risk forecasting model (Objective 5)
    run_forecasting_pipeline(db_path=db_path)

    return audit_result


def main():
    """CLI Entry point for data seeding."""
    parser = argparse.ArgumentParser(description="Seed CuraCascade database with realistic sample data.")
    parser.add_argument("--db", type=str, default=str(C.DB_PATH), help="Target database path.")
    parser.add_argument("--no-audit", action="store_true", help="Skip immediate audit generation.")
    args = parser.parse_args()

    print(f"--- Seeding CuraCascade Database: {args.db} ---")
    seed_database(
        db_path=args.db,
        drop_existing=True,
        run_audit_after_seed=not args.no_audit,
    )
    print("Seeding complete! You can now run 'streamlit run app.py' to explore the dashboard.")


if __name__ == "__main__":
    main()
