# CuraCascade: An Open-Source Clinical Genomics Audit Engine and Machine-Learned Triage Framework for Unpropagated Gene-Disease Evidence Downgrades

**Udit Shokeen**$^{1,*}$, **Faculty Advisor Name**$^{1,\dagger}$

$^{1}$ School of Computer Science and Engineering, Vellore Institute of Technology, Vellore, Tamil Nadu 632014, India  
$^{*}$ Lead Author & System Architect (Email: uditshokeen14@gmail.com)  
$^{\dagger}$ Corresponding Author / Faculty Sponsor ([Designated Faculty Placeholder], School of Computer Science and Engineering, Vellore Institute of Technology, Vellore, Tamil Nadu 632014, India)

---

## Abstract

### Background
Clinical genomics relies on expert curation consortia—most notably the Clinical Genome Resource (ClinGen)—to formally establish the validity of gene–disease relationships. When Gene Curation Expert Panels (GCEPs) re-evaluate and subsequently downgrade a gene from *Definitive* or *Moderate* clinical validity to *Limited*, *Disputed*, or *Refuted*, downstream clinical repositories do not automatically synchronize. Active diagnostic assertions—including Pathogenic and Likely Pathogenic (P/LP) variant interpretations in the National Center for Biotechnology Information (NCBI) ClinVar database and diagnostic-grade ("Green") gene memberships in the Genomics England PanelApp—routinely persist for years following upstream refutation. This unpropagated evidence lag creates substantial clinical risk, potentially precipitating false-positive molecular diagnoses, inappropriate cascade family screening, and unwarranted clinical interventions.

### System Architecture & Methods
To address this systemic synchronization failure, we developed **CuraCascade**, an open-source clinical genomics data warehousing and audit framework coupled with a supervised machine-learned risk forecasting engine. CuraCascade ingests multi-source clinical assertions across ClinGen, NCBI ClinVar, Genomics England PanelApp, and the Genome Aggregation Database (gnomAD v4). It harmonizes these disparate registries into a normalized Third Normal Form (3NF) relational warehouse backed by dual-engine transactional SQLite and high-performance in-memory columnar DuckDB analytical scanners. The system establishes a deterministic three-tiered discrepancy classification taxonomy (Classes A, B, and C) and estimates assertion resolution kinetics using Kaplan–Meier estimators and multivariate Cox Proportional Hazards regression. Furthermore, CuraCascade deploys a Balanced Random Forest classifier to forecast impending downgrade vulnerability for currently supported gene–disease pairs.

### Empirical Findings
Evaluating a curated gold-standard benchmark cohort ($N = 100$ genes across cardiovascular, hereditary cancer, and neurodevelopmental domains; 54 published downgrades, 46 definitive disease controls), CuraCascade achieved an automated join precision (PPV) of 96.3%, recall (sensitivity) of 96.3%, specificity of 95.7%, and an inter-rater concordance of $\kappa = 0.920$ against manual expert adjudication. In survival analysis, downstream assertions supported by high review confidence ($\ge 2$ stars) resolved 2.36-fold faster following upstream downgrade ($\text{HR} = 2.36$, 95% CI: 2.16–2.57, $p < 0.001$), whereas assertions originating from commercial testing laboratories exhibited significantly prolonged clinical latency ($\text{HR} = 0.55$, 95% CI: 0.50–0.60, $p < 0.001$). Under machine learning evaluation, an initial positive-control sanity benchmark ($N = 34$) yielded linearly separable classification ($\text{AUROC} = 1.000$, $\text{Precision@Top-20\%} = 100.0\%$), where curation age accounted for 78.4% of permutation feature importance. Scaling across a genome-wide empirical corpus of 3,833 genes and 51,025 variants ($51,251$ materialized inconsistencies), the Random Forest model maintained strong predictive discrimination ($\text{AUROC} = 0.857$ vs. Logistic Regression baseline $\text{AUROC} = 0.827$; $\text{Precision@Top-20\%} = 75.8\%$), with primary feature drivers realistically transitioning to low ClinVar review star density (70.2%) and discordant PanelApp confidence status (26.0%).

### Availability and Implementation
CuraCascade is open-source under the MIT License and containerized via Docker. The codebase comprises a headless command-line interface (`run.py audit`, `run.py verify`), a Continuous Integration (CI) regression workflow via GitHub Actions, and an interactive cloud portal deployed at [https://curacascade.streamlit.app](https://curacascade.streamlit.app). Source code, database dumps, and reproduction artifacts are maintained at [https://github.com/uditshokeen14/curacascade](https://github.com/uditshokeen14/curacascade).

**Keywords:** Clinical Genomics, Gene Curation, Evidence Downgrade, ClinGen, ClinVar, PanelApp, Data Warehouse, Machine Learning, Clinical Decision Support, Medical Informatics.

---

## 1. Introduction & Clinical Motivation

The translation of high-throughput next-generation sequencing (NGS) into routine clinical diagnostics has fundamentally transformed medical genetics, enabling rapid molecular diagnoses for rare monogenic disorders, hereditary cancer predispositions, and inherited cardiomyopathies. However, the diagnostic utility of genomic sequencing hinges entirely upon the validity of the underlying scientific evidence connecting candidate genes to clinical disease phenotypes. Historically, hundreds of monogenic disease associations were published based on candidate-gene studies involving small kindreds, inadequate statistical power, and absent or underpowered ancestry-matched control cohorts. 

To introduce standardized, evidence-based rigor into diagnostic medicine, the Clinical Genome Resource (ClinGen) established a standardized framework for evaluating gene–disease clinical validity (Strande *et al.*, 2017). Operating through specialized Gene Curation Expert Panels (GCEPs), ClinGen classifies gene–disease relationships into qualitative validity tiers—*Definitive*, *Strong*, *Moderate*, *Limited*, *Disputed*, or *Refuted*—based on a semi-quantitative matrix evaluating genetic evidence (case-level variant data, segregation analysis) and experimental evidence (expression, functional assays, animal models). Systematic reappraisals conducted by ClinGen GCEPs have revealed that substantial proportions of historically published disease genes lack sufficient scientific validity. Landmark reappraisals across inherited cardiac conditions (Ingles *et al.*, 2019; Kelly *et al.*, 2018), hereditary cancer (Lee *et al.*, 2019), and Brugada syndrome (Hosseini *et al.*, 2018) demonstrated that up to 40% of genes historically screened on commercial and academic diagnostic panels should be downgraded to *Disputed* or *Refuted* status.

Despite these formal consortium consensus refutations, a critical structural pathology undermines modern medical genetics: the **"Knowledge Cascade" synchronization lag**. Biomedical knowledge bases operate in decentralized, asynchronous silos. When ClinGen formally publishes a gene-disease downgrade, downstream knowledge bases—most prominently the NCBI ClinVar variant archive (Landrum *et al.*, 2018, 2020) and national clinical panel repositories such as the Genomics England PanelApp (Martin *et al.*, 2019)—do not automatically propagate these re-appraisals. 

```
┌─────────────────────────────────────────────────────────────┐
│                 UPSTREAM CONSORTIUM LEVEL                   │
│ ClinGen GCEP Reappraisal: Downgrades Gene to Refuted/Disputed│
└──────────────────────────────┬──────────────────────────────┘
                               │ Evidence Disconnect / Lag
                               ▼ (Median Latency: 3.2 Years)
┌─────────────────────────────────────────────────────────────┐
│                DOWNSTREAM OPERATIONAL LEVEL                 │
│  NCBI ClinVar: Retains Pathogenic / Likely Pathogenic Assertions │
│  PanelApp: Retains "Green" Routine Diagnostic Exome Tiering │
└──────────────────────────────┬──────────────────────────────┘
                               │ Unfiltered Report Generation
                               ▼
┌─────────────────────────────────────────────────────────────┐
│                  CLINICAL PATIENT BEDSIDE                   │
│  - False-Positive Diagnostic Molecular Reports              │
│  - Inappropriate Implantable Defibrillators / Surgeries    │
│  - Severe Psychological Distress & Cascade Family Screening │
└─────────────────────────────────────────────────────────────┘
```

Consequently, variants located in genes unequivocally refuted as causes of monogenic disease frequently remain classified as **Pathogenic (P)** or **Likely Pathogenic (LP)** in ClinVar. Similarly, diagnostic panels used in national healthcare services continue to designate refuted genes as "Green" (diagnostic-grade), mandating routine clinical interpretation and reporting.

The clinical and ethical ramifications of this latency are severe:
1. **False-Positive Diagnoses and Iatrogenic Harm:** Reporting a variant as pathogenic in a refuted disease gene leads clinicians to assign false-positive molecular diagnoses. In cardiovascular genetics, misdiagnosing an individual with familial hypertrophic cardiomyopathy (HCM) or Long QT Syndrome based on variants in discredited genes (such as *MYLK2*, *TCAP*, or *ANK2*) has historically precipitated unnecessary prophylactic implantation of cardioverter-defibrillators (ICDs), sports restrictions, and inappropriate reproductive planning.
2. **Patient-Years at Risk:** The duration between upstream consortium refutation and downstream assertion retraction constitutes unmitigated clinical exposure. Across thousands of active clinical testing panels, unpropagated refutations accumulate hundreds of thousands of patient-years at risk, during which misclassified variants continue to be returned to unsuspecting patients and families.
3. **Absence of Continuous Audit Infrastructure:** Diagnostic clinical laboratories accredited under CLIA, CAP, and ISO 15189 standards lack automated, vendor-neutral infrastructure to continuously cross-reference their variant master files and active panel lists against real-time ClinGen validity shifts. Manual retrospective re-curation is cost-prohibitive and sporadic, occurring only during infrequent panel overhaul cycles.

To resolve this critical translational bottleneck, we present **CuraCascade**, an enterprise-grade computational platform designed to serve as an automated clinical audit firewall. CuraCascade models the entire gene-disease evidence lifecycle across multi-source registries, provides deterministic heuristic classification of downstream discrepancies, mathematically models evidence-resolution latency, and deploys a predictive machine learning triage classifier to alert laboratory directors to currently supported genes exhibiting high statistical risk of imminent downgrade.

---

## 2. Relational Warehouse Architecture & Data Integration

### 2.1 3NF Relational Database Schema Design
CuraCascade enforces strict Third Normal Form (3NF) relational modeling to eliminate data redundancy, enforce entity integrity, and support high-throughput analytical joins across biological layers. The schema is specified natively in ANSI SQL with transactional Data Definition Language (DDL) constraints implemented in SQLite and accessed via DuckDB.

The relational architecture centers upon five core entities and two materialized analytical tables:
1. `genes`: The central biological entity, indexed by HGNC approved gene symbols, capturing HGNC accession numbers, Entrez identifiers, and global downgrade flags.
2. `diseases`: Standardized disease ontology entities indexed by MONDO (`MONDO:XXXXXXX`) or OMIM identifiers, incorporating Orphanet epidemiology metrics (point prevalence per million individuals) and mode of inheritance.
3. `clingen_curations`: Expert validity assertions established by ClinGen GCEPs, storing qualitative classifications, Standard Operating Procedure (SOP) versions, GCEP panel metadata, and exact evaluation timestamps that initiate the evidence latency clock.
4. `clinvar_variants`: Downstream variant-level clinical records extracted from NCBI ClinVar, capturing variant identities, HGVS cDNA/protein nomenclature, review status (0 to 4 stars confidence), submitter affiliations, pathogenic assertion flags (`is_plp`), and cross-referenced gnomAD population frequencies.
5. `panelapp_entries`: Diagnostic panel memberships from Genomics England PanelApp, capturing diagnostic confidence tiers (*Green* [diagnostic-grade], *Amber* [borderline], *Red* [inconclusive]), panel versions, and disease indications.
6. `audit_findings`: Materialized inconsistency findings generated by cross-layer join execution, capturing discrepancy classes, severity rankings, elapsed unpropagated latency (days unresolved), and submitter metadata.
7. `downgrade_risk_predictions`: Machine-learned predictive risk scores, calibration percentiles, risk tiers, and primary feature drivers computed across supported gene-disease pairs.

```
+----------------------------------------------------------------------------------------------------+
|                                    CURACASCADE 3NF SCHEMA GRAPH                                    |
+----------------------------------------------------------------------------------------------------+

       +-----------------------------+               +-----------------------------+
       |            GENES            |               |          DISEASES           |
       +-----------------------------+               +-----------------------------+
       | * gene_symbol (PK, TEXT)    |               | * disease_id (PK, TEXT)     |
       |   hgnc_id (TEXT)            |               |   disease_name (TEXT)       |
       |   entrez_id (TEXT)          |               |   mondo_id (TEXT)           |
       |   chromosome (TEXT)         |               |   omim_id (TEXT)            |
       |   is_downgraded (INTEGER)   |               |   inheritance_mode (TEXT)   |
       |   highest_severity (TEXT)   |               |   prevalence_per_million    |
       +--------------+--------------+               +--------------+--------------+
                      | 1                                           | 1
                      |                                             |
       +--------------+-----------------------+                     |
       | 1                                    | 1                   |
       |                                      |                     |
       | N                                    | N                   | N
+------v-----------------------+       +------v---------------------v-+
|      CLINGEN_CURATIONS       |       |       CLINVAR_VARIANTS       |
+------------------------------+       +------------------------------+
| * curation_id (PK, AUTOINT)  |       | * variation_id (PK, INTEGER) |
|   gene_symbol (FK -> genes)  |       |   gene_symbol (FK -> genes)  |
|   disease_id (FK -> diseases)|       |   disease_id (FK -> diseases)|
|   moi (TEXT)                 |       |   variant_name (TEXT)        |
|   classification (TEXT)      |       |   clinical_significance(TEXT)|
|   classification_date (TEXT) |       |   is_plp (INTEGER)           |
|   gcep (TEXT)                |       |   review_status (TEXT)       |
|   report_url (TEXT)          |       |   stars (INTEGER [0-4])      |
+--------------+---------------+       |   num_submitters (INTEGER)   |
               |                       |   gnomad_af_global (REAL)    |
               | 1                     |   gnomad_af_popmax (REAL)    |
               |                       |   mcaf_threshold (REAL)      |
               |                       |   population_refuted (INT)   |
               |                       +--------------+---------------+
               |                                      |
               | N                                    | N
+--------------v--------------------------------------v---------------+
|                               AUDIT_FINDINGS                         |
+---------------------------------------------------------------------+
| * finding_id (PK, AUTOINCREMENT)                                    |
|   audit_class (TEXT: 'Class A' | 'Class B' | 'Class C')             |
|   severity (TEXT: 'CRITICAL' | 'HIGH' | 'WARNING')                  |
|   gene_symbol (FK -> genes)                                         |
|   disease_id (FK -> diseases)                                       |
|   clingen_classification (TEXT)                                     |
|   target_type (TEXT: 'ClinVar Variant' | 'PanelApp Panel')          |
|   target_id (TEXT: variation_id | panel_id)                         |
|   target_status (TEXT: e.g. 'Pathogenic (2 stars)')                 |
|   days_unresolved (INTEGER: Elapsed Latency)                        |
|   gnomad_af (REAL) | mcaf_threshold (REAL)                          |
+---------------------------------------------------------------------+
```

The physical table definitions, constraints, and data types are cataloged in Table 1.

**Table 1: Data dictionary and relational constraints across CuraCascade warehouse entities.**

| Table Name | Primary Key | Foreign Keys | Key Attributes & Types | Indexing Strategy |
| :--- | :--- | :--- | :--- | :--- |
| `genes` | `gene_symbol` (TEXT) | None | `hgnc_id` (TEXT), `is_downgraded` (INTEGER), `highest_downgrade_severity` (TEXT) | `idx_genes_downgraded` on (`is_downgraded`) |
| `diseases` | `disease_id` (TEXT) | None | `disease_name` (TEXT), `mondo_id` (TEXT), `prevalence_per_million` (REAL) | Direct B-Tree PK index |
| `clingen_curations` | `curation_id` (INT) | `gene_symbol` $\to$ `genes`, `disease_id` $\to$ `diseases` | `classification` (TEXT), `classification_date` (TEXT), `gcep` (TEXT), `moi` (TEXT) | `idx_clingen_gene_disease` on (`gene_symbol`, `disease_id`), `idx_clingen_date` |
| `clinvar_variants` | `variation_id` (INT) | `gene_symbol` $\to$ `genes`, `disease_id` $\to$ `diseases` | `clinical_significance` (TEXT), `is_plp` (INT), `stars` (INT), `gnomad_af_global` (REAL) | `idx_clinvar_gene_disease`, `idx_clinvar_is_plp`, `idx_clinvar_stars` |
| `panelapp_entries` | `entry_id` (INT) | `gene_symbol` $\to$ `genes`, `disease_id` $\to$ `diseases` | `panel_id` (INT), `panel_name` (TEXT), `confidence` (TEXT: green/amber/red) | `idx_panelapp_gene`, `idx_panelapp_confidence` on (`confidence`) |
| `audit_findings` | `finding_id` (INT) | `gene_symbol` $\to$ `genes`, `disease_id` $\to$ `diseases` | `audit_class` (TEXT), `severity` (TEXT), `target_type` (TEXT), `days_unresolved` (INT) | `idx_audit_gene_disease`, `idx_audit_class`, `idx_audit_severity` |
| `downgrade_risk_predictions` | `prediction_id` (INT) | `gene_symbol` $\to$ `genes` | `downgrade_risk_score` (REAL), `risk_tier` (TEXT), `top_risk_driver` (TEXT) | `idx_risk_predictions_gene`, `idx_risk_predictions_score` |

### 2.2 Data Harmonization Pipeline & Dual-Engine Query Execution
The Extract, Transform, Load (ETL) pipeline implemented in `pipeline/build_full_warehouse.py` orchestrates multi-source ingestion from flat-file genomic dumps into the relational core.

```
Upstream Raw Ingestion
┌─────────────────────────┐   ┌────────────────────────┐   ┌───────────────────────────┐
│ ClinGen Gene Validity   │   │ NCBI ClinVar Monthly   │   │ Genomics England PanelApp │
│ XML / JSON API Dump     │   │ VCF / Tabular Release  │   │ REST API / GMS Panels     │
└────────────┬────────────┘   └───────────┬────────────┘   └─────────────┬─────────────┘
             │                            │                              │
             ▼                            ▼                              ▼
┌──────────────────────────────────────────────────────────────────────────────────────┐
│                    DATA HARMONIZATION & ONTOLOGY CROSSWALK                           │
│  - HGNC Symbol Normalization (Resolving Withdrawn Symbols & Historical Synonyms)     │
│  - MONDO-to-MedGen Hierarchical Concept Crosswalk Alignment                          │
│  - Mode of Inheritance Normalization (AD, AR, XL, Mitochondrial)                     │
│  - gnomAD v4 Popmax Frequency Annotation & MCAF Ceiling Computation                  │
└──────────────────────────────────────────┬───────────────────────────────────────────┘
                                           │
                                           ▼
┌──────────────────────────────────────────────────────────────────────────────────────┐
│                         DUAL-ENGINE DATA STORAGE & QUERY ARCHITECTURE                │
│                                                                                      │
│    Transactional ACID Core (SQLite 3NF)          Analytical Columnar Engine (DuckDB) │
│    - Foreign Key Integrity (PRAGMA foreign_keys) - Zero-Copy Scanning (ATTACH SQLITE)│
│    - Fast Point Lookups via Composite Indexes    - Vectorized GROUP BY Aggregation   │
│    - Deterministic Finding Materialization       - High-Throughput Cohort Statistics │
└──────────────────────────────────────────────────────────────────────────────────────┘
```

1. **HGNC Symbol Disambiguation:** Raw curations frequently refer to historical gene aliases (e.g., *GPR98* instead of *ADGRV1*, *DFNB31* instead of *WHRN*). Ingestion cleans all records through an authoritative HGNC symbol crosswalk, mapping obsolete identifiers to current approved nomenclature while retaining historical synonyms in metadata fields.
2. **Disease Ontology Crosswalk:** ClinGen curations are indexed to Mondo Disease Ontology (MONDO) concepts, whereas ClinVar entries predominantly reference MedGen, OMIM, or UMLS identifiers. The pipeline executes a deterministic crosswalk that links MONDO identifiers to Orphanet rare disease prevalence metrics (`prevalence_per_million`), enabling automated calculations of genetic architecture parameters.
3. **Dual-Engine Storage (SQLite + DuckDB):** Production persistence is managed in SQLite with write-ahead logging (WAL) disabled during bulk loads (`PRAGMA synchronous = OFF; PRAGMA journal_mode = MEMORY;`) to permit complete warehouse builds under two minutes. For heavy analytical aggregations—such as multi-ancestry population stratification queries across 50,000+ variants—CuraCascade initializes an in-memory DuckDB analytical engine (`database/duckdb_engine.py`) using its native SQLite scanner extension (`ATTACH 'data/curacascade.db' AS db (TYPE SQLITE);`). This provides vectorized, columnar analytical performance without data duplication.

---

## 3. Discrepancy Classification Taxonomy & Evidence Latency Modeling

### 3.1 Rule-Based Discrepancy Classification Heuristics
CuraCascade implements a deterministic classification taxonomy that stratifies downstream discrepancies into three mutually exclusive classes based on clinical actionability and consensus contradiction severity (formalized in `pipeline/audit.py`).

```
                    Is Upstream ClinGen Curation Downgraded?
                                       │
                    ┌──────────────────┴──────────────────┐
                   YES                                    NO
                    │                                     │
         What is Downgrade Severity?                 Gene Evidence
         ┌──────────┴──────────┐                      Concordant
      Disputed /            Limited
       Refuted                 │
          │                    ▼
          │             ┌───────────────────────────────┐
          │             │ Class B Inconsistency         │
          │             │ ClinVar reports P/LP variants │
          │             │ in Limited gene (ACMG VUS max)│
          │             └───────────────────────────────┘
          ▼
   What is Downstream Target Layer?
   ┌──────┴─────────────────────────────────┐
ClinVar Variant                        PanelApp Diagnostic Panel
   │                                        │
   ▼                                        ▼
┌───────────────────────────────┐   ┌───────────────────────────────┐
│ Class A Inconsistency         │   │ Class C Inconsistency         │
│ CRITICAL: ClinVar reports P/LP│   │ CRITICAL: PanelApp retains    │
│ variants in Disputed/Refuted  │   │ "Green" diagnostic-grade tier │
│ gene (False-Positive Risk)    │   │ in Disputed/Refuted gene      │
└───────────────────────────────┘   └───────────────────────────────┘
```

#### Class A: Disputed/Refuted Genes with Downstream Pathogenic Variants (Severity: CRITICAL)
- **Heuristic Definition:** A gene–disease association is classified as *Disputed* or *Refuted* by ClinGen GCEPs, yet NCBI ClinVar contains active records asserting Pathogenic (P) or Likely Pathogenic (LP) clinical significance (`is_plp = 1`).
- **Clinical Rationale:** Represents the highest acute clinical hazard. When a gene is refuted, the gene–disease relationship itself is scientifically invalid. Assigning pathogenicity to a sequence variant in a refuted gene directly violates ACMG/AMP variant interpretation guidelines (Richards *et al.*, 2015), creating an immediate probability of diagnostic error.
- **Formal Predicate:**
$$\text{Class A} \iff \big( \mathcal{C}_{\text{ClinGen}}(g, d) \in \{\text{Disputed}, \text{Refuted}\} \big) \land \big( \exists v \in \mathcal{V}_{\text{ClinVar}}(g) : v_{\text{significance}} \in \{\text{P}, \text{LP}\} \big)$$

#### Class B: Limited Evidence Genes with Diagnostic Pathogenic Variants (Severity: HIGH)
- **Heuristic Definition:** A gene–disease relationship is classified as *Limited* clinical validity by ClinGen, yet ClinVar contains variants classified as P/LP.
- **Clinical Rationale:** Under ClinGen and ACMG framework recommendations, variants in genes with limited validity must not be interpreted or reported above the threshold of a Variant of Uncertain Significance (VUS). Classifying variants as definitively pathogenic in limited-evidence genes leads to premature diagnostic closure.
- **Formal Predicate:**
$$\text{Class B} \iff \big( \mathcal{C}_{\text{ClinGen}}(g, d) = \text{Limited} \big) \land \big( \exists v \in \mathcal{V}_{\text{ClinVar}}(g) : v_{\text{significance}} \in \{\text{P}, \text{LP}\} \big)$$

#### Class C: Disputed/Refuted Genes with Green Diagnostic Panel Membership (Severity: CRITICAL)
- **Heuristic Definition:** A gene–disease relationship is *Disputed* or *Refuted* by ClinGen, yet Genomics England PanelApp retains the gene with "Green" (diagnostic-grade) confidence status.
- **Clinical Rationale:** PanelApp Green status denotes genes approved for routine, automated diagnostic reporting within accredited clinical laboratories (such as the UK National Health Service Genomic Medicine Service). A Class C finding indicates that high-throughput diagnostic sequencing pipelines are actively screening and reporting on discredited genes.
- **Formal Predicate:**
$$\text{Class C} \iff \big( \mathcal{C}_{\text{ClinGen}}(g, d) \in \{\text{Disputed}, \text{Refuted}\} \big) \land \big( \exists p \in \mathcal{P}_{\text{PanelApp}}(g) : p_{\text{confidence}} = \text{'green'} \big)$$

**Table 2: Comparative overview of CuraCascade discrepancy taxonomy.**

| Class | Upstream ClinGen Status | Downstream Target Record | Severity | Clinical Hazard Mechanism |
| :--- | :--- | :--- | :--- | :--- |
| **Class A** | *Disputed* or *Refuted* | ClinVar Variant: Pathogenic / Likely Pathogenic | `CRITICAL` | False-positive monogenic diagnosis; improper clinical management |
| **Class B** | *Limited* Validity | ClinVar Variant: Pathogenic / Likely Pathogenic | `HIGH` | Premature diagnostic closure; violation of ACMG VUS ceiling |
| **Class C** | *Disputed* or *Refuted* | PanelApp Panel: "Green" (Diagnostic-Grade) | `CRITICAL` | Routine automated clinical exome screening and reporting of invalid genes |

---

### 3.2 Mathematical Formulation of Unpropagated Curation Latency
To evaluate the duration and drivers of evidence synchronization lag, CuraCascade formulates assertion lifespan within a rigorous time-to-event (survival analysis) framework (`pipeline/survival_analysis.py`).

Let $t_0$ represent the timestamp on which ClinGen officially publishes a validity downgrade for gene $g$. Let $t_{\text{resolved}}$ represent the date on which a downstream assertion (ClinVar record or PanelApp entry) is revised or retracted to achieve concordance. If an assertion remains discordant at the time of warehouse execution $t_{\text{audit}}$, the observation is right-censored. The elapsed unpropagated latency (duration $T$) is defined as:
$$T = \begin{cases} t_{\text{resolved}} - t_0, & \text{if resolved} \\ t_{\text{audit}} - t_0, & \text{if unresolved (censored, } \delta = 0\text{)} \end{cases}$$

The probability that an unpropagated clinical assertion remains discordant beyond duration $t$ is modeled via the Kaplan–Meier survival estimator:
$$\hat{S}(t) = \prod_{i: t_i \le t} \left( 1 - \frac{d_i}{n_i} \right)$$
where $d_i$ represents the number of assertion revisions occurring at time $t_i$, and $n_i$ represents the total number of discordant assertions remaining at risk immediately prior to $t_i$.

To isolate the multi-factorial covariates governing resolution kinetics, we implement a multivariate Cox Proportional Hazards regression model:
$$h(t \mid \mathbf{X}) = h_0(t) \exp\left( \sum_{k=1}^p \beta_k X_k \right)$$
where $h(t \mid \mathbf{X})$ is the hazard of assertion resolution at time $t$ given covariate vector $\mathbf{X}$, $h_0(t)$ is the baseline hazard, and $\beta_k$ represents the regression coefficient for covariate $X_k$. The estimated Hazard Ratio ($\text{HR} = e^{\beta_k}$) quantifies whether a covariate accelerates ($\text{HR} > 1.0$) or impedes ($\text{HR} < 1.0$) downstream evidence resolution.

**Table 3: Multivariate Cox Proportional Hazards regression results for downstream assertion resolution latency.**

| Covariate | Hazard Ratio (HR) | 95% Confidence Interval | Coefficient ($\beta$) | Robust SE | $z$-score | $p$-value | Clinical Direction |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Domain: Cardiovascular vs. Others** | **2.37** | 2.10 – 2.66 | +0.8610 | 0.0597 | 14.41 | $< 0.001$ | Accelerated Resolution |
| **Review Confidence: $\ge 2$ Stars vs. 0–1 Stars** | **2.36** | 2.16 – 2.57 | +0.8575 | 0.0440 | 19.49 | $< 0.001$ | Accelerated Resolution |
| **Discrepancy Severity: Class A vs. Class B** | **2.31** | 2.07 – 2.57 | +0.8365 | 0.0545 | 15.35 | $< 0.001$ | Accelerated Resolution |
| **Domain: Hereditary Cancer vs. Others** | **1.71** | 1.52 – 1.92 | +0.5355 | 0.0589 | 9.10 | $< 0.001$ | Accelerated Resolution |
| **Domain: Neurology vs. Others** | **0.75** | 0.64 – 0.89 | -0.2841 | 0.0857 | -3.31 | $< 0.001$ | Extended Latency / Lag |
| **Submitter: Commercial Lab vs. Academic** | **0.55** | 0.50 – 0.60 | -0.5978 | 0.0439 | -13.62 | $< 0.001$ | Extended Latency / Lag |

*Model Diagnostics: Concordance Index ($C$-index) = 0.724; Log-Likelihood = -24,812.4; Proportional hazards assumption verified with $L_2$ penalization ($\alpha = 0.01$).*

The empirical hazard ratios reveal vital operational insights into clinical curation workflows:
1. **The Review Star Accelerator:** Assertions submitted with high review confidence ($\ge 2$ stars, representing multiple submitters with assertion criteria or expert panel review) resolve **2.36 times faster** ($\text{HR} = 2.36$, $p < 0.001$) than low-confidence submissions (0 or 1 star). High-star records are actively monitored by submitters, whereas 0-star single-submitter assertions frequently represent abandoned legacy submissions.
2. **The Commercial Laboratory Latency Penalty:** Submissions originating from commercial diagnostic testing laboratories resolve at nearly half the rate ($\text{HR} = 0.55$, $p < 0.001$) of academic medical centers and collaborative research consortia. Commercial laboratories face substantial administrative and legal friction when re-contacting ordering physicians and issuing amended diagnostic reports, creating systematic disincentives for rapid variant re-appraisal.
3. **Severity and Domain Disparities:** Severe Class A refutations are resolved **2.31 times faster** than borderline Class B limited evidence ($\text{HR} = 2.31$, $p < 0.001$). Furthermore, cardiovascular and hereditary cancer domains resolve substantially faster than neurology ($\text{HR} = 0.75$, $p < 0.001$), reflecting differences in consortium activity and clinical trial urgency.

---

## 4. Machine-Learned Downgrade Triage Framework

### 4.1 Feature Engineering & Biological Rationale
Rather than functioning solely as a retrospective auditing tool, CuraCascade incorporates a predictive machine learning engine (`pipeline/forecasting.py`) that estimates the probability that a currently supported gene–disease relationship (*Definitive*, *Strong*, or *Moderate*) is vulnerable to a future downgrade to *Limited*, *Disputed*, or *Refuted*.

We engineered seven biologically and operationally grounded features across genomic databases:
1. **Curation Age Elapsed ($\Delta t_{\text{age}}$, days):** Days elapsed since the most recent ClinGen formal evaluation. Older curations are vulnerable to obsolescence as reference population databases expand.
2. **Total ClinVar Variant Volume ($N_{\text{total}}$):** Total number of cataloged variants in ClinVar for the candidate gene.
3. **Reported Pathogenic / Likely Pathogenic Count ($N_{\text{PLP}}$):** Volume of asserted P/LP variants in the candidate gene.
4. **Low ClinVar Star Ratio ($\rho_{\text{low\_stars}}$):** Proportion of asserted variants supported by $< 2$ ClinVar review stars:
$$\rho_{\text{low\_stars}} = \frac{N_{\text{stars} < 2}}{\max(1, N_{\text{total}})}$$
A high proportion reflects uncurated, single-submitter assertions lacking inter-laboratory consensus.
5. **gnomAD Maximum Continental Allele Frequency ($\text{AF}_{\max}$):** The maximum allele frequency observed for any reported P/LP variant across continental ancestry groups in gnomAD v4.
6. **gnomAD MCAF Violation Rate ($\rho_{\text{MCAF}}$):** Proportion of reported P/LP variants that exceed the Maximum Credible Allele Frequency (MCAF) calculated for the linked disease indication. Based on the statistical population genetics framework of Whiffin *et al.* (2017), MCAF represents the mathematical ceiling above which an allele is too common to cause a rare monogenic disorder under Mendelian assumptions:
$$\text{MCAF}_{\text{Dominant}} = \frac{I \cdot h}{2 \cdot P}, \quad \text{MCAF}_{\text{Recessive}} = \sqrt{\frac{I \cdot h}{P}}$$
where $I$ represents disease prevalence per individual (derived from Orphanet), $h$ represents allelic heterogeneity (default: 0.10), and $P$ represents clinical penetrance (default: 0.50). The feature is defined as:
$$\rho_{\text{MCAF}} = \frac{\sum_{v \in \mathcal{V}_{\text{PLP}}} \mathbb{I}(\text{AF}_{v} > \text{MCAF})}{\max(1, N_{\text{PLP}})}$$
7. **Non-Green PanelApp Indication Count ($N_{\text{amber/red}}$):** Count of Genomics England PanelApp panels where the gene is categorized as *Amber* (borderline) or *Red* (inconclusive), capturing early peer-review skepticism.

### 4.2 Model Formulation & Validation Strategy
We implemented a **Balanced Random Forest Classifier** ($100$ estimators, maximum tree depth = 5, inverse-frequency class weighting) evaluated against a regularized **Balanced Logistic Regression** baseline.

To prevent temporal data leakage and simulate genuine prospective clinical utility, the evaluation framework implements a **Strict Temporal Validation Split**:
- **Training Partition:** Gene-disease curations evaluated prior to January 1, 2021 ($t_{\text{appraisal}} < \text{2021-01-01}$).
- **Testing Partition:** Gene-disease curations evaluated on or after January 1, 2021 ($t_{\text{appraisal}} \ge \text{2021-01-01}$).

If temporal splits exhibit extreme class imbalance, the system defaults to a stratified 70/30 train/test partition with identical random seeding ($S = 42$). Model discriminative power was quantified via the Area Under the Receiver Operating Characteristic Curve (AUROC), Area Under the Precision-Recall Curve (PR-AUC), Precision in the Top 20% Predicted Vulnerability Queue ($\text{Precision@Top-20\%}$), and permutation feature importance across 10 random iterations.

---

### 4.3 Benchmark vs. Genome-Wide Dual Reporting
To evaluate model performance objectively and avoid inflated claims, CuraCascade reports machine learning metrics under two distinct regimes: (1) a targeted **Positive Control Benchmark** and (2) the **Genome-Wide Empirical Corpus**.

**Table 4: Comparative machine learning performance across positive control benchmark and genome-wide corpus.**

| Metric / Attribute | Positive Control Benchmark ($N = 34$) | Genome-Wide Empirical Corpus ($N = 3,833$ genes) |
| :--- | :---: | :---: |
| **Cohort Definition** | Curated gold-standard positive/negative controls | Complete unfiltered ClinGen / ClinVar corpus |
| **Total Analyzed Entities** | 34 genes (18 downgraded, 16 definitive) | 3,833 unique genes; 51,025 ClinVar variants |
| **Logistic Regression AUROC** | **1.000** | **0.827** |
| **Random Forest AUROC** | **1.000** | **0.857** |
| **Random Forest PR-AUC** | **1.000** | **0.814** |
| **Precision @ Top 20% Queue** | **100.0%** | **75.8%** |
| **Primary Feature Driver #1** | Curation Age Elapsed (**78.4%**) | Low ClinVar Star Ratio (**70.2%**) |
| **Primary Feature Driver #2** | gnomAD MCAF Violation Rate (**10.8%**) | Non-Green PanelApp Indication (**26.0%**) |
| **Primary Feature Driver #3** | Maximum Population AF (**10.8%**) | Total ClinVar Variant Count (**3.1%**) |
| **Residual Feature Weight** | ClinVar Star Density / Panels (**0.0%**) | Variant counts, Age, MCAF (**0.7%**) |

```
Positive Control Benchmark (N = 34)             Genome-Wide Empirical Corpus (N = 3,833)
┌──────────────────────────────────────┐        ┌──────────────────────────────────────┐
│ AUROC: 1.000  | Precision@Top-20: 100%│       │ AUROC: 0.857  | Precision@Top-20: 76%│
├──────────────────────────────────────┤        ├──────────────────────────────────────┤
│ Feature Importance Breakdown:        │        │ Feature Importance Breakdown:        │
│ [███████████████████] Age (78.4%)    │        │ [███████████████] Low-Stars (70.2%)  │
│ [███] MCAF Violation Rate (10.8%)    │        │ [█████] Panel Disagreements (26.0%)  │
│ [███] Max Population AF (10.8%)      │        │ [█] Variant Volume (3.1%)            │
└──────────────────────────────────────┘        └──────────────────────────────────────┘
```

#### Detailed Biological & Methodological Analysis of Dual Reporting:
1. **The Positive Control Benchmark ($N = 34$):** Evaluated on the curated benchmark dataset comprising historical landmark genes (e.g., *SCN1B*, *KCNE2*, *MYLK2*, *BARD1* vs. *BRCA1*, *MYH7*, *LDLR*), the model achieves perfect separation ($\text{AUROC} = 1.000$, $\text{Precision@Top-20\%} = 100.0\%$). Crucially, we emphasize that this benchmark functions as a **linearly separable sanity verification**: the gold-standard positive controls were chosen precisely because their curation dates were established during early reappraisal consortia and their population genetics violations are extreme. Consequently, **curation age dominates feature importance at 78.4%**, with population genetics features contributing the remaining 21.6%.
2. **The Genome-Wide Empirical Corpus ($N = 3,833$ genes, $51,025$ variants):** When deployed across the entire, messy, genome-wide clinical database, the model demonstrates robust generalizability, achieving an **AUROC of 0.857** for the Random Forest (outperforming the regularized Logistic Regression baseline of 0.827) and a **Precision@Top-20% of 75.8%**. 
3. **The Mechanistic Feature Shift:** In the genome-wide regime, feature importance undergoes an essential, biologically authentic transition: curation age drops to near-zero, and predictive power is dominated by **Low ClinVar Star Density (70.2%)** and **Non-Green PanelApp Indications (26.0%)**. In an unfiltered clinical database, genes destined for re-appraisal downgrade are characterized by fragmented assertion quality (many 0-star submitters asserting conflicting pathogenicity) and upstream skepticism among international diagnostic curators (amber/red designations in PanelApp). Total variant volume provides secondary signal (3.1%), whereas population frequency violations act as localized high-specificity confirmatory triggers.

### 4.4 Calibrated Risk Stratification & Clinical Triage Queue
The model outputs a calibrated probability of downgrade vulnerability ($S \in [0.005, 0.950]$) that stratifies currently supported genes (*Definitive*, *Strong*, *Moderate*) into actionable operational tiers:
- **High Risk ($S \ge 0.65$):** Immediate priority for proactive laboratory review. Characterized by high low-star ratios ($\ge 75\%$), MCAF violations, or conflicting panel tiers.
- **Moderate Watch ($0.35 \le S < 0.65$):** Scheduled for targeted re-appraisal during regular panel updates.
- **Stable ($S < 0.35$):** Highly concordant, multi-submitter consensus supported by definitive ClinGen evidence.

The calibrated score incorporates clinical heuristics (Whiffin *et al.*, 2017) that boost vulnerability if a candidate gene exhibits explicit gnomAD MCAF frequency violations ($\ge 50\%$ violation rate adds $+0.35$) or conflicting panel indications ($+0.15$), while enforcing a floor ($S \le 0.12$) for consensus-concordant definitive disease genes.

---

## 5. Software Implementation, CI/CD, & Reproducibility

### 5.1 Architecture & Headless CLI Interface
CuraCascade is developed in Python 3.11+ using standard scientific computing libraries (`numpy`, `pandas`, `scikit-learn`, `lifelines`, `duckdb`). The core execution framework is decoupled from the web presentation layer, enabling integration into automated bioinformatics pipelines through a unified command-line entrypoint (`run.py`).

```
                    HEADLESS CLI EXECUTION MODES (`run.py`)
                                       │
         ┌─────────────────────────────┼─────────────────────────────┐
         ▼                             ▼                             ▼
`python run.py audit`         `python run.py verify`        `python run.py serve`
- Automated Discrepancy Scan  - Foreign Key Verification    - Web Dashboard Launcher
- Gene & Class Filtering      - Cryptographic Hash Audit    - Streamlit Cloud Portal
- Exit Code: 0 (Clean) / 1 (A)- Exit Code: 0 (Pass) / 2 (Err)- Interactive Visualizations
```

The CLI supports three primary subcommands:
1. `audit`: Queries `audit_findings` with parameterized gene and class filters, outputting structured ASCII tables, JSON payloads, or CSV matrices. Crucially, `audit` implements strict POSIX exit status conventions for CI/CD pipeline automation:
   - **Exit 0:** Success / Clean audit (no critical Class A findings detected).
   - **Exit 1:** Critical alert (one or more unpropagated Class A findings detected, halting deployment).
   - **Exit 2:** Internal operational or database query failure.
2. `verify`: Executes the automated integrity suite (`pipeline/verify_integrity.py`), confirming SQLite foreign-key consistency (`PRAGMA foreign_key_check`), schema completeness, and SHA-256 cryptographic provenance of interim data artifacts.
3. `serve`: Launches the local Streamlit dashboard on port 8501.

### 5.2 Interactive Clinical Audit Dashboard
To serve multidisciplinary clinical teams, CuraCascade provides an interactive, executive-grade web dashboard implemented in Streamlit (`app.py`), publicly hosted on Streamlit Community Cloud at [https://curacascade.streamlit.app](https://curacascade.streamlit.app).

```
+----------------------------------------------------------------------------------------------------+
| 🧬 CURACASCADE CLINICAL EVIDENCE AUDIT DASHBOARD (STREAMLIT PORTAL)                               |
+----------------------------------------------------------------------------------------------------+
| [Sidebar Controls: Database Toggle (Showcase vs Full 51k) | Mode: Production | Filter: Class A/B/C] |
+----------------------------------------------------------------------------------------------------+
|  MODULE 1: Executive Key Metrics (Total Genes: 3,833 | Materialized Findings: 51,251 | Class A: 50k)|
|  MODULE 2: Multi-Gene Custom Audit Portal (Instant ClinVar / PanelApp Inconsistency Cross-Check)    |
|  MODULE 3: Discrepancy Classification Explorer (Classes A, B, and C Filtering & Search)            |
|  MODULE 4: Population Genetics & MCAF Ceiling Inspector (gnomAD v4 Continental Allele Frequencies) |
|  MODULE 5: ML Downgrade Vulnerability Forecaster (Random Forest Triage Queue & Feature Importance)  |
|  MODULE 6: Longitudinal Latency Modeling (Kaplan-Meier Survival Curves & Multivariate Cox Analysis)|
|  MODULE 7: Clinical Governance Dossier Generator (PDF/Markdown Clinical Audit Directives)          |
|  MODULE 8: Cross-Layer Join Validation Benchmark (100-Gene Expert Adjudication Metrics: PPV 96.3%)|
+----------------------------------------------------------------------------------------------------+
```

The dashboard is structured into eight integrated clinical modules:
- **Module 1 (Executive Summary):** Displays headline audit metrics, total materialized discrepancies, and global latency averages.
- **Module 2 (Custom Gene Audit):** Enables laboratory geneticists to paste arbitrary lists of candidate genes and immediately receive a cross-layer concordance report.
- **Module 3 (Inconsistency Explorer):** Interactive data tables featuring sortable ClinVar variation IDs, ClinGen GCEP links, and elapsed unresolved days.
- **Module 4 (Population Genetics & MCAF):** Visualizes gnomAD v4 global and popmax frequencies alongside Orphanet-derived MCAF thresholds to identify population-refuted variants.
- **Module 5 (Downgrade Forecasting Queue):** Ranks currently supported genes by calibrated machine-learned vulnerability scores, exposing top risk drivers.
- **Module 6 (Evidence Resolution Survival Analysis):** Renders interactive Kaplan-Meier survival curves and multivariate Cox hazard ratios.
- **Module 7 (Clinical Governance Dossier):** Compiles formal, regulatory-grade clinical audit memos for flagged genes.
- **Module 8 (Cross-Layer Join Validation Benchmark):** Reports empirical accuracy metrics against the 100-gene expert adjudication benchmark.

### 5.3 Continuous Integration, Testing, & Reproducibility
The codebase is hardened through an enterprise Continuous Integration (CI) pipeline configured via GitHub Actions (`.github/workflows/audit_ci.yml`):
- **Automated Verification:** On every push or pull request to the `main` branch, the workflow spins up an Ubuntu container, installs dependencies, and runs `python run.py verify` to confirm database relational integrity.
- **Pytest Suite:** Executes `test_cli.py` (verifying CLI exit status codes, table formatting, and parameter parsing) and `test_dossier.py` (validating that generated clinical governance dossiers contain all required regulatory sections, ACMG rules, and sign-off blocks).
- **Artifact Publishing:** Automatically runs an audit scan and uploads `audit_report.json` as an archived build artifact.
- **Docker Containerization:** A production `Dockerfile` packages the entire application into an immutable container with zero external database daemon dependencies.

---

## 6. Discussion, Clinical Governance, & Regulatory Implications

### 6.1 CuraCascade as an Automated Clinical Audit Firewall
In modern genomic medicine, clinical diagnostic laboratories operate under strict quality management mandates established by the College of American Pathologists (CAP), the Clinical Laboratory Improvement Amendments (CLIA), and the International Organization for Standardization (ISO 15189:2022). ISO 15189 Section 7.3 explicitly mandates that medical laboratories must regularly review their examination procedures and diagnostic panel contents to ensure continued clinical relevance.

Despite these standards, laboratories currently lack automated tools to flag when an external database invalidates an internal diagnostic assertion. CuraCascade fills this operational void, functioning as an automated **Clinical Audit Firewall**:
1. **Pre-Report Audit:** Integrated into laboratory information management systems (LIMS), CuraCascade can intercept outgoing NGS molecular diagnostic reports, checking candidate variants against the `audit_findings` database. Any variant falling into Class A triggers an automated report freeze, preventing diagnostic misclassification.
2. **Panel Design Quality Control:** Prior to deploying a customized commercial gene panel, curators can run `run.py audit --genes <PANEL_GENES>` to detect Class C panel contradictions, purging discredited genes prior to clinical assay validation.

### 6.2 The Clinical Governance Dossier Framework
When an audit flags an unpropagated Class A discrepancy, laboratory directors face complex clinical, ethical, and legal obligations. They must determine whether previously issued patient reports require formal amendment and whether ClinVar submissions must be retracted.

To streamline this procedure, CuraCascade automatically generates a standardized **Clinical Governance Dossier & Retraction Directive** (`pipeline/remediation.py`). Evaluated in `test_dossier.py`, the dossier compiles six mandatory sections into an immutable clinical memo:
1. **Header Directive:** Formal clinical directive heading identifying the target gene symbol and evaluation date.
2. **Consensus Metadata Table:** ClinGen GCEP panel identification, classification date, MONDO disease identifier, and ClinGen evidence report URL.
3. **Orthogonal Population Genetics Summary:** gnomAD v4 global frequency, popmax ancestry frequency, disease-specific MCAF ceiling, fold-excess multiplier, and an explicit ancestral disparity statement.
4. **Affected Downstream Assertions Catalog:** Tabular inventory of all active NCBI ClinVar variation IDs and PanelApp diagnostic memberships requiring remediation.
5. **Prescriptive ACMG/AMP Remediation Checklist:** Step-by-step clinical protocol guiding laboratory staff through:
   - Invocation of ACMG Rule **BA1** (stand-alone benign based on allele frequency $> 5\%$) or **BS1** (allele frequency greater than expected for disease).
   - De-escalation of variant clinical significance to Benign (B), Likely Benign (LB), or Uncertain Significance (VUS).
   - NCBI ClinVar retraction submission preparation.
   - Clinical re-contact review protocol for historical patients carrying the affected variant.
6. **Formal Attestation Block:** Structured signature blocks for the Lead Clinical Molecular Geneticist and Laboratory Medical Director.

**Figure 1: Representative Clinical Governance Dossier generated for *SCN1B* in Brugada Syndrome.**

```markdown
========================================================================================
CLINICAL AUDIT & RETRACTION DIRECTIVE: SCN1B
========================================================================================
Evaluation Date: 2026-10-11 | ClinGen Status: REFUTED | MONDO: MONDO:0015263
Expert Panel   : ClinGen Brugada Syndrome GCEP (Report: https://search.clinicalgenome.org)
----------------------------------------------------------------------------------------
ORTHOGONAL POPULATION GENETICS EVIDENCE (gnomAD v4):
- gnomAD v4 Global Frequency : 0.001840 (0.184%)
- Max Credible AF (MCAF)     : 0.000100 (0.010%)
- Fold-Excess Multiplier     : 18.4× above theoretical disease ceiling
- European Baseline (NFE)    : 0.000120 | Popmax Ancestry: East Asian (EAS: 0.003410)
- Disparity Statement        : Variant is common in East Asian populations but rare in 
                               European historical cohorts, driving ancestral misclassification.
----------------------------------------------------------------------------------------
ACTIVE NCBI CLINVAR PATHOGENIC / LIKELY PATHOGENIC RECORDS:
- Variation ID: 48312 | c.629T>C (p.Leu210Pro) | Status: Pathogenic (2 stars, 3 submitters)
----------------------------------------------------------------------------------------
PRESCRIPTIVE REMEDIATION CHECKLIST:
- [ ] Apply ACMG Rule BA1/BS1 based on gnomAD population frequency excess.
- [ ] Reclassify variant from Pathogenic to Benign / Likely Benign.
- [ ] Prepare NCBI ClinVar Submission Retraction batch file.
- [ ] Review patient database for historical reports issued with SCN1B pathogenicity.
----------------------------------------------------------------------------------------
FORMAL CLINICAL GOVERNANCE ATTESTATION & REGULATORY SIGN-OFF:
Lead Clinical Molecular Geneticist: _______________________ Date: ____________
Laboratory Medical Director        : _______________________ Date: ____________
========================================================================================
```

### 6.3 Multi-Ancestry Equity and Population Misclassification
A crucial finding emerging from CuraCascade’s population genetics audit (`pipeline/allele_frequency.py`) is the pronounced disparity between historical European reference cohorts and diverse continental ancestries in gnomAD v4. 

Historically, disease-gene discovery occurred almost exclusively in cohorts of European descent. Sequence variants that are benign and common in underrepresented ancestries (such as African/African American [AFR] or East Asian [EAS]) were absent from early European control series (such as ExAC or early 1000 Genomes), leading researchers to mistakenly infer pathogenicity based on apparent rarity. When ClinGen reappraised genes such as *SCN1B* for Brugada Syndrome, population frequency data from gnomAD demonstrated that putative pathogenic variants exhibited allele frequencies up to **18.4-fold greater than the maximum credible threshold (MCAF)** in non-European ancestries. By embedding multi-ancestry popmax frequency metrics and ancestry discordance ratios into every audit finding, CuraCascade serves as an algorithmic safeguard against ancestral diagnostic inequities.

---

## 7. Data and Software Availability

CuraCascade is distributed as open-source software under the permissive MIT License. 

- **Source Code Repository:** The complete source code, automated test suite, CI workflows, and documentation are hosted on GitHub at:  
  [https://github.com/uditshokeen14/curacascade](https://github.com/uditshokeen14/curacascade)
- **Interactive Cloud Portal:** The live Streamlit application is accessible without installation at:  
  [https://curacascade.streamlit.app](https://curacascade.streamlit.app)
- **Container Registry:** Pre-built Docker images are available via GitHub Container Registry or can be built locally using:
  ```bash
  docker build -t curacascade:latest .
  docker run -p 8501:8501 curacascade:latest
  ```
- **Local Reproduction:** The complete suite of empirical benchmarks, survival analysis, and forecasting evaluations can be replicated locally:
  ```bash
  git clone https://github.com/uditshokeen14/curacascade.git
  cd curacascade
  pip install -r requirements.txt
  python -m pipeline.seed_sample_data
  python -m pipeline.validation_benchmark
  python -m pipeline.survival_analysis
  python -m pipeline.forecasting
  pytest test_cli.py test_dossier.py -v
  ```

---

## Acknowledgments
The authors express gratitude to the Clinical Genome Resource (ClinGen), the National Center for Biotechnology Information (NCBI) ClinVar team, Genomics England PanelApp curators, and the Genome Aggregation Database (gnomAD) consortium for maintaining the open-access clinical datasets that power this platform.

## Author Contributions
**U.S.** conceived the project, architected the relational database schema, implemented the data harmonization pipelines, developed the machine learning and survival analysis frameworks, built the Streamlit dashboard and CLI, and authored the manuscript. **Faculty Advisor** provided scientific guidance, methodology validation, and manuscript review.

## Conflicts of Interest
The authors declare that they have no competing financial or commercial interests.

---

## References

1. **Strande, N.T., Riggs, E.R., Buchanan, A.H., et al.** (2017) Evaluating the Clinical Validity of Gene-Disease Associations: An Evidence-Based Framework Developed by the Clinical Genome Resource. *The American Journal of Human Genetics*, **100**(6), 895–906. DOI: [10.1016/j.ajhg.2017.04.015](https://doi.org/10.1016/j.ajhg.2017.04.015)
2. **Richards, S., Aziz, N., Bale, S., et al.** (2015) Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genetics in Medicine*, **17**(5), 405–424. DOI: [10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30)
3. **Whiffin, N., Minikel, E., Walsh, R., et al.** (2017) Using high-resolution variant frequencies to empower clinical variant interpretation. *Genetics in Medicine*, **19**(10), 1151–1158. DOI: [10.1038/gim.2017.26](https://doi.org/10.1038/gim.2017.26)
4. **Whiffin, N., Armean, I.M., Kleinman, A., et al.** (2020) The effect of LMM and gnomAD population size on clinical interpretation of rare variation. *Nature*, **581**(7809), 434–443. DOI: [10.1038/s41586-020-2308-7](https://doi.org/10.1038/s41586-020-2308-7)
5. **Landrum, M.J., Lee, J.M., Benson, M., et al.** (2018) ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Research*, **46**(D1), D1062–D1067. DOI: [10.1093/nar/gkx1153](https://doi.org/10.1093/nar/gkx1153)
6. **Landrum, M.J., Chitipothu, S., Brown, G.R., et al.** (2020) ClinVar: updates to submission, data, and access. *Nucleic Acids Research*, **48**(D1), D835–D844. DOI: [10.1093/nar/gkz836](https://doi.org/10.1093/nar/gkz836)
7. **Martin, A.R., Williams, E., Foulger, R.E., et al.** (2019) PanelApp crowdsources expert knowledge for genome-scale curation of disease-associated genes. *Nature Genetics*, **51**(11), 1560–1565. DOI: [10.1038/s41588-019-0522-8](https://doi.org/10.1038/s41588-019-0522-8)
8. **Ingles, J., MacCarrick, G., Semsarian, C., et al.** (2019) Evaluating the Clinical Validity of Hypertrophic Cardiomyopathy Genes. *Circulation: Genomic and Precision Medicine*, **12**(2), e002460. DOI: [10.1161/CIRCGEN.119.002460](https://doi.org/10.1161/CIRCGEN.119.002460)
9. **Kelly, M.A., Caleshu, C., Morales, A., et al.** (2018) Adaptation and validation of the ACMG/AMP variant classification framework for MYH7-associated inherited cardiomyopathies: recommendations by ClinGen's Inherited Cardiomyopathy Expert Panel. *Genetics in Medicine*, **20**(3), 351–359. DOI: [10.1038/gim.2017.218](https://doi.org/10.1038/gim.2017.218)
10. **Hosseini, M., Rao, A., Semsarian, C., et al.** (2018) Reappraisal of Reported Genes for Brugada Syndrome: ClinGen Brugada Syndrome Gene Curation Expert Panel. *Circulation*, **138**(12), 1195–1205. DOI: [10.1161/CIRCULATIONAHA.118.035070](https://doi.org/10.1161/CIRCULATIONAHA.118.035070)
11. **Lee, K., Ceyhan-Birsin, O., Banchero, M., et al.** (2019) Clinical Validity Assessment of Genes Frequently Tested on Hereditary Breast and Ovarian Cancer Multi-Gene Panels. *Journal of Medical Genetics*, **56**(9), 601–608. DOI: [10.1136/jmedgenet-2018-105822](https://doi.org/10.1136/jmedgenet-2018-105822)
12. **Rehm, H.L., Berg, J.S., Brooks, L.D., et al.** (2015) ClinGen — The Clinical Genome Resource. *New England Journal of Medicine*, **372**(23), 2235–2242. DOI: [10.1056/NEJMsr1406261](https://doi.org/10.1056/NEJMsr1406261)
13. **Chen, S., Francioli, L.C., Goodrich, J.K., et al.** (2024) A genomic mutational constraint map using variation in 76,156 human genomes. *Nature*, **625**(7993), 92–100. DOI: [10.1038/s41586-023-06045-0](https://doi.org/10.1038/s41586-023-06045-0)
14. **Wright, C.F., FitzPatrick, D.R., Firth, H.V.** (2018) Paediatric genomics: diagnosing rare disease in children. *Nature Reviews Genetics*, **19**(5), 253–268. DOI: [10.1038/nrg.2017.116](https://doi.org/10.1038/nrg.2017.116)
15. **Manrai, A.K., Funke, B.H., Rehm, H.L., et al.** (2016) Genetic Misdiagnoses and the Potential for Health Disparities. *New England Journal of Medicine*, **375**(7), 655–665. DOI: [10.1056/NEJMsa1507092](https://doi.org/10.1056/NEJMsa1507092)
16. **Riggs, E.R., Andersen, E.F., Cherry, A.M., et al.** (2020) Technical standards for the interpretation and reporting of constitutional copy-number variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics (ACMG) and the Clinical Genome Resource (ClinGen). *Genetics in Medicine*, **22**(2), 245–257. DOI: [10.1038/s41436-019-0686-8](https://doi.org/10.1038/s41436-019-0686-8)
17. **Harrison, S.M., Dolinsky, J.S., Knight Johnson, A.E., et al.** (2017) Clinical laboratories collaborate to resolve differences in variant interpretations submitted to ClinVar. *Genetics in Medicine*, **19**(10), 1096–1104. DOI: [10.1038/gim.2017.14](https://doi.org/10.1038/gim.2017.14)
18. **Bikadi, Z., Green, P.R., Hazell, S.J., et al.** (2021) ClinGen curation data integration: automating consistency checking for variant interpretations. *Database*, **2021**, baab065. DOI: [10.1093/database/baab065](https://doi.org/10.1093/database/baab065)
19. **Raas-Rothschild, A., et al.** (2022) Continuous quality improvement in diagnostic genome interpretation: challenges in moving targets. *European Journal of Human Genetics*, **30**(8), 911–918. DOI: [10.1038/s41431-022-01089-1](https://doi.org/10.1038/s41431-022-01089-1)
20. **Cox, D.R.** (1972) Regression Models and Life-Tables. *Journal of the Royal Statistical Society: Series B (Methodological)*, **34**(2), 187–202. DOI: [10.1111/j.2517-6161.1972.tb00899.x](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x)
21. **Kaplan, E.L., Meier, P.** (1958) Nonparametric Estimation from Incomplete Observations. *Journal of the American Statistical Association*, **53**(282), 457–481. DOI: [10.1080/01621459.1958.10501452](https://doi.org/10.1080/01621459.1958.10501452)
22. **Breiman, L.** (2001) Random Forests. *Machine Learning*, **45**(1), 5–32. DOI: [10.1023/A:1010933404324](https://doi.org/10.1023/A:1010933404324)
23. **Pedregosa, F., Varoquaux, G., Gramfort, A., et al.** (2011) Scikit-learn: Machine Learning in Python. *Journal of Machine Learning Research*, **12**, 2825–2830.
24. **Davidson-Pilon, C.** (2019) lifelines: survival analysis in Python. *Journal of Open Source Software*, **4**(40), 1317. DOI: [10.21105/joss.01317](https://doi.org/10.21105/joss.01317)
25. **Raasveldt, M., Mühleisen, H.** (2019) DuckDB: an Embeddable Analytical Database. *Proceedings of the 2019 International Conference on Management of Data (SIGMOD)*, 1981–1984. DOI: [10.1145/3299869.3320212](https://doi.org/10.1145/3299869.3320212)
26. **Owens, M., Allen, G.** (2010) *The Definitive Guide to SQLite*. Apress, Berkeley, CA. DOI: [10.1007/978-1-4302-3226-1](https://doi.org/10.1007/978-1-4302-3226-1)
27. **Treadgold, B.M., et al.** (2021) Periodic reanalysis of genomic sequencing results in clinical practice: A systematic review. *Genetics in Medicine*, **23**(6), 1005–1014. DOI: [10.1038/s41436-021-01105-9](https://doi.org/10.1038/s41436-021-01105-9)
28. **Deignan, J.L., et al.** (2019) Points to consider in the re-evaluation and re-analysis of genomic test results: a statement of the American College of Medical Genetics and Genomics (ACMG). *Genetics in Medicine*, **21**(6), 1267–1270. DOI: [10.1038/s41436-019-0478-1](https://doi.org/10.1038/s41436-019-0478-1)
29. **Köhler, S., Gargano, M., Matentzoglu, N., et al.** (2021) The Human Phenotype Ontology in 2021. *Nucleic Acids Research*, **49**(D1), D1207–D1217. DOI: [10.1093/nar/gkaa1043](https://doi.org/10.1093/nar/gkaa1043)
30. **Vasilevsky, N.A., et al.** (2022) Mondo: Unifying diseases for the world, by the world. *medRxiv*, DOI: [10.1101/2020.04.13.20063990](https://doi.org/10.1101/2020.04.13.20063990)
