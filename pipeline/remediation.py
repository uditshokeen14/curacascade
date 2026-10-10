"""CuraCascade — Prescriptive Clinical Remediation & Decision Support Engine (Objective 7).

Transforms CuraCascade from an observational auditor into an active clinical governance
system by providing deterministic, actionable, step-by-step remediation protocols
for clinical diagnostic laboratories, ClinVar variant submitters, panel curators,
and medical geneticists.

Adheres strictly to established professional clinical genomics standards:
- Richards S, et al. (2015) ACMG/AMP Variant Interpretation Standards
- Strande NT, et al. (2017) ClinGen Gene-Disease Clinical Validity Framework
- Whiffin N, et al. (2017) Clinical Genomics Allele Frequency Filtering (MCAF)
- ClinGen Gene Curation Working Group SOP v8 (2022)
"""

from __future__ import annotations

import datetime
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Union

import pandas as pd

# Ensure repository root is on sys.path
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))



# ==============================================================================
# Standard Professional Guidelines Citations Registry
# ==============================================================================

STANDARD_CITATIONS: List[Dict[str, str]] = [
    {
        "framework": "ACMG/AMP Variant Interpretation Standards",
        "citation": "Richards S, et al. (2015). Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med. 17(5):405-424.",
        "doi": "10.1038/gim.2015.30",
    },
    {
        "framework": "ClinGen Gene-Disease Validity Framework",
        "citation": "Strande NT, et al. (2017). Evaluating the Clinical Validity of Gene-Disease Associations: An Evidence-Based Framework by ClinGen. Am J Hum Genet. 100(6):895-906.",
        "doi": "10.1016/j.ajhg.2017.04.015",
    },
    {
        "framework": "Population Frequency MCAF Framework",
        "citation": "Whiffin N, et al. (2017). Using high-resolution variant frequencies to empower clinical genome interpretation. Genet Med. 19(10):1151-1158.",
        "doi": "10.1038/gim.2017.26",
    },
    {
        "framework": "ClinGen Gene Curation SOP v8",
        "citation": "ClinGen Gene Curation Working Group (2022). Standard Operating Procedure for Clinical Validity of Gene-Disease Associations, Version 8. Clinical Genome Resource.",
        "url": "https://www.clinicalgenome.org/curation-activities/gene-disease-validity/",
    },
]


# ==============================================================================
# 1. Deterministic Remediation Decision Engine
# ==============================================================================

def get_prescriptive_resolution(
    audit_class: str,
    clingen_status: str,
    target_type: str = "ClinVar Variant",
    population_refuted: bool = False,
    stars: int = 1,
) -> str:
    """Compute a concise, actionable resolution badge for tabular audit views.

    Args:
        audit_class: Inconsistency class ('Class A', 'Class B', 'Class C', or empty).
        clingen_status: Authoritative ClinGen classification.
        target_type: 'ClinVar Variant' or 'PanelApp Panel'.
        population_refuted: True if gnomAD population AF exceeds MCAF.
        stars: ClinVar review stars rating.

    Returns:
        Concise badge string for table display.
    """
    cls_str = str(audit_class or "").strip()
    status_str = str(clingen_status or "").strip().title()
    t_type = str(target_type or "").strip()

    if population_refuted:
        return "🔴 Apply ACMG BA1/BS1 Benign (AF > MCAF)"

    if "Panel" in t_type or cls_str == "Class C":
        if status_str in ("Refuted", "Disputed") or cls_str == "Class C":
            return "🟣 Demote Panel to Red (Research Only)"
        elif status_str == "Limited" or cls_str == "Class B":
            return "🟡 Restrict to Amber Tier (Borderline)"
        return "🟢 Retain Green Diagnostic Tier"

    # ClinVar Variant assertions
    if status_str in ("Refuted", "Disputed") or cls_str == "Class A":
        return "🔴 Reclassify to VUS (Retract P/LP)"
    elif status_str == "Limited" or cls_str == "Class B":
        return "🟠 Cap at VUS (ACMG Guidance)"
    elif status_str in ("Definitive", "Strong", "Moderate"):
        return "🟢 Maintain Concordant Clinical Assertion"

    return "⚪ Review Assertion Consistency"


def generate_remediation_plan(
    gene_symbol: str,
    finding_class: str,
    clingen_status: str,
    variant_row: Optional[Union[Dict[str, Any], pd.Series]] = None,
    panel_row: Optional[Union[Dict[str, Any], pd.Series]] = None,
    disease_info: Optional[Dict[str, Any]] = None,
    gcep_info: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Generate a structured, 4-dimensional clinical remediation and governance protocol.

    Args:
        gene_symbol: Target HGNC gene symbol (e.g., 'ANK2', 'MYLK2', 'BARD1').
        finding_class: Audit class ('Class A', 'Class B', 'Class C', or composite).
        clingen_status: Authoritative ClinGen classification ('Disputed', 'Refuted', 'Limited', etc.).
        variant_row: Optional ClinVar variant dictionary or series.
        panel_row: Optional PanelApp panel dictionary or series.
        disease_info: Optional dictionary containing disease name and ID.
        gcep_info: Optional dictionary containing GCEP title, appraisal date, and report URL.

    Returns:
        Structured dictionary containing protocols across all four action dimensions.
    """
    sym = str(gene_symbol).strip().upper()
    status = str(clingen_status).strip().title()
    f_class = str(finding_class or "").strip()

    # Extract variant metadata
    v_data: Dict[str, Any] = {}
    if variant_row is not None:
        v_data = variant_row.to_dict() if isinstance(variant_row, pd.Series) else dict(variant_row)

    # Extract panel metadata
    p_data: Dict[str, Any] = {}
    if panel_row is not None:
        p_data = panel_row.to_dict() if isinstance(panel_row, pd.Series) else dict(panel_row)

    # Context values
    v_name = v_data.get("variant_name") or f"{sym} clinical variant"
    v_id = str(v_data.get("variation_id") or v_data.get("target_id") or "")
    v_stars = int(v_data.get("stars", 1) or 1)
    v_sig = str(v_data.get("clinical_significance") or "Pathogenic")

    p_name = p_data.get("panel_name") or f"{sym} diagnostic gene panel"
    p_id = str(p_data.get("panel_id") or "")
    p_conf = str(p_data.get("confidence") or "green").lower()

    d_name = (
        (disease_info or {}).get("disease_name")
        or v_data.get("disease_name")
        or p_data.get("disease_name")
        or "Associated Inherited Condition"
    )
    d_id = (
        (disease_info or {}).get("disease_id")
        or v_data.get("disease_id")
        or "MONDO:GENERAL"
    )

    gcep_name = (
        (gcep_info or {}).get("gcep")
        or v_data.get("clingen_gcep")
        or "ClinGen Gene Curation Expert Panel"
    )
    c_date = (
        (gcep_info or {}).get("classification_date")
        or v_data.get("clingen_date")
        or "Authoritative Appraisal"
    )
    report_url = (
        (gcep_info or {}).get("report_url")
        or v_data.get("clingen_report_url")
        or "https://search.clinicalgenome.org/kb/gene-validity"
    )

    # Check if population allele frequency exceeds MCAF threshold
    gnomad_af_g = float(v_data.get("gnomad_af_global", 0.0) or 0.0)
    gnomad_af_p = float(v_data.get("gnomad_af_popmax", 0.0) or 0.0)
    mcaf = float(v_data.get("mcaf_threshold", 0.0001) or 0.0001)
    pop_refuted = (
        v_data.get("population_refuted_flag", 0) == 1
        or (gnomad_af_g > mcaf or gnomad_af_p > mcaf)
    ) and ("PATHOGENIC" in v_sig.upper() or v_data.get("is_plp", 0) == 1)

    max_af = max(gnomad_af_g, gnomad_af_p)
    mcaf_ratio = max_af / max(1e-8, mcaf) if mcaf > 0 else 0.0

    # --------------------------------------------------------------------------
    # 1. ClinVar Submitter Protocol
    # --------------------------------------------------------------------------
    submitter_protocol: Dict[str, Any] = {}

    if pop_refuted:
        rule_code = "BA1 (Stand-Alone Benign)" if max_af >= 0.05 else "BS1 (Strong Benign)"
        submitter_protocol = {
            "action_title": "Stand-Alone Benign Reclassification (ACMG BA1 / BS1)",
            "action_badge": f"Rule: ACMG {rule_code.split()[0]} Benign",
            "badge_color": "#ef4444",
            "target_classification": "Benign / Likely Benign",
            "acmg_rules": [
                f"ACMG {rule_code}: Population AF ({max_af:.2e}) significantly exceeds disease MCAF ({mcaf:.2e}) by {mcaf_ratio:.1f}x.",
                "ClinGen Allele Frequency Threshold Standard (Whiffin et al. 2017).",
                "Retract Pathogenic / Likely Pathogenic assertion immediately.",
            ],
            "checklist": [
                f"Access NCBI ClinVar Submission Portal for Variation ID {v_id if v_id else sym}.",
                f"Change clinical significance from '{v_sig}' to 'Benign' or 'Likely Benign'.",
                f"Cite gnomAD v4 popmax allele frequency ({max_af:.2e}) exceeding MCAF threshold ({mcaf:.2e}).",
                f"Attach ClinGen GCEP determination report URL: {report_url}.",
                "Notify contributing clinicians and medical directors of variant downgrading.",
            ],
            "rationale": (
                f"Variant occurs in gnomAD v4 at a population frequency of {max_af:.2e}, exceeding the "
                f"Maximum Credible Allele Frequency (MCAF = {mcaf:.2e}) by {mcaf_ratio:.1f}-fold for {d_name}. "
                f"Under ACMG/AMP rules and Whiffin et al. (2017), this frequency constitutes stand-alone or strong "
                "evidence against disease causality."
            ),
            "evidence_anchor": f"{gcep_name} ({c_date})",
            "report_url": report_url,
        }

    elif status in ("Disputed", "Refuted") or "Class A" in f_class:
        submitter_protocol = {
            "action_title": "Formal Assertion Retraction & Reclassification to VUS",
            "action_badge": "Action Required: Reclassify to VUS",
            "badge_color": "#ef4444",
            "target_classification": "Variant of Uncertain Significance (VUS)",
            "acmg_rules": [
                "ClinGen SOP v8: Disputed/Refuted gene validity precludes Pathogenic assignment.",
                "ACMG/AMP 2015 Framework: Lack of robust human gene-disease validity overrides variant-level computational scores.",
                "Reclassify to VUS or retract assertion from public ClinVar repository.",
            ],
            "checklist": [
                f"Initiate revision of ClinVar submission SCV/RCV records for {v_name} (Variation ID {v_id}).",
                f"Downgrade assertion from '{v_sig}' to 'Variant of Uncertain Significance (VUS)'.",
                f"Reference authoritative ClinGen {status} appraisal by {gcep_name} ({c_date}).",
                f"Embed permanent ClinGen Evidence Dossier hyperlink: {report_url}.",
                "Update internal laboratory variant knowledge base to block automated P/LP calling.",
            ],
            "rationale": (
                f"ClinGen Expert Panel ({gcep_name}) formally determined the {sym}–{d_name} association "
                f"to be '{status}' on {c_date}. Continued reporting of pathogenic status in ClinVar causes false "
                "positive diagnostic calls and misinforms patient management."
            ),
            "evidence_anchor": f"{gcep_name} ({c_date})",
            "report_url": report_url,
        }

    elif status == "Limited" or "Class B" in f_class:
        submitter_protocol = {
            "action_title": "ACMG Guidance Cap: Restrict Variant Classification to VUS",
            "action_badge": "ACMG Guidance Cap: Maximum VUS",
            "badge_color": "#f59e0b",
            "target_classification": "Variant of Uncertain Significance (VUS)",
            "acmg_rules": [
                "ACMG/ClinGen 2017 Rule: Limited-evidence genes MUST NOT carry Pathogenic or Likely Pathogenic assertions.",
                "Cap classification at VUS until replicated human genetic and functional proof emerges.",
            ],
            "checklist": [
                f"Review all clinical reports in laboratory database asserting {v_name} as Pathogenic.",
                "Cap variant classification at 'Variant of Uncertain Significance (VUS)' in ClinVar.",
                f"Document ClinGen 'Limited' evidence appraisal status ({c_date}) in assertion notes.",
                f"Provide link to ClinGen Dossier: {report_url}.",
                "Flag gene for re-evaluation upon publication of future ClinGen GCEP updates.",
            ],
            "rationale": (
                f"Evidence supporting {sym} in {d_name} is graded as 'Limited'. Current ACMG and ClinGen "
                "consensus explicitly contraindicates Pathogenic or Likely Pathogenic assertions for genes "
                "lacking established definitive validity."
            ),
            "evidence_anchor": f"{gcep_name} ({c_date})",
            "report_url": report_url,
        }

    else:
        # Definitive / Strong / Moderate concordant
        submitter_protocol = {
            "action_title": "Concordant Clinical Assertion Maintenance",
            "action_badge": "Status: Concordant Curation",
            "badge_color": "#10b981",
            "target_classification": v_sig,
            "acmg_rules": [
                "Adheres to ACMG/AMP 2015 criteria for supported gene-disease relationship.",
                "ClinGen validity supported as Definitive/Strong.",
            ],
            "checklist": [
                "Maintain existing ClinVar submission in active standing.",
                "Conduct routine biennial review of newly submitted submitter assertions.",
                "Monitor gnomAD population allele frequencies for emerging drift.",
            ],
            "rationale": f"Gene-disease relationship is clinically validated as {status}. Assertions align with standards.",
            "evidence_anchor": f"{gcep_name} ({c_date})",
            "report_url": report_url,
        }

    # --------------------------------------------------------------------------
    # 2. Diagnostic Panel Maintenance (PanelApp & Diagnostic Labs)
    # --------------------------------------------------------------------------
    panel_governance: Dict[str, Any] = {}

    if status in ("Disputed", "Refuted") or "Class C" in f_class or "Class A" in f_class:
        panel_governance = {
            "action_title": "Immediate Panel Demotion to Red Tier & BED File Exclusion",
            "panel_badge": "Panel Advisory: Demote to Red Tier (Research Only)",
            "badge_color": "#ef4444",
            "target_tier": "Red (Research Only / Do Not Report)",
            "filtering_advisory": (
                f"Immediately remove genomic coordinates for {sym} from routine diagnostic exome / targeted "
                f"panel filtering BED files for {d_name}."
            ),
            "checklist": [
                f"Genomics England PanelApp: Demote {sym} on Panel '{p_name}' (ID {p_id}) from Green to Red.",
                f"Bioinformatics Pipelines: Exclude {sym} coding exons from clinical filtering BED targets.",
                f"Diagnostic Sign-Off: Prohibit reporting variants in {sym} as primary diagnostic findings for {d_name}.",
                "Route incidental discoveries to research protocols with secondary consent requirements.",
            ],
            "rationale": (
                f"Clinical gene validity for {sym} is '{status}'. Retention on diagnostic-grade 'Green' lists "
                "generates inappropriate diagnostic reports and clinical misattribution."
            ),
        }

    elif status == "Limited" or "Class B" in f_class:
        panel_governance = {
            "action_title": "Panel Restriction to Amber Tier (Borderline Confidence)",
            "panel_badge": "Panel Advisory: Restrict to Amber Tier",
            "badge_color": "#f59e0b",
            "target_tier": "Amber (Borderline / Insufficient Evidence)",
            "filtering_advisory": (
                f"Restrict {sym} to secondary or candidate reporting tier in {d_name} panels. Do not issue primary diagnostic findings."
            ),
            "checklist": [
                f"Adjust Genomics England PanelApp rating for {sym} on Panel '{p_name}' to Amber.",
                "Ensure diagnostic pipeline requires senior molecular geneticist sign-off before reporting.",
                "Annotate diagnostic reports with explicit disclosure of limited human validity evidence.",
            ],
            "rationale": (
                f"{sym} evidence is 'Limited'. While not refuted, validity is insufficient to warrant unconstrained "
                "first-line diagnostic reporting."
            ),
        }

    else:
        panel_governance = {
            "action_title": "Diagnostic-Grade Panel Inclusion Verified",
            "panel_badge": "Panel Advisory: Maintain Green Tier",
            "badge_color": "#10b981",
            "target_tier": "Green (Diagnostic-Grade)",
            "filtering_advisory": f"Retain {sym} in routine diagnostic exome and panel BED targets for {d_name}.",
            "checklist": [
                f"Confirm active Green status on diagnostic panel '{p_name}'.",
                "Ensure 100% exon coverage standards (≥30x) are met in NGS runs.",
            ],
            "rationale": f"{sym} is verified as {status} for {d_name}. Green diagnostic tier is fully warranted.",
        }

    # --------------------------------------------------------------------------
    # 3. Clinical Governance & Patient Safety Advisory
    # --------------------------------------------------------------------------
    patient_safety: Dict[str, Any] = {}

    if status in ("Disputed", "Refuted") or "Class A" in f_class or pop_refuted:
        patient_safety = {
            "alert_level": "CRITICAL",
            "alert_badge": "🚨 Critical Clinical Governance Alert",
            "badge_color": "#ef4444",
            "advisory_title": "Historical Index Case Audit & Immediate Cascade Testing Halt",
            "guidance_summary": (
                f"Medical geneticists and molecular pathology directors must conduct an immediate retrospective "
                f"audit of all historical patient cases where a diagnosis of {d_name} was established based on {sym}."
            ),
            "action_items": [
                f"Audit Historical Cases: Identify all index patients previously issued positive reports for {sym}.",
                "Halt Cascade Testing: Immediately discontinue predictive genetic testing in asymptomatic relatives.",
                f"Re-Evaluate Clinical Interventions: Re-assess invasive or irreversible management decisions "
                f"(e.g., prophylactic surgeries, pacemaker/ICD implants, intensive surveillance) predicated solely on {sym}.",
                f"Issue Amended Molecular Reports: Provide corrected clinical reports citing ClinGen {gcep_name} ({c_date}).",
                "Genetic Counseling Re-Contact: Offer updated genetic counseling to affected families to alleviate unwarranted anxiety.",
            ],
        }

    elif status == "Limited" or "Class B" in f_class:
        patient_safety = {
            "alert_level": "HIGH",
            "alert_badge": "⚠️ High Clinical Governance Advisory",
            "badge_color": "#f59e0b",
            "advisory_title": "Patient Re-Evaluation & Restrictive Cascade Protocol",
            "guidance_summary": (
                f"For patients carrying candidate variants in {sym}, inform treating clinicians of the Limited "
                "evidence status. Do not use variants as sole criteria for irreversible clinical interventions."
            ),
            "action_items": [
                "Review clinical documentation of patients with reported variants in this gene.",
                "Advise against predictive testing of asymptomatic family members based on VUS assertions.",
                "Recommend broader clinical exome re-analysis to identify bona fide pathogenic causes.",
                "Document limited validity in patient medical record and electronic health systems (EHR).",
            ],
        }

    else:
        patient_safety = {
            "alert_level": "STANDARD",
            "alert_badge": "✅ Validated Standard of Care",
            "badge_color": "#10b981",
            "advisory_title": "Routine Clinical Care & Cascade Testing Supported",
            "guidance_summary": (
                f"Pathogenic assertions in {sym} represent established causal etiology for {d_name}. "
                "Standard clinical actionability and predictive cascade family testing apply."
            ),
            "action_items": [
                "Continue standard-of-care clinical management protocols.",
                "Offer predictive genetic counseling and testing to at-risk biological relatives.",
                "Adhere to established disease-specific surveillance and treatment guidelines.",
            ],
        }

    # Concise overall prescriptive badge
    prescriptive_badge = get_prescriptive_resolution(
        audit_class=f_class,
        clingen_status=status,
        target_type="PanelApp Panel" if "Panel" in f_class else "ClinVar Variant",
        population_refuted=pop_refuted,
        stars=v_stars,
    )

    return {
        "gene_symbol": sym,
        "finding_class": f_class,
        "clingen_status": status,
        "disease_name": d_name,
        "disease_id": d_id,
        "prescriptive_badge": prescriptive_badge,
        "submitter_protocol": submitter_protocol,
        "panel_governance": panel_governance,
        "patient_safety": patient_safety,
        "guidelines_citations": STANDARD_CITATIONS,
        "generated_at": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC"),
    }


# ==============================================================================
# 2. Clinical Remediation Memo Generator (Markdown)
# ==============================================================================

def build_remediation_memo_markdown(
    gene_symbol: str,
    gene_data: Optional[Dict[str, Any]] = None,
    remediation_plan: Optional[Dict[str, Any]] = None,
) -> str:
    """Format remediation directives into a standardized, professional clinical memorandum.

    Args:
        gene_symbol: HGNC gene symbol.
        gene_data: Gene metadata dictionary (name, HGNC ID, disease, curations, variants, etc.).
        remediation_plan: Output dictionary from generate_remediation_plan().

    Returns:
        Structured Markdown text suitable for formal laboratory export and audit filing.
    """
    sym = str(gene_symbol).strip().upper()
    g_data = gene_data or {}
    r_plan = remediation_plan or generate_remediation_plan(sym, "Class A", "Disputed")

    today_str = datetime.date.today().strftime("%B %d, %Y")
    ref_id = f"CURACASCADE-REM-{sym}-{datetime.date.today().strftime('%Y%m%d')}"

    gene_name = g_data.get("gene_name") or f"{sym} Human Gene"
    hgnc_id = g_data.get("hgnc_id") or "HGNC:RECORD"
    disease_name = r_plan.get("disease_name") or g_data.get("disease_name") or "Associated Phenotype"
    disease_id = r_plan.get("disease_id") or "MONDO:RECORD"
    clingen_status = r_plan.get("clingen_status") or "Disputed"
    finding_class = r_plan.get("finding_class") or "Clinical Inconsistency"

    sub_proto = r_plan.get("submitter_protocol", {})
    panel_gov = r_plan.get("panel_governance", {})
    safety_adv = r_plan.get("patient_safety", {})
    citations = r_plan.get("guidelines_citations", STANDARD_CITATIONS)

    # Build memo content
    lines = [
        "# CLINICAL REMEDIATION MEMORANDUM & DECISION-SUPPORT PROTOCOL",
        "**CuraCascade Clinical Governance & Genomic Audit Advisory Committee**  ",
        f"**Protocol Reference:** `{ref_id}`  ",
        f"**Date of Issuance:** {today_str}  ",
        "**Target Stakeholders:** Molecular Genetic Pathology Directors, Clinical Variant Curators, Medical Geneticists, PanelApp Administrators  ",
        f"**Subject:** **URGENT CLINICAL REMEDIATION DIRECTIVE: {sym} ({disease_name})**  ",
        "",
        "---",
        "",
        "## 1. Executive Summary & Authoritative ClinGen Determination",
        "",
        f"CuraCascade has identified critical downstream discordance between authoritative ClinGen Gene-Disease Validity appraisals and downstream repositories (NCBI ClinVar and Genomics England PanelApp) for **{sym}** in the context of **{disease_name}** (`{disease_id}`).",
        "",
        f"- **Gene Symbol:** **{sym}** ({gene_name}, HGNC: `{hgnc_id}`)",
        f"- **Authoritative ClinGen Classification:** **{clingen_status.upper()}**",
        f"- **Audit Taxonomy Finding:** **{finding_class}**",
        f"- **ClinGen Evidence Anchor:** {sub_proto.get('evidence_anchor', 'ClinGen Expert Panel')}",
        f"- **ClinGen Evidence Dossier:** [{sub_proto.get('report_url', 'ClinGen Dossier')}]({sub_proto.get('report_url', 'https://search.clinicalgenome.org/kb/gene-validity')})",
        "",
        "> **Clinical Rationale:**  ",
        f"> {sub_proto.get('rationale', 'Curation discordance detected requiring immediate resolution.')}",
        "",
        "---",
        "",
        "## 2. NCBI ClinVar Submitter Protocol: Variant Retraction & Downgrade Directives",
        "",
        f"### Action Required: **{sub_proto.get('action_title', 'Variant Reclassification')}**",
        f"- **Target Classification:** `{sub_proto.get('target_classification', 'VUS')}`",
        f"- **Governing Rule:** `{sub_proto.get('action_badge', 'ACMG Guidance')}`",
        "",
        "#### Procedural Checklist for Clinical Laboratories & Variant Curators:",
    ]

    for step in sub_proto.get("checklist", []):
        lines.append(f"- [ ] {step}")

    lines.extend([
        "",
        "#### ACMG / ClinGen Applicable Rules:",
    ])
    for rule in sub_proto.get("acmg_rules", []):
        lines.append(f"- • {rule}")

    lines.extend([
        "",
        "---",
        "",
        "## 3. Diagnostic Panel Maintenance: Genomics England PanelApp & Lab Filtering Directives",
        "",
        f"### Directive: **{panel_gov.get('action_title', 'Panel Update')}**",
        f"- **Target Panel Tier:** `{panel_gov.get('target_tier', 'Red Tier')}`",
        f"- **Filtering Directive:** {panel_gov.get('filtering_advisory', 'Exclude from routine diagnostic filtering.')}",
        "",
        "#### Bioinformatic & Panel Curator Directives:",
    ])

    for step in panel_gov.get("checklist", []):
        lines.append(f"- [ ] {step}")

    lines.extend([
        "",
        "---",
        "",
        "## 4. Patient Safety & Clinical Governance Directive: Retrospective Case Audit & Cascade Protocol",
        "",
        f"### Status: **{safety_adv.get('alert_badge', 'Clinical Alert')}**",
        f"**{safety_adv.get('advisory_title', 'Retrospective Clinical Review')}**",
        "",
        f"{safety_adv.get('guidance_summary', 'Immediate retrospective case review required.')}",
        "",
        "#### Mandatory Clinical Action Items for Medical Geneticists & Clinical Teams:",
    ])

    for act in safety_adv.get("action_items", []):
        lines.append(f"- [ ] **{act.split(':')[0]}:** {':'.join(act.split(':')[1:]).strip() if ':' in act else act}")

    lines.extend([
        "",
        "---",
        "",
        "## 5. Regulatory Compliance & Guidelines Citations",
        "",
        "This remediation protocol is derived deterministically in compliance with standard clinical genetics recommendations:",
        "",
    ])

    for cit in citations:
        lines.append(f"- **{cit['framework']}:** {cit['citation']}")

    lines.extend([
        "",
        "---",
        "",
        "*Memorandum generated automatically by CuraCascade Clinical Governance Engine.*  ",
        f"*Audit Timestamp:* `{r_plan.get('generated_at', datetime.datetime.now().isoformat())}`  ",
        "*System Authority: CuraCascade Cross-Layer Clinical Audit Architecture (3NF SQLite)*",
    ])

    return "\n".join(lines)


# ==============================================================================
# Clinical Governance Dossier Compiler
# ==============================================================================

def compile_clinical_governance_dossier(
    gene_symbol: str,
    gene_meta: Optional[Dict[str, Any]] = None,
    curations_df: Optional[pd.DataFrame] = None,
    variants_df: Optional[pd.DataFrame] = None,
    panels_df: Optional[pd.DataFrame] = None,
    findings_df: Optional[pd.DataFrame] = None,
    ancestry_profile: Optional[Dict[str, Any]] = None,
    remediation_plan: Optional[Dict[str, Any]] = None,
) -> str:
    """Compile a comprehensive Clinical Governance Dossier and Retraction Directive (.md).

    Generates a formal, regulatory-grade markdown report containing:
    1. Header: CLINICAL AUDIT & RETRACTION DIRECTIVE: <GENE_SYMBOL>
    2. Metadata Table: Evaluation date, ClinGen classification, MONDO ID, expert panel name.
    3. Orthogonal Population Genetics Evidence: gnomAD v4 global frequency, MCAF ceiling,
       fold-excess multiplier, and ancestral disparity statement.
    4. Downstream Affected Assertions: Table of all active ClinVar P/LP records and
       PanelApp diagnostic memberships.
    5. Prescriptive Remediation Checklist: Step-by-step ACMG/AMP reclassification protocol.
    6. Formal Clinical Sign-off Block: Signatures for Molecular Geneticist and
       Laboratory Medical Director.

    Args:
        gene_symbol: HGNC gene symbol (e.g. 'ANK2', 'SCN1B').
        gene_meta: Gene metadata dictionary (name, HGNC ID, etc.).
        curations_df: DataFrame of ClinGen curations for the gene.
        variants_df: DataFrame of ClinVar variants for the gene.
        panels_df: DataFrame of PanelApp entries for the gene.
        findings_df: DataFrame of audit findings flagged for the gene.
        ancestry_profile: Optional dictionary of ancestry breakdown metrics.
        remediation_plan: Optional remediation plan dictionary from generate_remediation_plan().

    Returns:
        Structured Markdown text for clinical governance archiving and laboratory compliance.
    """
    sym = str(gene_symbol).strip().upper()
    g_meta = gene_meta or {}
    gene_name = g_meta.get("gene_name") or f"{sym} Curated Human Gene"
    hgnc_id = g_meta.get("hgnc_id") or "HGNC:RECORD"

    now = datetime.datetime.now()
    today_str = now.strftime("%B %d, %Y")
    today_compact = now.strftime("%Y%m%d")

    # 1. Resolve Primary Clinical Curation Metadata
    eval_date = "N/A"
    clingen_status = "Disputed"
    mondo_id = "MONDO:RECORD"
    disease_name = "Associated Clinical Phenotype"
    expert_panel = "ClinGen Gene Curation Expert Panel"
    report_url = "https://search.clinicalgenome.org/kb/gene-validity"

    if curations_df is not None and not curations_df.empty:
        c_row = curations_df.iloc[0]
        eval_date = str(c_row.get("classification_date") or "N/A")
        clingen_status = str(c_row.get("classification") or "Disputed")
        mondo_id = str(c_row.get("disease_id") or "MONDO:RECORD")
        disease_name = str(c_row.get("disease_name") or "Associated Clinical Phenotype")
        expert_panel = str(c_row.get("gcep") or "ClinGen GCEP")
        report_url = str(c_row.get("report_url") or "https://search.clinicalgenome.org/kb/gene-validity")
    elif remediation_plan:
        eval_date = str(remediation_plan.get("classification_date") or "N/A")
        clingen_status = str(remediation_plan.get("clingen_status") or "Disputed")
        mondo_id = str(remediation_plan.get("disease_id") or "MONDO:RECORD")
        disease_name = str(remediation_plan.get("disease_name") or "Associated Clinical Phenotype")
        expert_panel = str(remediation_plan.get("gcep") or "ClinGen GCEP")

    # Audit Findings Summary
    finding_classes = "Consistent (No Inconsistencies)"
    if findings_df is not None and not findings_df.empty:
        finding_classes = ", ".join(sorted(findings_df["audit_class"].unique()))
    elif remediation_plan and remediation_plan.get("finding_class"):
        finding_classes = remediation_plan.get("finding_class")

    # 2. Resolve Population Genetics & Ancestral Discordance Evidence
    anc_prof = ancestry_profile or {}
    global_af = 0.0
    mcaf_ceiling = 0.0001
    top_var_name = "N/A"
    top_var_id = "N/A"

    if variants_df is not None and not variants_df.empty:
        v_sorted = variants_df.copy()
        if "population_refuted_flag" in v_sorted.columns:
            v_sorted = v_sorted.sort_values(
                by=["population_refuted_flag", "stars"], ascending=[False, False]
            )
        top_v = v_sorted.iloc[0]
        top_var_name = str(top_v.get("variant_name") or "Index Variant")
        top_var_id = str(top_v.get("variation_id") or "N/A")
        global_af = float(top_v.get("gnomad_af_global", 0.0) or 0.0)
        mcaf_ceiling = float(top_v.get("mcaf_threshold", 0.0001) or 0.0001)

    if anc_prof:
        global_af = float(anc_prof.get("af_global", global_af) or global_af)
        mcaf_ceiling = float(anc_prof.get("mcaf_threshold", mcaf_ceiling) or mcaf_ceiling)

    nfe_af = float(anc_prof.get("af_nfe", 0.0) or 0.0)
    afr_af = float(anc_prof.get("af_afr", 0.0) or 0.0)
    amr_af = float(anc_prof.get("af_amr", 0.0) or 0.0)
    eas_af = float(anc_prof.get("af_eas", 0.0) or 0.0)
    jpn_af = float(anc_prof.get("af_jpn", 0.0) or 0.0)
    disc_ratio = float(anc_prof.get("discordance_ratio", 1.0) or 1.0)
    risk_flag = int(anc_prof.get("ancestry_misclassification_risk", 0) or 0)

    # Continental cohort maximum
    non_euro_map = {
        "African / African American (AFR)": afr_af,
        "Admixed American / Latino (AMR)": amr_af,
        "East Asian Reference (EAS)": eas_af,
        "Japanese Cohort (JPN - ToMMo)": jpn_af,
    }
    highest_non_euro_cohort, highest_non_euro_af = max(
        non_euro_map.items(), key=lambda it: it[1]
    )

    max_observed_af = max(global_af, nfe_af, highest_non_euro_af)
    fold_excess = round(max_observed_af / max(1e-8, mcaf_ceiling), 1)

    # Formulation of Ancestral Disparity Statement
    if risk_flag == 1 or (highest_non_euro_af > mcaf_ceiling and nfe_af <= mcaf_ceiling):
        anc_disparity_statement = (
            f"CRITICAL DIAGNOSTIC EQUITY ALERT: The index variant demonstrates acute ancestral frequency discordance "
            f"(Ancestry Discordance Ratio = {disc_ratio:.1f}×). In historical European reference cohorts, the variant "
            f"appears rare (gnomAD NFE = {nfe_af:.2e} ≤ MCAF ceiling {mcaf_ceiling:.2e}), which historically led variant "
            f"curators to infer pathogenicity. However, in non-European ancestral reference cohorts, this variant is "
            f"observed at high population frequency ({highest_non_euro_cohort} AF = {highest_non_euro_af:.2e}, exceeding "
            f"the disease MCAF ceiling by {round(highest_non_euro_af / max(1e-8, mcaf_ceiling), 1):.1f}×). Clinical "
            f"laboratories relying on European-biased reference databases face an acute risk of issuing false-positive "
            f"pathogenic misdiagnoses to patients of non-European ancestry. Population frequency data provides independent, "
            f"orthogonal evidence that this variant represents a benign ancestral polymorphism under ACMG/AMP Rule BA1."
        )
    elif max_observed_af > mcaf_ceiling:
        anc_disparity_statement = (
            f"GLOBAL POPULATION REFUTATION: The variant allele frequency in gnomAD v4 ({max_observed_af:.2e}) "
            f"exceeds the calculated Maximum Credible Allele Frequency threshold (MCAF = {mcaf_ceiling:.2e}) by "
            f"{fold_excess:.1f}× across reference cohorts. Under the Whiffin et al. (2017) clinical genomics framework "
            f"and ACMG/AMP Rule BA1/BS1, this variant occurs at a population frequency incompatible with rare, "
            f"high-penetrance disease causality for {disease_name}. Population frequency alone provides definitive "
            f"orthogonal refutation of pathogenicity independent of historical case reports."
        )
    else:
        anc_disparity_statement = (
            f"CONCORDANT POPULATION RARITY: Observed allele frequencies across continental cohorts in gnomAD v4 remain "
            f"below the calculated MCAF ceiling ({mcaf_ceiling:.2e}), indicating population rarity. However, authoritative "
            f"ClinGen GCEP curation has established that the {sym}–{disease_name} relationship is classified as "
            f"{clingen_status.upper()}. In accordance with ACMG/ClinGen recommendations, variant-level population rarity "
            f"cannot establish pathogenicity in a gene lacking confirmed clinical validity."
        )

    # 3. Assemble Markdown Document
    dossier: List[str] = [
        f"# CLINICAL AUDIT & RETRACTION DIRECTIVE: {sym}",
        "",
        "**Document Authority:** CuraCascade Clinical Evidence Governance Framework  ",
        f"**Directive Protocol Identifier:** `DIRECTIVE-GOV-{sym}-{today_compact}`  ",
        f"**Effective Issuance Date:** {today_str}  ",
        "**Classification Standards:** ACMG/AMP (2015) • ClinGen Framework (2017/2022) • Whiffin et al. (2017)  ",
        "**Target Stakeholders:** Molecular Genetic Pathologists, Laboratory Medical Directors, Clinical Curators, Diagnostic Panel Committees  ",
        "",
        "---",
        "",
        "## 1. Clinical Evidence & Curation Metadata",
        "",
        "The following clinical curation metadata establishes the authoritative gene-disease validity determination appraised by expert consortia:",
        "",
        "| Clinical Metadata Field | Registry Value & Evidence Detail |",
        "| :--- | :--- |",
        f"| **Target Gene Symbol** | **{sym}** ({gene_name}) |",
        f"| **HGNC Accession** | `{hgnc_id}` |",
        f"| **Evaluation Date** | **{eval_date}** |",
        f"| **ClinGen Classification** | **{clingen_status.upper()}** |",
        f"| **Disease Indication** | {disease_name} |",
        f"| **MONDO Disease Ontology ID** | `{mondo_id}` |",
        f"| **Expert Panel Name (GCEP)** | **{expert_panel}** |",
        f"| **Authoritative Evidence Dossier** | [{report_url}]({report_url}) |",
        f"| **CuraCascade Audit Finding Class** | **{finding_classes}** |",
        "",
        "---",
        "",
        "## 2. Orthogonal Population Genetics Evidence",
        "",
        "### Statistical Frequency Evaluation against Disease Pathogenic Ceiling (MCAF)",
        "",
        "Under the Whiffin et al. (2017) clinical framework and ACMG/AMP BA1/BS1 rules, variants that occur in general population cohorts at frequencies exceeding the Maximum Credible Allele Frequency (MCAF) calculated from disease prevalence cannot plausibly cause high-penetrance rare genetic disease.",
        "",
        "| Population Genetics Parameter | Metric Value | Clinical Interpretation |",
        "| :--- | :--- | :--- |",
        f"| **Index Variant Evaluated** | `{top_var_name}` (ClinVar `{top_var_id}`) | Representative downstream assertion |",
        f"| **gnomAD v4 Global Frequency** | **{global_af:.2e}** ({global_af * 100:.4f}%) | General unselected population baseline |",
        f"| **Max Credible AF (MCAF Ceiling)** | **{mcaf_ceiling:.2e}** ({mcaf_ceiling * 100:.4f}%) | Calculated maximum credible pathogenic ceiling |",
        f"| **Fold-Excess Multiplier** | **{fold_excess:.1f}× MCAF** | Ratio of observed population frequency to pathogenic limit |",
        f"| **Population Genetics Verdict** | **{'🔴 Refuted (Exceeds MCAF)' if max_observed_af > mcaf_ceiling else '⚪ Rarity Compatible'}** | Standalone BA1/BS1 filtering criterion |",
        "",
        "### Multi-Ancestral Frequency Distribution Across Continental Cohorts",
        "",
        "| Continental Cohort | Region Category | Cohort Allele Frequency | Comparison to MCAF Ceiling |",
        "| :--- | :--- | :---: | :--- |",
        f"| **European Baseline (NFE)** | Historical Reference Baseline | {nfe_af:.2e} | {'🔴 Exceeds MCAF' if nfe_af > mcaf_ceiling else '🟢 Below MCAF'} |",
        f"| **African / African American (AFR)** | Underrepresented Cohort | {afr_af:.2e} | {'🔴 Exceeds MCAF' if afr_af > mcaf_ceiling else '🟢 Below MCAF'} |",
        f"| **Admixed American / Latino (AMR)** | Underrepresented Cohort | {amr_af:.2e} | {'🔴 Exceeds MCAF' if amr_af > mcaf_ceiling else '🟢 Below MCAF'} |",
        f"| **East Asian Reference (EAS)** | Continental Reference Cohort | {eas_af:.2e} | {'🔴 Exceeds MCAF' if eas_af > mcaf_ceiling else '🟢 Below MCAF'} |",
        f"| **Japanese Cohort (JPN - ToMMo)** | Continental Reference Cohort | {jpn_af:.2e} | {'🔴 Exceeds MCAF' if jpn_af > mcaf_ceiling else '🟢 Below MCAF'} |",
        "",
        "### Ancestral Disparity Statement",
        "",
        f"> **{anc_disparity_statement}**",
        "",
        "---",
        "",
        "## 3. Downstream Affected Assertions & Contradiction Inventory",
        "",
        "### Active NCBI ClinVar Pathogenic / Likely Pathogenic Records",
        "",
    ]

    # ClinVar Variants Table
    if variants_df is not None and not variants_df.empty:
        dossier.extend([
            "| ClinVar ID | Variant Name / cDNA | Clinical Significance | Stars | Review Status | Submitters | Last Evaluated | gnomAD AF | MCAF Ceiling | Population Status |",
            "| :---: | :--- | :--- | :---: | :--- | :---: | :---: | :---: | :---: | :--- |",
        ])
        for _, vr in variants_df.iterrows():
            stars_cnt = int(vr.get("stars", 0) or 0)
            star_glyphs = "★" * stars_cnt + "☆" * (4 - stars_cnt)
            v_af = float(vr.get("gnomad_af_global", 0.0) or 0.0)
            v_mcaf = float(vr.get("mcaf_threshold", 0.0001) or 0.0001)
            is_ref = vr.get("population_refuted_flag", 0) == 1 or v_af > v_mcaf
            pop_stat = f"🔴 {v_af / max(1e-8, v_mcaf):.1f}× MCAF" if is_ref else "🟢 Compatible"
            dossier.append(
                f"| [{vr['variation_id']}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{vr['variation_id']}/) "
                f"| {vr['variant_name']} "
                f"| **{vr['clinical_significance']}** "
                f"| {star_glyphs} "
                f"| {vr.get('review_status', 'N/A')} "
                f"| {vr.get('num_submitters', 1)} "
                f"| {vr.get('last_evaluated', 'N/A')} "
                f"| {v_af:.2e} "
                f"| {v_mcaf:.2e} "
                f"| {pop_stat} |"
            )
    else:
        dossier.append("*No active ClinVar assertions registered in the target database for this gene.*")

    dossier.extend([
        "",
        "### Genomics England PanelApp Diagnostic Panel Memberships",
        "",
    ])

    # PanelApp Panels Table
    if panels_df is not None and not panels_df.empty:
        dossier.extend([
            "| Panel ID | Panel Name | Version | Confidence Tier | Mode of Inheritance | Associated Phenotypes | Mandatory Governance Action |",
            "| :---: | :--- | :---: | :--- | :--- | :--- | :--- |",
        ])
        for _, pr in panels_df.iterrows():
            conf_raw = str(pr.get("confidence", "Green")).strip()
            conf_badge = f"🟢 {conf_raw.title()} (Diagnostic-Grade)" if conf_raw.lower() == "green" else f"🟡 {conf_raw.title()}"
            action_req = "Demote to Red (Contraindicated)" if conf_raw.lower() == "green" else "Demote to Amber / Review"
            dossier.append(
                f"| {pr['panel_id']} "
                f"| {pr['panel_name']} "
                f"| {pr.get('panel_version', '1.0')} "
                f"| {conf_badge} "
                f"| {pr.get('mode_of_inheritance', 'N/A')} "
                f"| {pr.get('phenotypes', 'N/A')} "
                f"| **{action_req}** |"
            )
    else:
        dossier.append("*No active PanelApp panel entries registered in the target database for this gene.*")

    dossier.extend([
        "",
        "---",
        "",
        "## 4. Prescriptive Remediation Checklist: Step-by-Step ACMG/AMP Reclassification Protocol",
        "",
        "Clinical laboratories, variant curators, and bioinformatics pipeline administrators must execute the following multi-phase remediation protocol:",
        "",
        "### Phase 1: Bioinformatic Pipeline Gating & Diagnostic Filtering",
        f"- [ ] **Step 1.1:** Immediately lock and flag **{sym}** in Laboratory Information Management Systems (LIMS) to halt automated pathogenic variant calling.",
        f"- [ ] **Step 1.2:** Exclude **{sym}** genomic coordinates from diagnostic target BED files for clinical exome and panel assays targeting `{disease_name}`.",
        f"- [ ] **Step 1.3:** Configure bioinformatic variant caller annotation rules to automatically flag any future assertion referencing **{sym}** with authoritative ClinGen status: `{clingen_status.upper()}`.",
        "",
        "### Phase 2: ACMG/AMP Reclassification Protocol",
        f"- [ ] **Step 2.1:** Execute formal variant-level re-evaluation under ACMG/AMP 2015 and ClinGen 2022 guidelines for all active {sym} assertions.",
        f"- [ ] **Step 2.2:** **Apply ACMG Rule BA1/BS1:** Where population allele frequency in gnomAD v4 exceeds the calculated disease MCAF threshold ({mcaf_ceiling:.2e}), assign Standalone Benign (BA1) or Strong Benign (BS1) criterion.",
        f"- [ ] **Step 2.3:** **Enforce Gene-Validity Ceiling:** For rare variants where population frequency does not trigger BA1, apply ACMG/ClinGen recommendations capping maximum classification at **Variant of Uncertain Significance (VUS)** based on {clingen_status.upper()} gene validity.",
        f"- [ ] **Step 2.4:** Formulate updated classification evidence summary citing official ClinGen GCEP determination dated `{eval_date}`.",
        "",
        "### Phase 3: NCBI ClinVar Submission Retraction & Registry Synchronization",
        f"- [ ] **Step 3.1:** Generate updated SCV submission batch in NCBI ClinVar Submitter Portal retracting Pathogenic/Likely Pathogenic status for all submitted variants.",
        f"- [ ] **Step 3.2:** Reclassify active variants to **VUS** (or Benign/Likely Benign if BA1 criteria met) with submitter comments referencing ClinGen Evidence Dossier: `{report_url}`.",
        f"- [ ] **Step 3.3:** Notify the ClinVar Curation Consortium of corrected SCV accessions.",
        "",
        "### Phase 4: Diagnostic Panel Maintenance (Genomics England PanelApp & In-House Panels)",
        f"- [ ] **Step 4.1:** Submit formal PanelApp curator ticket demoting **{sym}** from Green (diagnostic-grade) to **Red** (research-only / contraindicated) tier.",
        f"- [ ] **Step 4.2:** In-house diagnostic panels: Remove gene from routine reporting panels; notify ordering clinicians regarding updated test requisition forms.",
        "",
        "### Phase 5: Retrospective Patient Review & Cascade Testing Directive",
        f"- [ ] **Step 5.1:** Perform institutional query across clinical archives to identify all historical probands reported with positive {sym} findings since `{eval_date}`.",
        f"- [ ] **Step 5.2:** Re-issue amended diagnostic genetics reports to referring clinical geneticists and genetic counselors with updated classification.",
        f"- [ ] **Step 5.3:** **IMMEDIATE CASCADE TESTING HALT:** Issue urgent clinical notification instructing medical teams to halt all predictive cascade carrier testing and risk-reducing prophylactic interventions in asymptomatic relatives based solely on {sym} assertions.",
        "",
        "---",
        "",
        "## 5. Formal Clinical Governance Attestation & Regulatory Sign-Off",
        "",
        "```",
        "+====================================================================================================================+",
        "| FORMAL CLINICAL LABORATORY GOVERNANCE ATTESTATION & REGULATORY SIGN-OFF                                            |",
        "|                                                                                                                    |",
        "| I hereby attest that I have reviewed the multi-database evidence, authoritative ClinGen validity appraisal,       |",
        "| orthogonal population genetics parameters, and downstream contradiction inventory presented in this directive.    |",
        "|                                                                                                                    |",
        "| The prescriptive remediation checklist and reclassification protocols specified herein represent formal approved  |",
        "| institutional policy in accordance with ACMG/AMP, ClinGen, CAP, CLIA (42 CFR § 493), and ISO 15189 quality         |",
        "| management standards.                                                                                              |",
        "|                                                                                                                    |",
        "|                                                                                                                    |",
        "| ____________________________________________________      ____________________________________________________     |",
        "| Lead Clinical Molecular Geneticist, PhD, FACMG            Laboratory Medical Director, MD, FCAP / FRCPath          |",
        "| Board Certified Clinical Molecular Genetics               Director, Molecular Diagnostics & Genomic Pathology      |",
        "|                                                                                                                    |",
        "| License / NPI: _____________________________________      License / NPI: _____________________________________     |",
        "| Date:   ____ / ____ / ________                            Date:   ____ / ____ / ________                           |",
        "|                                                                                                                    |",
        "| Clinical Diagnostic Laboratory: __________________________________________________________________________________ |",
        "| Institutional CLIA / CAP / ISO Accreditation Identifier: _________________________________________________________ |",
        "+====================================================================================================================+",
        "```",
        "",
        "---",
        f"*Dossier compiled by CuraCascade Clinical Governance & Genomic Evidence Audit Engine at {now.isoformat()} UTC.*  ",
        "*System Provenance: CuraCascade 3NF Relational Cross-Database Clinical Infrastructure.*",
    ])

    return "\n".join(dossier)



# ==============================================================================
# 3. CLI Verification Entry Point
# ==============================================================================

def main():
    """CLI testing and demonstration entry point."""
    if sys.platform == "win32":
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except Exception:
            pass

    print("--- CuraCascade Objective 7: Clinical Remediation Engine Test ---\n")

    test_cases = [
        ("ANK2", "Class A", "Disputed", {"variation_id": 51902, "variant_name": "NM_001148.5(ANK2):c.4216G>A (p.Glu1406Lys)", "stars": 2, "clinical_significance": "Pathogenic", "gnomad_af_global": 0.0034, "mcaf_threshold": 0.00005, "population_refuted_flag": 1}),
        ("BARD1", "Class B", "Limited", {"variation_id": 127902, "variant_name": "NM_000465.4(BARD1):c.1664T>C (p.Leu555Pro)", "stars": 1, "clinical_significance": "Likely pathogenic"}),
        ("MYLK2", "Class A", "Disputed", {"variation_id": 14210, "variant_name": "NM_005597.4(MYLK2):c.355G>A (p.Glu119Lys)", "stars": 2, "clinical_significance": "Pathogenic"}),
        ("MYH7", "", "Definitive", {"variation_id": 14092, "variant_name": "NM_000257.4(MYH7):c.1208G>A (p.Arg403Gln)", "stars": 3, "clinical_significance": "Pathogenic"}),
    ]

    for sym, f_cls, status, v_row in test_cases:
        plan = generate_remediation_plan(
            gene_symbol=sym,
            finding_class=f_cls,
            clingen_status=status,
            variant_row=v_row,
        )
        print(f"=== {sym} ({status} | {f_cls}) ===")
        print(f"  Badge:               {plan['prescriptive_badge']}")
        print(f"  Submitter Action:    {plan['submitter_protocol']['action_title']}")
        print(f"  Panel Governance:    {plan['panel_governance']['action_title']}")
        print(f"  Patient Safety:      {plan['patient_safety']['advisory_title']}")
        print()

    # Generate test memo for ANK2
    ank2_plan = generate_remediation_plan("ANK2", "Class A", "Disputed", test_cases[0][3])
    memo = build_remediation_memo_markdown("ANK2", {"gene_name": "Ankyrin 2", "hgnc_id": "HGNC:493"}, ank2_plan)
    print(f"Sample Memo Generated for ANK2 ({len(memo)} chars). First 300 chars:")
    print(memo[:300])
    print("...\nDone!")


if __name__ == "__main__":
    main()
