"""Unit test for Clinical Governance Dossier compilation."""

from __future__ import annotations

from pipeline.audit import get_gene_history
from pipeline.remediation import compile_clinical_governance_dossier, generate_remediation_plan
from pipeline.allele_frequency import get_variant_ancestry_profile


def test_clinical_governance_dossier_all_sections():
    gene = "SCN1B"
    history = get_gene_history(gene)
    assert history, "History for SCN1B should not be empty"

    g_meta = history.get("gene", {})
    curations = history.get("curations")
    variants = history.get("variants")
    panels = history.get("panels")
    findings = history.get("findings")

    prof = None
    if not variants.empty:
        prof = get_variant_ancestry_profile(int(variants.iloc[0]["variation_id"]))

    rem_plan = generate_remediation_plan(
        gene_symbol=gene,
        finding_class="Class A, Class C",
        clingen_status="Refuted",
        variant_row=variants.iloc[0].to_dict() if not variants.empty else None,
        panel_row=panels.iloc[0].to_dict() if not panels.empty else None,
    )

    dossier = compile_clinical_governance_dossier(
        gene_symbol=gene,
        gene_meta=g_meta,
        curations_df=curations,
        variants_df=variants,
        panels_df=panels,
        findings_df=findings,
        ancestry_profile=prof,
        remediation_plan=rem_plan,
    )

    # 1. Header: CLINICAL AUDIT & RETRACTION DIRECTIVE: <GENE_SYMBOL>
    assert f"CLINICAL AUDIT & RETRACTION DIRECTIVE: {gene}" in dossier

    # 2. Metadata Table: Evaluation date, ClinGen classification, MONDO ID, expert panel name
    assert "Evaluation Date" in dossier
    assert "ClinGen Classification" in dossier
    assert "REFUTED" in dossier
    assert "MONDO:0015263" in dossier
    assert "ClinGen Brugada Syndrome GCEP" in dossier

    # 3. Orthogonal Population Genetics Evidence: gnomAD v4 global frequency, MCAF ceiling, fold-excess multiplier, ancestral disparity statement
    assert "gnomAD v4 Global Frequency" in dossier
    assert "Max Credible AF (MCAF Ceiling)" in dossier
    assert "Fold-Excess Multiplier" in dossier
    assert "Ancestral Disparity Statement" in dossier
    assert "European Baseline" in dossier

    # 4. Downstream Affected Assertions: Table of all active ClinVar P/LP records and PanelApp diagnostic memberships
    assert "Active NCBI ClinVar Pathogenic / Likely Pathogenic Records" in dossier
    assert "48312" in dossier  # SCN1B ClinVar variation ID
    assert "Genomics England PanelApp Diagnostic Panel Memberships" in dossier
    assert "Brugada syndrome" in dossier

    # 5. Prescriptive Remediation Checklist: Step-by-step ACMG/AMP reclassification protocol
    assert "Prescriptive Remediation Checklist: Step-by-Step ACMG/AMP Reclassification Protocol" in dossier
    assert "- [ ]" in dossier
    assert "Apply ACMG Rule BA1/BS1" in dossier
    assert "NCBI ClinVar Submission Retraction" in dossier

    # 6. Formal Clinical Sign-off Block: Signatures for Molecular Geneticist and Laboratory Medical Director
    assert "Formal Clinical Governance Attestation & Regulatory Sign-Off" in dossier
    assert "Lead Clinical Molecular Geneticist" in dossier
    assert "Laboratory Medical Director" in dossier
