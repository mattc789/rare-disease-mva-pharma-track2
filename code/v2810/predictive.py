"""V2.8.6 patient-first, scenario-based pharmacological mechanism predictor.

Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch

Research triage only: this is NOT a clinically calibrated efficacy model.
The score is neither a probability nor a predicted improvement in symptoms.
It does not change the evidence grades A/B/C/D produced by v285.evidence.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import csv
import hashlib
import json
import math
from pathlib import Path
import re

VERSION = 'V288_PATIENT_FIRST_MECHANISM_2'

CORE_CELL_CYCLE = frozenset({'CDK1', 'CDK2', 'CDK5', 'PLK1', 'AURKA', 'AURKB', 'CDC20', 'BUB1B', 'CHEK1'})
SEVERE_LOF = frozenset({'premature_stop', 'frameshift'})
EFFECT_LABELS = {
    'premature_stop': 'premature stop / possible transcript decay or truncation',
    'frameshift': 'reading-frame disruption; downstream protein function uncertain',
    'splice_altering': 'candidate aberrant RNA splicing or isoform expression',
    'missense': 'candidate folding, activity or protein-interaction defect',
    'inframe_indel': 'possible altered protein stability, activity or interaction',
    'uncertain_effect': 'mechanism not established',
}
EVIDENCE_TIER_MEANINGS = {
  'A': 'Provisional published functional disease-model rescue evidence from at least two separate sources; no patient benefit is proven.',
  'B': 'Direct functional drug-target data or one candidate-specific functional experiment; patient rescue remains unconfirmed.',
  'C': 'Mechanism of action directly connected to a target, without proved functional disease rescue.',
  'D': 'Indirect network, transcriptomic or chemical-similarity hypothesis only.',
}
PREDICTION_LABELS = {
    'DIRECTIONALLY_CONSISTENT': 'Drug action is consistent with an allele-informed research objective, but functional patient rescue is unproven.',
    'HYPOTHESIS_ONLY': 'There is an indirect mechanistic hypothesis; no reliable patient-specific benefit has been demonstrated.',
    'CONFLICTING_DIRECTION': 'The proposed action may worsen the assumed gene-function deficit, or a critical mechanism conflicts with the proposed rescue.',
    'CASE_LINK_UNRESOLVED': 'There is no sufficiently attributable drug-to-affected-gene pathway in the available research model.',
    'INSUFFICIENT_EVIDENCE': 'There is no suitable evidence to estimate even directional mechanism compatibility.',
    'NOT_ELIGIBLE_REGISTRY': 'No eligible single-ingredient medicine approval could be confirmed in the pinned registry snapshot.',
}
# Explicitly transparent, uncalibrated heuristic weighting. These are *not* scientific priors.
WEIGHTS = {'gene_link_direct': 14.0, 'gene_link_indirect': 5.0,
           'variant_fit_direct': 13.0, 'variant_fit_indirect': 5.0,
           'measured_functional_assay': 17.0, 'direct_pharmacological_moa': 13.0,
           'network_1_hop': 6.0, 'network_2_hop': 2.5, 'network_3_hop': 0.5,
           'transcriptomic_expression': 1.5, 'structural_similarity': 0.5,
           'evidence_grade_A': 13.0, 'evidence_grade_B': 8.0,
           'evidence_grade_C': 2.0, 'curated_exposure': 7.0,
           'high_risk_cell_cycle': -7.0, 'mechanism_uncertainty': -4.0,
           'missing_residual_function': -2.0, 'phase_unresolved': -1.0,
           'unsupported_context': -4.0}


def _norm(value):
    return re.sub(r'\s+', ' ', re.sub(r'[^A-Za-z0-9]+', ' ', str(value or '').upper())).strip()


def _realnumber(value):
    try:
        f = float(value)
        return f if math.isfinite(f) and f > 0 else None
    except (TypeError, ValueError, OverflowError):
        return None


def _clamp(value, minimum=0.0, maximum=100.0):
    return round(max(minimum, min(maximum, float(value))), 3)


def _registry_default(candidate, registry_index):
    # The function is intentionally imported lazily: independent fake-data tests
    # do not depend on FDA or the connected a specific development workstation workspace.
    from . import registry_output
    return registry_output.for_candidate(candidate, registry_index or {})


def registry_eligibility(approval):
    """Require an actual single-ingredient register listing; combos stay excluded."""
    a = approval if isinstance(approval, dict) else {}
    state = str(a.get('overall_state') or '').lower()
    records = a.get('records') or []
    listing = state == 'single_ingredient_product_listed' or a.get('single_product_listed') is True
    positive_records = []
    for record in records:
        context = str(record.get('product_context') or '').lower()
        status = str(record.get('status') or '').lower()
        if context != 'single_ingredient_product':
            continue
        if any(term in status for term in ('withdrawn', 'revoked', 'suspended', 'cancelled', 'expired', 'rejected')):
            continue
        positive_records.append(record)
    eligible = bool(listing and positive_records)
    return {
        'eligible_for_research_search': eligible,
        'registry_state': state or 'unknown',
        'regulator_names': sorted(set(a.get('registry_names') or [])),
        'jurisdictions': sorted(set(a.get('jurisdictions') or [])),
        'approved_specific_forms': sorted({str(record.get('specific_ingredient_form') or record.get('component') or '')
                                           for record in positive_records if record.get('specific_ingredient_form') or record.get('component')}),
        'qualifying_register_records': len(positive_records),
        'marketed_now_verified': False,
        'indication_for_this_genetic_disease_verified': False,
        'warning': ('Single-ingredient register listing in a pinned snapshot is not proof of current availability, '
                    'pediatric indication, safe exposure or approval for the patient’s disease.'),
    }


def _variant_profiles(patient_model):
    profiles = {}
    for p in patient_model.get('primary_gene_hypotheses', []):
        gene = str(p.get('gene') or '').upper().strip()
        if not gene:
            continue
        profiles[gene] = {
            'gene': gene,
            'effects': {str(v.get('effect_class') or 'uncertain_effect')
                        for v in p.get('variant_effect_models', [])},
            'variants': list(p.get('variant_effect_models', [])),
            'phase': (p.get('phase') or {}).get('state') or 'unresolved',
            'measured_disease_mechanism': p.get('actual_disease_mechanism_confirmed') is True,
            'objectives': list(p.get('functional_objectives') or []),
        }
    return profiles


def _route_gene_link(route, models):
    """No invented gene associations when multiple genes are possible."""
    if not models:
        return None, 'no_genetic_hypothesis', 0.0
    target = _norm(route.get('target_concept'))
    targets = {_norm(x.get('gene')) for x in (route.get('targets') or [])}
    family = str(route.get('family') or '')
    for gene in models:
        if target == gene or (gene in targets and family == 'direct_gene'):
            return gene, 'direct_gene_target', WEIGHTS['gene_link_direct']
    source_gene=str(route.get('case_gene') or '').upper()
    if source_gene in models:
        return source_gene,'public_route_case_gene_provenance',WEIGHTS['gene_link_indirect']
    text = str(route.get('mechanism_concept') or '')
    matches = [gene for gene in models if re.search(r'(?<![A-Za-z0-9])' + re.escape(gene) + r'(?![A-Za-z0-9])', text, re.I)]
    if len(matches) == 1:
        return matches[0], 'gene_explicitly_named_in_route_hypothesis', WEIGHTS['gene_link_indirect']
    if len(models) == 1 and family in ('literature_discovered', 'signed_upstream', 'proteostasis',
                                      'variant_rescue', 'splicing', 'direct_gene'):
        gene = next(iter(models))
        return gene, 'provisional_single_gene_association_not_direct_target', 1.5
    return None, 'not_attributable_to_specific_gene', 0.0


def _functional_match(route, profile):
    """Mutation-specific *compatibility hypotheses*, not known rescue biology."""
    effects = profile['effects'] if profile else set()
    family = str(route.get('family') or '')
    variable = str(route.get('variable') or '').lower()
    target = str(route.get('target_concept') or '').lower()
    direction = str(route.get('desired_action') or '').lower()
    issues = []
    if not effects:
        return 0.0, 'genetic_effect_not_interpretable', ['no_variant_effect_class']
    if direction not in ('increase', 'decrease'):
        return 0.0, 'mechanistic_direction_unknown', ['drug_action_direction_unknown']
    if direction == 'decrease' and effects.intersection(SEVERE_LOF | {'missense', 'splice_altering', 'inframe_indel'}):
        if family in ('direct_gene', 'proteostasis', 'splicing', 'variant_rescue'):
            return 0.0, 'action_conflicts_with_possible_loss_of_function', ['decreasing_defective_gene_function_not_a_rescue']
    if family == 'variant_rescue':
        if 'readthrough' in target or 'full_length_translation' in variable:
            if 'premature_stop' not in effects:
                return 0.0, 'readthrough_not_applicable_without_nonsense_variant', ['readthrough_not_general_frameshift_rescue']
            return 13.0, 'nonsense_stop_readthrough_only_if_transcript_and_frame_permit', [
                'premature_stop_context_and_NMD_status_unverified']
        if 'nmd' in target or 'nmd' in variable:
            if not effects.intersection(SEVERE_LOF):
                return 0.0, 'nmd_therapy_not_genotypically_supported', ['NMD_mechanism_unconfirmed']
            return 3.0, 'NMD_hypothesis_not_functional_rescue', ['stabilized_mutant_transcript_may_be_nonfunctional']
    if family == 'splicing':
        if 'splice_altering' in effects:
            return 13.0, 'splice_hypothesis_relevant', ['specific_exon_and_transcript_response_unverified']
        return 0.0, 'no_splice_variant_evidence', ['RNA_isoform_dysfunction_not_established']
    if family == 'proteostasis':
        if effects.intersection({'missense', 'inframe_indel'}):
            return 13.0, 'missense_proteostasis_plausible_if_residual_protein_exists', [
                'protein_folding_and_residual_abundance_unmeasured']
        return 0.0, 'no_clear_residual_protein_substrate', ['protein_stabilization_requires_residual_function']
    if family == 'direct_gene':
        if 'activity' in variable or 'complex_function' in variable:
            if effects.intersection({'missense', 'inframe_indel'}):
                return 13.0, 'possible_missense_activity_rescue', ['activity_increase_not_measured_in_patient']
            return 2.0, 'protein_activity_may_not_survive_severe_variants', ['residual_protein_function_unknown']
        if 'abundance' in variable or 'folding' in variable:
            if effects.intersection({'missense', 'inframe_indel', 'splice_altering'}):
                return 9.0, 'possibly_restorable_abundance_conditional_on_residual_protein', [
                    'increasing_nonfunctional_protein_does_not_restore_function']
            if effects.intersection(SEVERE_LOF):
                return 1.0, 'protein_abundance_not_equivalent_to_full_length_translation', [
                    'truncated_or_null_protein_may_be_incapable_of_rescue']
        return 5.0, 'gene_targeting_relevant_but_specific_effect_unknown', ['allele_function_not_established']
    if family == 'signed_upstream':
        return 5.0, 'conditional_upstream_compensation', ['network_sign_not_disease_specific',
                                                         'upstream_modulation_may_have_many_off_targets']
    if family == 'literature_discovered':
        return 5.0, 'literature_molecular_hypothesis_only', ['paper_model_may_not_transfer_to_patient_variant']
    return 1.0, 'general_indirect_disease_association', ['functional_connection_unverified']


def _signal_strength(pair):
    signs = pair.get('signals') or []
    counts = Counter(s.get('provider') or 'unknown' for s in signs)
    hop = min((int(s.get('depth')) for s in signs
               if s.get('provider') == 'ChEMBL_MoA_SIGNOR' and s.get('depth') is not None), default=None)
    measured = [s for s in signs if s.get('provider') == 'ChEMBL_measured_activity'
                and s.get('direct_target_assay_verified') is True and _realnumber(s.get('value_nM')) is not None]
    moa = pair.get('route') or {}
    weights = []
    if hop is not None:
        weights.append(WEIGHTS.get('network_' + str(hop) + '_hop', 0.0) if hop > 0
                       else WEIGHTS['direct_pharmacological_moa'])
    if measured:
        weights.append(WEIGHTS['measured_functional_assay'])
    if counts['perturbational_transcriptomics']:
        weights.append(WEIGHTS['transcriptomic_expression'])
    if counts['rdkit_morgan_similarity']:
        weights.append(WEIGHTS['structural_similarity'])
    return {
        'network_hops': hop,
        'direct_measured_assays': len({s.get('assay_id') for s in measured if s.get('assay_id')}) or len(measured),
        'measured_potency_nM': min((_realnumber(s.get('value_nM')) for s in measured), default=None),
        'channel_counts': dict(counts),
        'independently_verified_disease_rescue': False,
        'signal_component': min(sum(weights), 32.0),
    }


def _load_exposure(path):
    """Opt-in human-curated PK only. Never retrieve or invent plasma exposures."""
    if not path:
        return {}
    file = Path(path)
    if not file.exists():
        return {}
    data = json.loads(file.read_text(encoding='utf-8'))
    if data.get('schema') != 'V286_CURATED_PK_1':
        raise ValueError('Curated PK file has unsupported schema')
    out = {}
    for row in data.get('entries') or []:
        compound = _norm(row.get('compound'))
        c = _realnumber(row.get('unbound_cmax_nM'))
        if (not compound or not c or str(row.get('data_type')) != 'measured_unbound_plasma_cmax'
            or str(row.get('units')) != 'nM' or not str(row.get('source_reference') or '').strip()
            or not str(row.get('population') or '').strip() or row.get('source_checked') is not True):
            continue
        out.setdefault(compound, []).append({**row, 'unbound_cmax_nM': c})
    return out


def _exposure_estimate(pair, pk_index):
    labels = pair.get('member_labels') or [pair.get('ingredient')]
    eligible = [record for label in labels for record in pk_index.get(_norm(label), [])]
    assay = [s for s in pair.get('signals') or []
             if s.get('provider') == 'ChEMBL_measured_activity'
             and s.get('direct_target_assay_verified') is True
             and _realnumber(s.get('value_nM')) is not None]
    if not eligible:
        return {'state': 'HUMAN_UNBOUND_EXPOSURE_UNKNOWN',
                'functional_target_engagement': None,
                'clinical_or_paediatric_exposure': 'not_estimated'}
    if not assay:
        return {'state': 'RELEVANT_DIRECT_FUNCTIONAL_POTENCY_UNKNOWN',
                'functional_target_engagement': None,
                'reference_population': eligible[0].get('population'),
                'clinical_or_paediatric_exposure': 'not_estimated'}
    record = eligible[0]
    c = record['unbound_cmax_nM']
    potency = min(_realnumber(s.get('value_nM')) for s in assay)
    if potency is None:
        return {'state': 'FUNCTIONAL_POTENCY_UNKNOWN', 'functional_target_engagement': None}
    # Hill=1 *illustration* in which apparent Emax fraction = C/(C+EC50).
    # It is NOT a patient-level model of target engagement or treatment effect.
    frac = c / (c + potency)
    low = (c / 3) / (c / 3 + 3 * potency)
    high = (3 * c) / (3 * c + potency / 3)
    return {'state': 'ILLUSTRATIVE_FUNCTIONAL_ASSAY_VS_HUMAN_PK_CONTEXT',
            'functional_target_engagement': round(frac, 4),
            'sensitivity_low_high_illustrative': [round(low, 4), round(high, 4)],
            'model': 'nominal_C_over_C_plus_EC50_Hill1_sensitivity_threefold_not_CI',
            'assay_potency_nM': potency, 'reported_free_plasma_cmax_nM': c,
            'reference_population': record['population'],
            'source_reference': record['source_reference'],
            'clinical_or_paediatric_exposure': 'not_estimated',
            'limitations': ['Measured unbound plasma exposure may not equal tissue or intracellular exposure',
                            'EC50/IC50 response is assay-specific, not whole-patient pharmacodynamics',
                            'Dose, age, distribution and protein binding cannot be extrapolated to this child',
                            'Prediction does not estimate restoration of chromosome segregation or clinical benefit']}


def predict_pair(pair, patient_model, approval, case_context=None, pk_index=None):
    """One auditable and patient-conditioned mechanistic counterfactual."""
    p = pair
    route = p.get('route') or {}
    models = _variant_profiles(patient_model)
    gene, link, link_score = _route_gene_link(route, models)
    model = models.get(gene) if gene else None
    issues = []
    compatibility, match_type, match_issues = _functional_match(route, model)
    issues.extend(match_issues)
    regulatory = registry_eligibility(approval)
    signal = _signal_strength(p)
    grade = str(p.get('tier') or 'D')
    state = str(route.get('state') or '').lower()
    if grade in ('A', 'B', 'C'):
        score_grade = WEIGHTS['evidence_grade_' + grade]
    else:
        score_grade = 0.0
    cell_cycle_targets = sorted({str(x.get('gene') or '').upper()
                                 for x in route.get('targets') or []} & CORE_CELL_CYCLE)
    if cell_cycle_targets:
        issues.append('cell_cycle_target_requires_chromosome_instability_and_oncologic_risk_review')
    if state == 'hold' or (route.get('scientific_gate') or {}).get('hard_hold'):
        issues.append('causal_route_on_HOLD_not_predictively_eligible')
    if model and model.get('phase') not in ('trans_supported_by_vcf', 'homozygous_observed'):
        issues.append('allelic_phase_not_confirmed_as_recessive_causal_evidence')
    if model and not model.get('measured_disease_mechanism'):
        issues.append('patient_variant_function_has_not_been_measured')
    patient_residual = (case_context or {}).get('residual_protein_status', 'unknown')
    if str(patient_residual).lower() in ('unknown', 'not_provided', 'none'):
        issues.append('residual_functional_protein_unknown')
    exposure = _exposure_estimate(pair, pk_index or {})
    if exposure['state'] != 'ILLUSTRATIVE_FUNCTIONAL_ASSAY_VS_HUMAN_PK_CONTEXT':
        issues.append('human_unbound_exposure_or_direct_assay_unavailable')
    else:
        issues.append('illustrative_pk_assay_context_not_patient_specific')
    sci_issues = (route.get('scientific_gate') or {}).get('issues') or []
    issues.extend(str(x) for x in sci_issues)
    direct = signal['direct_measured_assays'] > 0
    signal_score = signal['signal_component']
    evidence_grade_score = score_grade
    exposure_score = WEIGHTS['curated_exposure'] if exposure['state'].startswith('ILLUSTRATIVE_') else 0.0
    uncertain = 0.0
    if patient_residual in (None, 'unknown', 'not_provided'):
        uncertain += WEIGHTS['missing_residual_function']
    if model and model.get('phase') in ('unresolved', 'not_assessable'):
        uncertain += WEIGHTS['phase_unresolved']
    if link.startswith('provisional_'):
        uncertain += WEIGHTS['mechanism_uncertainty']
    if sci_issues:
        uncertain += WEIGHTS['unsupported_context']
    risk_score = WEIGHTS['high_risk_cell_cycle'] if cell_cycle_targets else 0.0
    score_raw = link_score + compatibility + signal_score + evidence_grade_score + exposure_score + risk_score + uncertain
    # No weak evidence or missing phase can be made into a high-confidence
    # prediction by adding many signatures or 3-step network paths.
    score = _clamp(score_raw)
    plausible_gene = gene is not None
    if not regulatory['eligible_for_research_search']:
        hypothesis_status = 'NOT_ELIGIBLE_REGISTRY'
        score = 0.0
    elif state == 'hold' or (route.get('scientific_gate') or {}).get('hard_hold'):
        hypothesis_status = 'CONFLICTING_DIRECTION'
        score = 0.0
    elif 'conflicts_with' in match_type or 'not_applicable' in match_type:
        hypothesis_status = 'CONFLICTING_DIRECTION'
        score = 0.0
    elif not plausible_gene:
        hypothesis_status = 'CASE_LINK_UNRESOLVED'
        score = min(score, 5.0)
    elif not signal['channel_counts']:
        hypothesis_status = 'INSUFFICIENT_EVIDENCE'
        score = min(score, 5.0)
    elif (link == 'direct_gene_target' and direct and compatibility >= 5
          and not sci_issues):
        hypothesis_status = 'DIRECTIONALLY_CONSISTENT'
        score = min(score, 65.0)
    else:
        hypothesis_status = 'HYPOTHESIS_ONLY'
        score = min(score, 45.0)
    # Two explicit sensitivity scenarios; not statistical credibility intervals.
    scenario_conservative = _clamp(score_raw - signal_score + (17.0 if direct else 0.0)
                                   - max(0, signal_score - 1.5) * 0.3)
    scenario_exploratory = _clamp(score_raw + min(4.0, max(0.0, signal_score) * 0.3))
    if hypothesis_status in ('NOT_ELIGIBLE_REGISTRY', 'CONFLICTING_DIRECTION'):
        scenario_conservative = 0.0
        scenario_exploratory = 0.0
    if hypothesis_status == 'CASE_LINK_UNRESOLVED':
        scenario_conservative = min(scenario_conservative, 5.0)
        scenario_exploratory = min(scenario_exploratory, 5.0)
    if hypothesis_status == 'HYPOTHESIS_ONLY':
        scenario_conservative = min(scenario_conservative, 45.0)
        scenario_exploratory = min(scenario_exploratory, 45.0)
    if hypothesis_status == 'DIRECTIONALLY_CONSISTENT':
        scenario_conservative = min(scenario_conservative, 65.0)
        scenario_exploratory = min(scenario_exploratory, 65.0)
    variant_refs = [v.get('variant_ref') for v in (model or {}).get('variants',[]) if v.get('variant_ref')]
    effects = sorted((model or {}).get('effects') or [])
    if not regulatory['eligible_for_research_search']:
        issues.append('not_confirmed_as_authorised_single_ingredient_product')
    if grade in ('A','B'):
        issues.append('evidence_grade_is_research_strength_not_patient_efficacy')
    result = {
        'schema': VERSION,
        'candidate': p.get('ingredient'),
        'scientific_group_id': p.get('scientific_group_id'),
        'member_labels': p.get('member_labels') or [p.get('ingredient')],
        'route_id': p.get('route_id'),
        'route_target': route.get('target_concept'),
        'route_family': route.get('family'),
        'desired_action': route.get('desired_action'),
        'primary_gene_hypothesis': gene,
        'affected_variant_refs': variant_refs,
        'affected_variant_effect_classes': effects,
        'phase_assessment': (model or {}).get('phase', 'unknown'),
        'gene_association_method': link,
        'functional_match_model': match_type,
        'mechanistic_priority_score': score,
        'score_is_calibrated_efficiency_probability': False,
        'score_role': 'uncalibrated_ordinal_research_triage_only',
        'conservative_mechanistic_scenario': scenario_conservative,
        'exploratory_mechanistic_scenario': scenario_exploratory,
        'scenario_span': round(abs(scenario_exploratory - scenario_conservative), 3),
        'prediction_class': hypothesis_status,
        'prediction_explanation': PREDICTION_LABELS[hypothesis_status],
        'research_evidence_tier': grade,
        'why_this_evidence_tier': EVIDENCE_TIER_MEANINGS.get(grade, 'Research evidence not classified'),
        'scientific_evidence_blockers': sorted(set(str(x) for x in (p.get('evidence_blockers') or []))),
        'research_evidence_tier_unchanged': True,
        'source_signal_summary': signal,
        'registry': regulatory,
        'exposure_target_assay_context': exposure,
        'clinical_suitability_verified': False,
        'predicted_probability_of_clinical_benefit': None,
        'predicted_degree_of_patient_rescue': None,
        'potential_risks_and_gaps': sorted(set(issues)),
        'mechanistic_risk_flags': cell_cycle_targets,
        'next_discriminating_experiments': _experiments(effects, route, signal, exposure),
        'score_components': {
            'gene_attribution': round(link_score, 2),
            'allele_mechanism_alignment': round(compatibility, 2),
            'pharmacology_and_network': round(signal_score, 2),
            'research_grade': round(evidence_grade_score, 2),
            'curated_pk_assay_context': round(exposure_score, 2),
            'safety_caution': round(risk_score, 2),
            'uncertainty_penalty': round(uncertain, 2),
        },
        'patient_sequence_included': False,
        'no_personal_genotype_sent_to_AI': True,
    }
    return result


def _experiments(effects, route, signal, exposure):
    tests = [
        'Verify VCF allelic phase and clinically interpret each implicated variant (ACMG/ClinGen).',
        'Measure disease-relevant protein abundance, activity and cellular phenotype in appropriate patient-derived cells.',
        'Measure candidate-specific dose-response in disease-relevant cells and check target expression.',
    ]
    if 'premature_stop' in effects:
        tests.append('Verify transcript abundance, NMD behavior and possibility of full-length translation for the nonsense allele.')
    if 'frameshift' in effects:
        tests.append('Confirm frameshift transcript/protein consequences; generic nonsense readthrough is not presumed useful.')
    if 'missense' in effects or 'inframe_indel' in effects:
        tests.append('Determine residual stability, localization and interacting partners of each missense/in-frame allele.')
    if 'splice_altering' in effects:
        tests.append('Quantify aberrant isoforms by RNA assays and assess productive splicing under perturbation.')
    if signal['network_hops'] is not None and signal['network_hops'] > 0:
        tests.append('Test every signed causal intermediate rather than inferring rescue from pathway proximity.')
    if exposure['state'] != 'ILLUSTRATIVE_FUNCTIONAL_ASSAY_VS_HUMAN_PK_CONTEXT':
        tests.append('Establish target-relevant potency and verified unbound human exposure before evaluating feasibility.')
    tests.append('Review genotype-dependent toxicity, oncologic/chromosomal-stability effects and appropriate pediatric pharmacology.')
    return tests


def _safe_summary(p):
    """Short result object for ranking and pipeline, patient alleles only hashed."""
    return {k:p.get(k) for k in (
        'candidate','scientific_group_id','route_id','research_evidence_tier',
        'primary_gene_hypothesis','affected_variant_refs','prediction_class',
        'mechanistic_priority_score','conservative_mechanistic_scenario',
        'exploratory_mechanistic_scenario','score_is_calibrated_efficiency_probability',
        'potential_risks_and_gaps')}


def rank_pairs(pairs, patient_model, registry_index, case_context=None,
               *, resolver=None, curated_pk_path=None, exclude_unregistered=True):
    """Score all pairs, then pass *only register-eligible* pairs to drug review.

    Excluded rows remain visible in the predictor audit summary; no evidence
    or clinical grade is silently upgraded. No external HTTP/API calls here.
    """
    resolve = resolver or _registry_default
    pk = _load_exposure(curated_pk_path)
    eligible = []
    rejected = []
    counts = Counter()
    for pair in pairs:
        approval = resolve(pair.get('ingredient'), registry_index)
        result = predict_pair(pair, patient_model, approval, case_context, pk)
        row = deepcopy(pair)
        row['patient_first_score'] = result['mechanistic_priority_score']
        row['patient_first_prediction'] = _safe_summary(result)
        row['patient_first_prediction_class'] = result['prediction_class']
        row['patient_first_hypothesis_status'] = result['functional_match_model']
        counts[result['prediction_class']] += 1
        if result['registry']['eligible_for_research_search'] or not exclude_unregistered:
            eligible.append(row)
        else:
            rejected.append({'candidate':pair.get('ingredient'),
                             'reason':'No verified single-ingredient regulatory listing',
                             'registry_state':result['registry']['registry_state']})
    eligible.sort(key=lambda r:(-r['patient_first_score'],
                                {'A':0,'B':1,'C':2,'D':3}.get(r.get('tier'),4),
                                -(float(r.get('score') or 0)),
                                str(r.get('ingredient') or '')))
    summary = {
        'schema': VERSION,
        'all_pairs_received': len(pairs),
        'regulator_eligible_route_drug_pairs': len(eligible),
        'excluded_unregistered_or_combo_pairs': len(rejected),
        'prediction_class_counts': dict(counts),
        'prediction_score_type': 'heuristic_research_triage_not_probability',
        'patient_genotypes_exported': False,
        'clinical_prediction_validated': False,
        'functional_proof_required_for_patient_candidacy': True,
        'curated_pk_entries_loaded': sum(map(len, pk.values())),
        'top_candidate_route_pairs': [r['patient_first_prediction'] for r in eligible[:15]],
        'excluded_examples': rejected[:20],
    }
    return eligible, summary


def build_reports(path, grouped, pairs, patient_model, registry_index, case_context=None,
                  safety_by_name=None, *, resolver=None, curated_pk_path=None, top_n=30):
    """Write auditable, machine-readable and Markdown counterfactual predictions."""
    out = Path(path)
    out.mkdir(parents=True, exist_ok=False)
    resolve = resolver or _registry_default
    pk = _load_exposure(curated_pk_path)
    safety_by_name = safety_by_name or {}
    prediction_rows = []
    by_drug = defaultdict(list)
    for pair in pairs:
        approved = resolve(pair.get('ingredient'), registry_index)
        result = predict_pair(pair, patient_model, approved, case_context, pk)
        result['expert_triage']=pair.get('expert_triage') or {}
        result['expert_investigation_bucket']=(pair.get('expert_triage') or {}).get('investigation_bucket','T4')
        if not result['registry']['eligible_for_research_search']:
            continue
        src_safety = safety_by_name.get(pair.get('ingredient')) or {}
        result['public_label_safety_state'] = src_safety.get('state', 'NOT_ASSESSED')
        result['public_label_safety_flags'] = sorted(src_safety.get('flags',[]) or [])
        if src_safety.get('state') in ('HOLD', 'CLOSED'):
            result['potential_risks_and_gaps'] = sorted(set(result['potential_risks_and_gaps']) | {
                'public_product_label_additional_safety_review_required'})
        # A separate *research follow-up* priority reflects regulatory label warnings.
        # It does NOT rewrite the mechanistic score or constitute risk-benefit clearance.
        safety_penalty = (10.0 if 'boxed_warning_requires_review' in result['public_label_safety_flags']
                          else 5.0 if result['public_label_safety_flags'] else 0.0)
        result['research_followup_safety_penalty'] = safety_penalty
        result['research_followup_priority'] = _clamp(result['mechanistic_priority_score'] - safety_penalty)
        by_drug[result['scientific_group_id'] or result['candidate']].append(result)
        prediction_rows.append(result)
    # A single approved compound can modulate multiple routes, in opposing directions.
    all_best = []
    for key, routes in by_drug.items():
        routes.sort(key=lambda p:(-p['mechanistic_priority_score'],p['route_id'] or ''))
        best = routes[0]
        actions = defaultdict(set)
        for r in routes:
            if r.get('route_target'):
                actions[str(r['route_target']).upper()].add(str(r.get('desired_action') or ''))
        conflicting = sorted(target for target, directions in actions.items()
                             if 'increase' in directions and 'decrease' in directions)
        best['alternative_route_count'] = len(routes) - 1
        best['cross_route_direction_conflicts'] = conflicting
        if conflicting:
            best['potential_risks_and_gaps'] = sorted(set(best['potential_risks_and_gaps']) | {
                'cross_route_action_direction_conflict_requires_resolution'})
        all_best.append(best)
    all_best.sort(key=lambda p:(
                                {'T1':0,'T2':1,'T3':2,'T4':3,'HOLD':4}.get(
                                  p.get('expert_investigation_bucket'),5),
                                -p['research_followup_priority'],
                                {'A':0,'B':1,'C':2,'D':3}.get(p['research_evidence_tier'],4),
                                p['candidate'] or ''))
    with (out/'PATIENT_FIRST_ALL_ROUTE_PREDICTIONS.jsonl').open('w',encoding='utf-8') as f:
        for row in prediction_rows:
            f.write(json.dumps(row, ensure_ascii=False,sort_keys=True)+'\n')
    with (out/'PATIENT_FIRST_CANDIDATES.jsonl').open('w',encoding='utf-8') as f:
        for row in all_best:
            f.write(json.dumps(row,ensure_ascii=False,sort_keys=True)+'\n')
    top = all_best[:top_n]
    (out/'PATIENT_FIRST_TOP.json').write_text(json.dumps(top,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    fields=['candidate','scientific_group_id','research_evidence_tier','prediction_class',
            'mechanistic_priority_score','research_followup_priority','primary_gene_hypothesis','route_target','route_id',
            'registry_names','jurisdictions','pk_context_status','safety_review_status',
            'predicted_probability_of_clinical_benefit','limits']
    with (out/'PATIENT_FIRST_RANKING.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for row in all_best:
            w.writerow({'candidate':row['candidate'], 'scientific_group_id':row['scientific_group_id'],
             'research_evidence_tier':row['research_evidence_tier'],
             'prediction_class':row['prediction_class'],
             'mechanistic_priority_score':row['mechanistic_priority_score'],
             'research_followup_priority':row['research_followup_priority'],
             'primary_gene_hypothesis':row['primary_gene_hypothesis'],
             'route_target':row['route_target'], 'route_id':row['route_id'],
             'registry_names':'; '.join(row['registry']['regulator_names']),
             'jurisdictions':'; '.join(row['registry']['jurisdictions']),
             'pk_context_status':row['exposure_target_assay_context']['state'],
             'safety_review_status':row['public_label_safety_state'],
             'predicted_probability_of_clinical_benefit':'UNKNOWN',
             'limits':'; '.join(row['potential_risks_and_gaps'])})
    classes = dict(Counter(p['prediction_class'] for p in all_best))
    grades = dict(Counter(p['research_evidence_tier'] for p in all_best))
    summary = {
      'schema':VERSION, 'approved_single_ingredient_drug_groups_considered':len(all_best),
      'route_drug_pairs_analyzed':len(prediction_rows),
      'gene_hypotheses_in_case_model':len(_variant_profiles(patient_model)),
      'patient_genetic_variants_are_functional_hypotheses_not_proven':True,
      'mechanistic_prediction_status_counts':classes,
      'research_evidence_tiers_in_approved_prediction_groups':grades,
      'clinical_suitability_verified':False,
      'clinical_efficacy_probability_computed':False,
      'uncertainty_model_calibrated_against_independent_disease_data':False,
      'curated_pk_entries_loaded':sum(map(len,pk.values())),
      'approval_does_not_imply_indication_for_genetic_disease':True,
      'top_candidate_names':[p['candidate'] for p in top[:15]],
      'top_report_file':'PATIENT_FIRST_REPORT.md',
      'machine_results':'PATIENT_FIRST_CANDIDATES.jsonl',
      'patient_sequence_included':False,
    }
    (out/'PREDICTOR_SUMMARY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    lines=[
      '# Rare Disease MVA Pharma — Patient-first mechanism predictions (V2.8.8)',
      'Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch',
      '',
      '**Research predictions only. No drug on this list is thereby clinically appropriate, safe, or effective.**',
      'The previous numeric screening score is NOT a validated efficacy estimate.',
      'Investigation buckets T1–T4/HOLD, grounded in directness and source reliability, take precedence over weak numerical scores.',
      '',
      '## Interpretation',
      '- **Scientific evidence A/B/C/D** measures available published/assay evidence.',
      '- **Patient-first prediction** evaluates hypothetical direction of effect relative to this VCF-derived mechanism, where possible.',
      '- **Regulatory eligibility** requires a listed single-ingredient product in a pinned official source; product status and local availability still need review.',
      '- **Exposure unknown** means human/pediatric therapeutic target engagement cannot be estimated.',
      '',
      f"Approved research groups with a prediction record: {summary['approved_single_ingredient_drug_groups_considered']}.",
      f"Drug–route hypotheses analyzed: {summary['route_drug_pairs_analyzed']}.",
      '',
      '## Ranked candidate mechanisms and falsification requirements',
    ]
    for idx,r in enumerate(top,1):
        reg=', '.join(r['registry']['regulator_names']) or 'none'
        evid=r['research_evidence_tier']
        lines.extend(['',f"### {idx}. {r['candidate']} — {r.get('expert_investigation_bucket','T4')} — evidence {evid}",
          f"**Uncalibrated old heuristic (audit only):** {r['mechanistic_priority_score']:.1f}/100 — not used to establish probable efficacy",
          f"**Hypothetical gene/mechanism:** {r['primary_gene_hypothesis'] or 'not attributable'} → {r['route_target'] or 'unknown'} ({r['desired_action'] or 'unknown'})",
          f"**VCF-informed hypothesis:** {', '.join(r['affected_variant_effect_classes']) or 'unknown'}; phase: {r['phase_assessment']}",
          f"**Prediction:** {r['prediction_explanation']}",
          f"**Why evidence level {evid}:** {r['why_this_evidence_tier']}",
          f"**Approved research-product forms:** {', '.join(r['member_labels']) or 'unknown'}",
          f"**Approval-source registries:** {reg}. No disease-specific indication is inferred.",
          f"**Exposure/assay:** {r['exposure_target_assay_context']['state']}",
          f"**Safety review:** {r['public_label_safety_state']} (not patient clearance).",
          '**Main limitations:**'])
        for item in r['potential_risks_and_gaps'][:8]:
            lines.append('- '+item.replace('_',' '))
        if r['scientific_evidence_blockers']:
            lines.append('**Why this evidence grade cannot yet be upgraded:**')
            lines.extend('- '+x.replace('_',' ') for x in r['scientific_evidence_blockers'][:7])
        lines.append('**Next falsifiable investigations:**')
        lines.extend('- '+x for x in r['next_discriminating_experiments'][:5])
    lines += ['', '## Methodology and validity',
      'This model uses variant-effect *classes*, allele phase, available route directions, pharmacological target evidence, and official drug-register listings.',
      'It does not simulate a whole human physiology, an allele-specific protein structure, patient pharmacokinetics, or disease-specific clinical outcome.',
      'Its 0–100 score is an explicitly chosen NONCLINICAL screening heuristic and is not used to claim benefit.',
      'The expert triage distinguishes direct protein-target mechanisms from indirect RNA signatures.',
      'Only an independently source-verified exposure table allows an illustrative assay-versus-plasma comparison, never an efficacy prediction.',
      'When a finding is exploratory, negative or missing data do not imply the drug is ineffective.',
      'Source-bound AI literature claims require independent expert appraisal and functional experiments.',
      'Do not disclose patient genomic reports or personally identifying records in publications or public repositories.'
    ]
    (out/'PATIENT_FIRST_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return summary, top
