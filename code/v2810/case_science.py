"""Patient-first genetic uncertainty and mechanism objective model, local only.

Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch

This does NOT diagnose a child or recommend a medicine. It identifies
which experimentally verifiable phenotypes are needed to select research
interventions. Patient-specific variant strings are never sent to AI.
"""
from __future__ import annotations
import hashlib,re
from collections import Counter

GATE_VERSION="CASE_SCIENCE_V287_1"

def clinvar_interpretation(cv):
    if not cv:
        return {"state":"not_available","weight":0,"review":"unknown"}
    sig=str(cv.get("CLNSIG") or "").lower().replace("_"," ").replace("-"," ").strip()
    review=str(cv.get("CLNREVSTAT") or "").lower().replace("_"," ").strip()
    if not sig:
        return {"state":"not_available","weight":0,"review":review or "unknown"}
    # ClinVar review status "no_conflicts" is CONCORDANCE (two-star review),
    # not a conflict. Substring searching for "conflict" rejected it erroneously.
    review_conflict=("conflicting" in review or
        ("conflicts" in review and "no conflicts" not in review))
    if "conflict" in sig or review_conflict:
        return {"state":"conflicting","weight":0,"review":review}
    parts=[re.sub(r"\s+"," ",v.strip()) for v in re.split(r"[|/,;]",sig) if v.strip()]
    types=set()
    for p in parts:
        if p in ("pathogenic","likely pathogenic"):
            types.add("pathogenic")
        elif p in ("benign","likely benign"):
            types.add("benign")
        elif "uncertain" in p:
            types.add("uncertain")
        else:
            types.add("unknown")
    if len(types)!=1:
        return {"state":"conflicting_or_mixed","weight":0,"review":review}
    t=next(iter(types))
    if t in ("unknown","uncertain"):
        return {"state":t,"weight":0,"review":review}
    strong=("reviewed by expert panel" in review or "practice guideline" in review)
    moderate="multiple submitters" in review and "no conflicts" in review
    single=("criteria provided" in review and "no assertion" not in review
            and "no criteria" not in review and "not provided" not in review)
    if t=="pathogenic":
        w=10 if strong else 6 if moderate else 3 if single else 0
        return {"state":"reported_pathogenic","weight":w,"review":review or "not_supplied"}
    w=-10 if strong or moderate else -5 if single else 0
    return {"state":"reported_benign","weight":w,"review":review or "not_supplied"}

def presence_from_row(cells):
    """Only explicit negation counts as absent, and ambiguous status stays unknown."""
    data=[str(x or "").strip() for x in cells]
    first=data[0] if data else ""
    notes=" ".join((data[3:] if len(data)>3 else [])).lower()
    lower=first.lower()
    negative_patterns=[
      r"\babsent\b",r"\bnot present\b",r"\bdenies?\b",
      r"\breported absent\b",r"\bnegative for\b",r"\bno evidence of\b",
      r"^\s*no\s+\w",r"\bwithout\s+\w",
    ]
    if any(re.search(p,lower) for p in negative_patterns) or re.search(
        r"\b(?:explicitly denied|absent|not present|negative for|not observed)\b",notes):
        return "reported_absent"
    if re.search(r"\b(?:unknown|not assessed|uncertain|not evaluated)\b",notes):
        return "unknown"
    return "reported_present"

def phase_assessment(variants):
    """Only phase calls in same PS block are comparable; no parental proof implied."""
    if not variants:
        return {"state":"not_assessable","evidence":"no eligible variants",
                "trans_alleles_established":False}
    homo=[v for v in variants if str(v.get("GT") or "") in ("1/1","1|1")]
    if homo:
        # A homozygous call does not phase additional heterozygous variants.
        if len(homo)<len(variants):
            return {"state":"homozygous_with_other_unphased_variants",
                    "evidence":"Homozygous call present; phase of remaining variants unknown",
                    "trans_alleles_established":False}
        return {"state":"homozygous_observed","evidence":"VCF homozygous genotype",
                "trans_alleles_established":False}
    hetero=[]
    for v in variants:
        gt=str(v.get("GT") or "")
        if gt in ("0|1","1|0"):
            tags=v.get("phase_tags") or {}
            ps=tags.get("PS") or tags.get("PID")
            if ps: hetero.append((str(v.get("chrom") or ""),str(ps),gt))
    if len(variants)==2 and len(hetero)==2 and hetero[0][:2]==hetero[1][:2]:
        if hetero[0][2]==hetero[1][2]:
            return {"state":"cis_supported","evidence":"same haplotype in the same VCF phase block",
                    "trans_alleles_established":False}
        return {"state":"trans_supported_by_vcf","evidence":"opposite haplotypes in same VCF phase block; orthogonal confirmation not performed",
                "trans_alleles_established":True}
    return {"state":"unresolved","evidence":"insufficient comparable VCF phase tags / more than two alleles",
            "trans_alleles_established":False}

def _variant_group(v):
    c=str(v.get("consequence") or "").lower()
    h=str(v.get("hgvsp") or "").lower()
    if "frameshift" in c or "fs" in h:return "frameshift"
    if "stop_gained" in c or "nonsense" in c or re.search(r"(?:ter|\*)",h):return "premature_stop"
    if "splice" in c:return "splice_altering"
    if "missense" in c:return "missense"
    if "inframe" in c:return "inframe_indel"
    return "uncertain_effect"

def _safe_key(v):
    info="|".join(str(v.get(k) or "") for k in ("gene","chrom","pos","ref","alt","hgvsc","hgvsp"))
    return hashlib.sha256(info.encode("utf-8")).hexdigest()[:14]

def build_case_science(private):
    """Local genotype-based hypotheses. Only aggregate rationale is exportable."""
    by_gene=private.get("genes") or {}
    reports=[]
    for gene,variants in sorted(by_gene.items()):
        phase=phase_assessment(variants)
        effects=[]
        for v in variants:
            varianttype=_variant_group(v)
            cv=clinvar_interpretation(v.get("clinvar"))
            hypotheses=[]
            if varianttype in ("premature_stop","frameshift"):
                hypotheses=["possible NMD/transcript loss","possible truncated or partially functional protein"]
            elif varianttype=="splice_altering":
                hypotheses=["altered splicing / altered isoforms","possible reduced protein abundance"]
            elif varianttype in ("missense","inframe_indel"):
                hypotheses=["possible altered folding/stability","possible change in protein-specific activity",
                    "possible localization or protein-interaction defect"]
            else:
                hypotheses=["undetermined consequence; independent variant interpretation required"]
            effects.append({
                "variant_ref":"SHA256:"+_safe_key({**v,"gene":gene}),
                "effect_class":varianttype,
                "clinical_assertion":cv["state"],
                "clinvar_review":cv["review"],
                "clinvar_heuristic_weight":cv["weight"],
                "hypotheses_not_proven":hypotheses,
                "measured_functional_impact":False,
                "mutation_specific_therapy_direction_known":False,
                "requires": ["validated variant interpretation","variant-specific functional assay",
                    "allele-specific protein abundance/activity", "relevant cellular disease phenotype"],
            })
        reports.append({"gene":gene,"variant_count":len(variants),
          "phase":phase,"variant_effect_models":effects,
          "functional_objectives":[
              {"endpoint":"functional protein abundance","state":"measurement_required"},
              {"endpoint":"protein-specific biochemical and interaction activity","state":"measurement_required"},
              {"endpoint":"mitotic/checkpoint function or disease-specific cellular phenotype","state":"measurement_required"},
              {"endpoint":"rescue of disease-relevant functional defect under drug exposure","state":"measurement_required"},
          ],
          "actual_disease_mechanism_confirmed":False,
          "caveat":"Variant class alone does not prove causal mechanism or indicate a medicine."})
    present=[p for p in private.get("phenotype_features",[]) if p.get("presence")=="reported_present" and p.get("subject")=="proband"]
    absent=[p for p in private.get("phenotype_features",[]) if p.get("presence")=="reported_absent" and p.get("subject")=="proband"]
    uncertain=[p for p in private.get("phenotype_features",[]) if p.get("presence") not in ("reported_present","reported_absent")]
    return {"schema":GATE_VERSION,"primary_gene_hypotheses":reports,
       "phenotype_counts":{"present":len(present),"absent":len(absent),"unknown":len(uncertain)},
       "diagnostic_certainty":private.get("diagnostic_certainty","unknown"),
       "population_frequency_data":"not_loaded_unless_separate_local_source_provided",
       "clinical_variant_classification":"not_full_ACMG_AMP",
       "genomic_allele_strings_exported":False,
       "clinical_recommendation":False}

def case_suitability_constraints(private):
    """A conservative baseline relevant for every reported candidate."""
    domains=[]
    present=" ".join(str(p.get("term") or "") for p in private.get("phenotype_features",[])
      if p.get("presence")=="reported_present" and p.get("subject")=="proband").lower()
    for domain,words in {
       "genome_stability_or_cancer":("cancer","tumor","tumour","aneuploid","neoplasm"),
       "renal":("renal","kidney","nephro"),
       "growth":("growth","short stature","failure to thrive"),
       "hepatic":("liver","hepatic")}.items():
        if any(word in present for word in words):domains.append(domain)
    return {"assessment":"not_patient_specific_clearance",
      "risk_domains_to_review":domains,
      "required_missing_or_unverified":[
        "patient age and weight at proposed exposure",
        "medications and interactions",
        "renal and hepatic function", 
        "formulation, dose and pharmacokinetics",
        "paediatric indication and safety", 
        "long-term tumor and chromosomal instability effects",
        "potential worsening of the defective cell-cycle pathway",
      ],
      "drug_exposure_suitability":"unverified",
      "clinical_recommendation":False}
