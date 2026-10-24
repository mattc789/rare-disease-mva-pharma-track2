# Rare Disease MVA Pharma

## From patient-first discovery to testable pharmacology

### Bortezomib and trametinib: two evidence-led proposals within a conditional LZTR1–RIT1 checkpoint-modifier framework

**Scientific report and integrated methods disclosure**

Rare Disease, Real Kid — MVA Hackathon 2026  
**Track 2 · Drug Repurposing**

**Team:** ExaMedix Technologies  
**Research lead and developer:** Matteo Caprani  
**Software:** V3.1 Challenge Edition  
**Recorded scientific execution:** 22 October 2026  
**Report:** 23 October 2026 · Expanded mechanism-and-evidence edition

**The contribution**

A completed patient-first investigation connects an authorized genomic and phenotypic analysis to two specific pharmacological questions, attributable experimental literature, an inspectable software workflow and a staged validation program. The expanded report explains both the selection process and the biological experiments needed to establish relevance to the case.

**The proposals**

**Bortezomib:** investigate whether a recoverable LZTR1 deficit can be improved at the level of useful protein function, explicitly testing the competing effect of proteasome inhibition on RIT1 turnover.

**Trametinib:** investigate whether a measured LZTR1-associated signaling defect can be compensated through MEK inhibition, with functional and chromosome-segregation outcomes assessed separately from phospho-ERK.

**Research setting:** independent expert assessment and preclinical investigation. These are separate hypotheses, not a combination regimen or instructions to administer medicines to the participant.

Repository: https://github.com/mattc789/rare-disease-mva-pharma-track2

Original work: CC BY 4.0. Controlled case data and third-party resources retain their own access conditions and licenses.

<!-- PAGEBREAK -->

## Executive summary and methods abstract

**Purpose.** Rare Disease MVA Pharma proposes bortezomib and trametinib as two approved-medicine hypotheses for further investigation. The submission combines an executed computational workflow with a mechanism-specific experimental plan, addressing the Track 2 objective of identifying plausible research candidates rather than demonstrating an existing treatment. [1]

**Methods.** The V3.1 workflow starts with authorized genomic and phenotypic inputs, performs local annotation and case-focused population-frequency checks, and defines biological rescue objectives before searching medicines. Regulatory normalization, signed pathways, pharmacology, measured bioactivity, public perturbational-expression signatures, chemical similarity and source-bound literature review contribute distinct evidence channels. A separate program checks technical execution and provenance. Post-run review reconciles the recorded candidate dispositions and adds primary literature connecting LZTR1 regulation of RIT1 to mitotic-checkpoint function. [T1–T9]

**Results.** The completed run processed 7,536 medicinal entities, constructed 9,916 candidate paths and evaluated 7,266 route–drug pairs. Two final proposals were retained. Bortezomib has a published LZTR1-abundance observation in USP9X-perturbed melanoma-related cells; trametinib has a published cardiac-phenotype and phospho-ERK response in a specific dominant LZTR1-mutant mouse model. All 12 implemented technical-verification checks passed, with provenance for 107 distinct API responses; the project also records 280 successful synthetic tests. The retained list contains 87 BUB1B-associated exploratory entries, whose dispositions are reported without treating absent retrieval signals as negative biological experiments. [6–7; T1–T9]

**Mechanistic synthesis.** BUB1B remains the established disease-gene anchor. Published LZTR1–RIT1 regulation and RIT1–checkpoint experiments motivate a conditional modifier hypothesis: an actual LZTR1 defect might contribute to a BUB1B-associated checkpoint phenotype. This is a post-run research interpretation, not a demonstrated interaction in the participant. The decisive comparison is LZTR1 correction while holding the BUB1B background constant. Candidate testing then distinguishes protein abundance from substrate regulation, signaling correction from chromosome fidelity, and genuine functional improvement from cytostasis or selective cell loss. [5,39–40]

**Conclusion.** The submission delivers two traceable pharmacological propositions, a reusable research system and a practical sequence of causal and functional tests. The preserved strict outcome is `NO_VALIDATED_LEAD`; exposure qualification and safety review remain open. These boundaries identify the evidence needed for advancement while preserving the scientific value of the completed investigation. Neither medicine is proposed for administration on the basis of this work.

### At a glance

| Completed contribution | Next scientific decision |
|---|---|
| Executed case-to-candidate workflow with source-linked outputs | Establish whether the LZTR1 branch contributes to the measured case phenotype. |
| Two nominated medicines with distinct published experimental anchors | Test each separately in a qualified model at justified exposure. |
| Documented alternatives, uncertainty states and technical verification | Confirm useful functional improvement with independent replication. |

<!-- PAGEBREAK -->

## Reader’s guide and evaluation map

### How the argument is organized

| Part | Sections | Pages |
|---|---|---:|
| Disease context and the executed discovery method | 1–4 | 4–7 |
| Recorded results and BUB1B candidate dispositions | 5–7 | 8–10 |
| LZTR1–RIT1 convergence and case-specific attribution | 8–9 | 11–12 |
| Bortezomib: rationale and decisive experiment | 10–11 | 13–14 |
| Trametinib: rationale and decisive experiment | 12–13 | 15–16 |
| Shared validation, causal comparisons and endpoints | 14–17 | 17–20 |
| Reproduction, scalability and scientific contribution | 18–20 | 21–23 |
| Methods, governance, rights and provenance | Appendices A–D | 24–27 |
| Numbered references | References | 28–29 |

### The four official judging criteria

The organizers assign **Scientific Rigor 35%, Potential Impact 25%, Innovation 25% and Scalability 15%**. The report maps each criterion to concrete, inspectable deliverables; the judging panel assigns the scores. [1]

| Criterion | Evidence offered for assessment |
|---|---|
| **Scientific Rigor** | Variant-mechanism distinctions; documented selection decisions; primary-study context; competing mechanisms; genetic attribution and functional validation design. |
| **Potential Impact** | A specific way to determine whether an additional checkpoint-related axis matters in the case, and whether either nominated intervention warrants further laboratory investment. |
| **Innovation** | An implemented rescue-objective-first integration of computational discovery, source-level criticism and technical verification, extended by an explicit modifier-testing hypothesis. |
| **Scalability** | Reusable adapters, structured evidence objects, local genomic processing, reference provenance and a reviewer-operated installation route, with multi-case evaluation defined as the next test. |

### Four evidence labels used throughout

**Recorded result** identifies an output of the preserved software run or its supplied read-only documentation audit. **Published observation** identifies an experiment performed by the cited authors. **Research interpretation** identifies the reasoning connecting those observations. **Proposed experiment** identifies work still to be performed. This vocabulary makes the positive contribution and its current evidence level visible together.

References **[1–43]** identify external sources. **[T1–T9]** identify technical artifacts and the sanitized audit listed in Appendix D. Historical execution metrics come from those records; the current revision did not rerun the research. The new LZTR1–RIT1 synthesis and expanded laboratory design are dated post-run additions. Original candidate scores, evidence labels and hashes remain unchanged.

The public report uses gene-level findings and broad consequence classes. It excludes exact case alleles, genomic coordinates, individual identifiers, the source phenotype narrative, credentials and private API response identifiers.

<!-- PAGEBREAK -->

## 1. The biological problem and the research objective

### 1.1 Start with the disease, then define the useful intervention

MVA provides a demanding setting for drug repurposing because the relevant outcome is accurate chromosome segregation in viable tissue, rather than simply changing the abundance of a protein or suppressing cellular proliferation. Biallelic BUB1B variants are an established cause of constitutional aneuploidy and cancer predisposition. BUB1B encodes BUBR1, a component of the mitotic checkpoint complex, alongside MAD2, BUB3 and CDC20. [5,40–41]

This disease anchor determines the experimental objective: preserve or improve useful cell function while reducing a demonstrated disease-relevant defect. A compound that changes a signaling marker may be informative; a compound that changes chromosome fidelity would answer a more proximal question. These outcomes should be connected experimentally rather than treated as synonyms.

The proposal therefore distinguishes three possible achievements. **Molecular correction** improves a measured transcript, protein or signaling abnormality. **Cellular functional improvement** improves a defined output such as segregation fidelity or a model-qualified tissue phenotype. **Clinical benefit** requires a later body of evidence beyond this computational submission. The current work develops candidates and experiments at the first two levels.

### 1.2 Why variant mechanism matters

The same gene can be affected through loss of productive transcript, unstable protein, impaired interactions, altered localization or another functional change. BUBR1 studies show why these distinctions matter: selected MVA-associated missense proteins had low abundance and rapid turnover, while experimental normalization of some mutant constructs restored an engineered checkpoint response. That supports investigating recoverable function in selected alleles; it does not assign the same mechanism to every missense change. [41]

Consequently, the program uses the case analysis to generate provisional mechanisms and rescue objectives. It does not equate a missense label with instability, a splice annotation with a measured isoform, or a premature-stop annotation with a quantified absence of protein. This preserves multiple explanations until the appropriate evidence distinguishes them. [T3,T8]

### 1.3 Why consider a modifier alongside the primary disease anchor?

A patient-first investigation may identify an additional axis that could influence the same cellular phenotype or a separable component of it. Investigating such an axis is justified when there is a plausible mechanistic connection and a direct way to test its incremental contribution. It does not require replacing the established disease explanation.

That is the role proposed here for LZTR1. The expanded analysis asks whether altered LZTR1-dependent regulation could contribute to checkpoint dysfunction in a BUB1B-associated background. The literature supports the component links through RIT1; their operation together in this case is the research question. Sections 8–9 make the distinction explicit, and Sections 14–16 specify the comparisons that would resolve it. [39–40]

**The intended result is a qualified research decision:** identify a measurable, pharmacologically addressable contribution, or refocus the investigation when the proposed contribution is not supported. Both outcomes make the next experiment more informative.

<!-- PAGEBREAK -->

## 2. From authorized inputs to rescue objectives

### 2.1 Case intake and genomic evidence

The recorded investigation used the organizer-authorized single-proband GRCh38 compressed VCF and its matching structured clinical DOCX. Genomic annotation and phenotype-informed prioritization were performed locally. The workflow then carried up to three leading gene hypotheses into downstream investigation. This defined search boundary focuses downstream investigation while leaving alternative diagnoses available for separate assessment. [T1,T3,T8]

| Prioritized axis | Consequence classes represented in the analysis | Mechanistic question retained |
|---|---|---|
| **BUB1B** | Missense and premature-stop classes | What are the actual transcript, residual-product, checkpoint and segregation effects? |
| **LZTR1** | Premature-stop and splice-altering classes | Is productive transcript or functional protein altered, and does that change substrate regulation or downstream function? |
| **POU6F2** | Frameshift and missense classes | Does the proposed transcription-factor defect contribute to a measurable phenotype with an investigable intervention point? |

These broad annotations organize the hypotheses; the number of functionally distinct alleles and their phase remain to be established. Actual RNA products, allelic relationships and functional consequences remain model-qualification questions. The public account deliberately avoids exposing unrelated case-level findings.

The run passed checks for actual offline Ensembl VEP work and six required gnomAD chromosome/cohort resources covering the selected queries. These checks establish that specified computational work occurred. They do not substitute for orthogonal variant confirmation, direct phasing or new read-level analysis. No successful Pangolin prediction was generated; the splice-related branch therefore particularly benefits from direct transcript work. [T3,T7]

### 2.2 Translate each mechanism into an intervention objective

The software defines the desired biological change before searching for medicine names. Examples include increasing functional gene product, restoring full-length translation, correcting productive splicing, stabilizing a recoverable protein, restoring an interaction or compensating for a downstream deficit. Each objective carries requirements that must be satisfied before a hypothesis can be advanced. [T8]

For a possible protein-stability defect, the useful objective is **more functional residual protein**, not indiscriminate accumulation. For a downstream signaling defect, the useful objective is **normalization of a causally relevant output**, not inhibition merely because an inhibitor exists. For an unresolved mechanism, the immediate objective is measurement. This ordering connects the drug search to a falsifiable biological need.

### 2.3 Preserve uncertainty without stopping discovery

The workflow can explore a conditional route while recording what would invalidate it. This is valuable in a rare-disease setting where measured patient-cell data may not yet exist. It separates the permission to investigate an idea from the evidence needed to act on it.

The resulting patient-state and objective records supply the inputs for pharmacological discovery. They organize the authorized data into research hypotheses for subsequent clinical and functional qualification. Phase, causal contribution, expressed protein and intervention direction must be established where they determine the candidate’s mechanism. The staged experimental program follows precisely those dependencies.

<!-- PAGEBREAK -->

## 3. Complementary discovery channels

### 3.1 A broad search with identifiable evidence types

The discovery system combines several channels because useful pharmacology may appear as a drug–target interaction, an assay measurement, an expression signature, a chemical relationship or a statement in primary literature. Each channel contributes a different kind of support; agreement between channels is not automatically independent biological replication. [T4,T8]

| Channel | Implemented operation | Contribution to the research question |
|---|---|---|
| **Regulatory normalization** | Product records are organized into medicinal entities with eligibility attributes. | Defines the approved-medicine search space while distinguishing products from molecules. |
| **Pharmacology and signed pathways** | ChEMBL actions are combined with directional SIGNOR relationships and bounded propagation. | Proposes a route whose intervention sign can be assessed against a rescue objective. |
| **Measured bioactivity** | Existing assay values, endpoints, concentration and direction are retrieved and qualified. | Adds experimental pharmacology context from prior work. |
| **Perturbational expression** | Public LINCS/GEO signatures accessed through Enrichr are matched locally to intended up/down changes. | Adds a discovery signal from public perturbation data, not patient RNA measurements. |
| **Chemical similarity** | RDKit Morgan fingerprints, radius 2 and 2,048 bits, with a Tanimoto threshold of 0.60 against qualified references. | Extends discovery using explicit chemical similarity, not a predicted binding affinity. |
| **Direct and cross-source literature** | Source-linked drug, gene, intervention and effect statements are extracted and adjudicated. | Identifies the experiment supporting a hypothesis and the links still requiring validation. |

### 3.2 Why literature is a separate channel

An interaction database can establish that a drug acts on a target without establishing whether the action helps the disease mechanism. Conversely, a relevant experimental observation may not appear as a conventional drug–target edge. The direct-literature channel addresses that gap by seeking an attributable intervention on a prioritized gene product or a mutant-model phenotype. [T6,T8]

This mattered for bortezomib: its nomination was supported by a published LZTR1-abundance observation recovered through the literature channel rather than by a positive initial chemical/network signal. The final interpretation preserves that provenance. It does not portray every nominated medicine as having been independently rediscovered by all channels.

### 3.3 Capabilities used, and their scientific scope

The implemented computation includes local annotation, source queries, direction matching, fingerprints, structured reasoning and evidence checking. AlphaFold context, where retrieved, is structural context with model-confidence limitations. The executed default result did not include systematic docking, molecular dynamics, a patient-specific drug-bound structure, a whole-cell simulation or a qualified quantitative pharmacokinetic model. [T7–T8]

The delivered capability is a source-linked, mechanism-organized search and prioritization process. Multiple channels provide complementary opportunities to identify testable hypotheses. Source appraisal and the appropriate functional experiment then determine which routes warrant further investigation.

<!-- PAGEBREAK -->

## 4. From discovered paths to final proposals

### 4.1 Three decisions, not one score

The workflow separates **discovery**, **final-proposal eligibility** and the stricter **validated-lead** decision. A candidate can be worth scientific follow-up without already meeting a demonstrated-rescue standard. Equally, a high discovery score does not establish clinical usefulness. [T4–T8]

The heuristic ranking considers evidence tier, objective match, direct target assays, network distance, source-bound observations, coverage of supported rescue objectives and evidence blockers. It is an ordinal research-priority device. It was not calibrated as a probability of response or validated against a clinical outcome dataset.

A same-gene observation is particularly useful for specifying a molecular experiment. Its relevance to the participant still depends on the variant mechanism, model and tissue. Network closeness and biological transferability are therefore reported separately in the scientific interpretation: a short computational path does not erase a cross-disease inference.

### 4.2 Source-bound review and accountable curation

AI assists reasoning, extraction and source criticism under structured schemas and deterministic checks. The direct-rescue channel checks that the named medicine, gene and claimed effect remain tied to an available source passage; classification cannot legitimately add an experiment absent from that source. The workflow also records opposite effects, unresolved context and rejected extractions. [T6,T8]

The record contains an instructive example: a bortezomib-related extraction was rejected because the selected exact quotation did not contain a recognized candidate name or alias. This shows a concrete verification control, while also identifying an alias/passage-retrieval limitation. Repeated appearances of one paper are counted as one publication rather than multiple replications.

Automated discovery and ranking were followed by AI-assisted full-text checking, numerical reconciliation and author-led interpretation. For trametinib, full-text inspection recovered ventricular phospho-ERK evidence that the automated passage set had omitted. The report incorporates the observation while retaining the raw classification and score. [7; T6]

### 4.3 Why these two medicines emerged

**Bortezomib** met the literature channel’s criterion for a drug-associated change in a prioritized gene product. **Trametinib** supplied an intervention-associated phenotype response in a model carrying a mutation in that gene. Those are different experimental anchors and justify different follow-up questions. [6–7]

The implemented criteria prioritized these two types of same-gene intervention evidence. Case relevance remains the first laboratory question, and superiority over other MVA hypotheses has not been established. The recorded BUB1B dispositions and the added LZTR1–RIT1 connection together motivate a focused attribution experiment. This biological synthesis is explicitly a post-run interpretation.

The adjusted priorities remain **29 for bortezomib and 27 for trametinib**, with trametinib’s earlier value recorded as 32. The two-point difference is not an efficacy estimate and does not determine laboratory order ahead of model qualification. [T5]

<!-- PAGEBREAK -->

## 5. Recorded results and screening coverage

### 5.1 The completed computational funnel

| Recorded measure | Result | What the number represents |
|---|---:|---|
| Regulatory product records | 233,464 | Product records represented in the normalized catalogue. |
| Medicinal entities processed | 7,536 | Entries in the initial screening ledger. |
| At least one initial signal | 1,212 | Entries with a qualifying discovery signal. |
| No qualifying signal retrieved | 2,322 | Entries labeled `NO_SIGNAL` by the attempted channels. |
| Insufficient data | 3,855 | Entries without enough usable evidence for a supported screen result. |
| Incompatible modality | 147 | Entries outside the applicable screening modality. |
| Candidate paths constructed | 9,916 | Computational drug-to-objective paths. |
| Final biological routes | 104 | Routes after downstream expansion. |
| Route–drug pairs evaluated | 7,266 | Expanded assessment units. |
| Retained single-agent hypothesis rows | 100 | Bounded final exploratory list. |
| Final nominated medicines | 2 | Bortezomib and trametinib. |
| Qualifying combinations / validated leads | 0 / 0 | Separate strict decision outcomes. |

The four initial categories sum to 7,536. Initial screening used 98 routes, excluding four held routes; downstream expansion produced 104 final routes. These denominators describe different stages and should not be combined as independent candidate counts. [T4]

### 5.2 Signals, known pharmacology and source review

The run recorded **14,745 mechanism/network**, **1,013 measured-activity**, **4,042 perturbational-expression** and **52 chemical-similarity** signals. A medicine or route may contribute to more than one category. Known pharmacology was profiled for 6,240 entities: 975 had recorded mechanisms, 190 had multiple known targets and 22 had compatible route-mechanism matches. [T4]

Twenty final source-critique jobs completed. Among the retained 100 entries, 15 had individual final-review annotations: three conflicted, one partially supported, nine without source-bound support and two with counterevidence only. The other 85 remained `NOT_ATTACKED` in that artifact. Jobs and annotated entries are different units; the twenty jobs do not represent complete final review of one hundred medicines. Adagrasib and sotorasib were held for source conflicts. [T5–T6]

### 5.3 What coverage contributes to interpretation

The catalogue establishes search breadth, while the state labels show evidence depth. The **51.2% insufficient-data fraction** identifies a substantial retrieval/knowledge boundary; it is not a set of negative efficacy results. This distinction is part of the scientific output because it tells a reviewer where better data could change the answer.

The outcome is a focused pair of research propositions and an auditable record of how that pair was reached. The strict label `NO_VALIDATED_LEAD` remains consistent with Track 2’s hypothesis-generating purpose: the investigation produced proposals, while the evidence required for validated functional rescue remains to be generated. [1; T4–T6]

<!-- PAGEBREAK -->

## 6. Why BUB1B strategies were not final nominations

### 6.1 What the preserved hypothesis file actually contains

The supplied read-only audit of `TRACK2_HYPOTHESES_V31_FINAL.json` records **87 BUB1B-associated entries** among the retained 100. All 87 were `EXPLORATORY_ONLY`: four had `P0_NO_SOURCE_BOUND_SUPPORT`, one had `P0_COUNTEREVIDENCE_ONLY`, and 82 had `NOT_ATTACKED`. The remaining rows were twelve LZTR1-associated entries and one POU6F2 entry. Four LZTR1 rows had a plausible-hypothesis label, but only two medicines appear in the final submission-level list. [T9]

These observations show that BUB1B was extensively represented in the computational search. They do not establish that every potentially useful BUB1B medicine received a complete individual literature review. The retained list is a subset of the wider search, and final-critique status is distinct from earlier screening.

### 6.2 Selected dispositions and their reasons

| Strategy | Recorded evidence | Interpretation of non-promotion |
|---|---|---|
| **Arimoclomol citrate; sodium and glycerol phenylbutyrate** | `NO_SIGNAL`; zero mechanism/network and measured-activity records for these entries; transcriptomics not attempted; a structural record existed. | The attempted channels did not produce a qualifying signal. This is a retrieval outcome, not demonstrated lack of biological activity. |
| **Sirolimus** | BUB1B row at tier D, exploratory, without final individual critique. | The retained path did not supply the verified intervention-to-function evidence needed for promotion. |
| **Everolimus** | Selected BUB1B graph paths through CSNK2A1/ODC1 remained tier D with direction/function blockers. | These indirect paths did not qualify; this does not evaluate every possible mTOR-based rationale. |
| **Neratinib, erlotinib, gefitinib** | ERBB2-associated, multi-step tier-D paths; reviewed without sufficient source-bound support. | A target interaction did not establish beneficial BUB1B-relevant function along the proposed route. |
| **Cyclophosphamide** | Counterevidence-only state and a source-evidence submission hold. | The BUB1B-restoration proposition did not advance; standard oncological use is a separate question. |

The audit also records an isolated niacin perturbational-expression signal. That observation alone does not assign the signal to BUB1B or establish why niacin was absent from final nominations. It is preserved as a bounded ledger fact rather than expanded into an unsupported exclusion explanation. [T9]

**Selection conclusion.** BUB1B-associated strategies did not meet the implemented promotion criteria within the evidence retrieved and processed by this run. The two LZTR1 nominees carried specific same-gene experimental observations that supported an identifiable next test. This explains the software outcome without converting it into a biological ranking of the genes or a claim that all BUB1B alternatives were disproven.

<!-- PAGEBREAK -->

## 7. Scientific comparison with BUB1B-directed alternatives

### 7.1 An allele-dependent restoration strategy remains credible

BUBR1 restoration has a meaningful experimental foundation. Suijkerbuijk et al. found increased turnover of selected MVA-associated missense proteins; HSP90 inhibition depleted unstable constructs, whereas proteasome perturbation affected their accumulation. Expression-normalized I909T and L1012P constructs restored an engineered checkpoint readout. These results support a recoverable-abundance mechanism in selected experimental alleles. They do not establish drug responsiveness of the participant’s missense allele. [41]

An arimoclomol or phenylbutyrate strategy would therefore benefit from the same discipline applied to our nominees: establish the actual defect, show that useful function remains recoverable, then test the medicine rather than assuming that its broad proteostasis properties reproduce the published genetic/expression intervention. The run’s `NO_SIGNAL` states do not close this scientific question.

### 7.2 NAD–SIRT2 and mTOR are distinct hypotheses

North et al. linked SIRT2 to BUBR1 stabilization, studied NMN-mediated abundance changes and reported benefits of genetic SIRT2 elevation in hypomorphic mice. These interventions are not interchangeable with niacin or nicotinamide exposure, and the experiments do not measure rescue of this participant’s genotype. The study’s nicotinamide inhibition experiments further show why precursor identity and biological direction must be checked rather than inferred from a shared vitamin family. [42]

Sieben et al. reported allele- and tissue-dependent phenotypes and an mTORC1 association in selected BubR1-mutant muscle contexts. That supports measuring whether the pathway is active in a relevant model; it is not a direct demonstration that an mTOR inhibitor corrects this case’s checkpoint defect. [43]

### 7.3 Apply one evidentiary standard to every route

| Research route | Evidence needed to justify the first drug experiment |
|---|---|
| **BUBR1 stabilization** | An abundance/stability deficit with residual function that improves when protein quantity is normalized. |
| **NAD/SIRT2 modulation** | A responsive, directionally appropriate stabilization mechanism for the specific tested precursor and allele context. |
| **mTOR compensation** | A measured, causally relevant pathway abnormality and a functional endpoint appropriate to that downstream hypothesis. |
| **LZTR1 abundance intervention** | A relevant LZTR1 contribution, recoverable product and preserved substrate-regulating function under intervention. |
| **LZTR1 downstream compensation** | A measured signaling abnormality whose correction improves a relevant function, with checkpoint effects assessed independently. |

All five routes share the same evidentiary standard. The LZTR1 nominees enter with specific published observations; case-specific qualification determines whether those observations support further experimentation.

### 7.4 Why the LZTR1 exploration adds value

The choice is not “BUB1B or LZTR1.” The primary disease anchor and the modifier hypothesis can be tested together while maintaining distinct causal claims. The LZTR1 branch is useful because the newly appraised RIT1 literature provides a specific point of convergence on checkpoint function and suggests a discriminating experiment. It widens the research question without treating the original BUB1B mechanism as resolved or irrelevant. [39–40]

Direct BUB1B strategies remain appropriate comparators in a future systematic evaluation. Their broader comparative assessment is prospective; the final nominations and the historical selection record remain unchanged.
<!-- PAGEBREAK -->

## 8. The experimental connection: LZTR1, RIT1 and the checkpoint

### 8.1 LZTR1 regulates a protein with direct mitotic effects

Castel et al. demonstrated CUL3/LZTR1-dependent proteasomal regulation of RIT1. LZTR1 loss increased RIT1 abundance in experimental cells, and fibroblasts from children with biallelic LZTR1-associated Noonan syndrome showed higher RIT1 than parental and unrelated controls. This provides a human-cell precedent for the regulatory relationship, although it concerns a different disorder and different variants. [39]

Cuevas-Navarro et al. then demonstrated direct RIT1 interactions with MAD2 and p31comet. In their systems, LZTR1 knockout or stabilized pathogenic RIT1 accelerated mitotic progression; pathogenic RIT1 promoted chromosome errors and aneuploidy. RIT1 also reduced MAD2–CDC20 and BUBR1–CDC20 association, connecting its action to the mitotic checkpoint complex rather than only to a generic growth pathway. [40]

These complementary studies provide the biological reason to examine LZTR1 alongside BUB1B. Their authorship overlaps; independent replication of a combined disease mechanism remains a separate validation step.

### 8.2 The point of convergence

| Evidence layer | Supported relationship | Contribution to this proposal |
|---|---|---|
| Established MVA biology | BUB1B dysfunction can compromise BUBR1-dependent checkpoint and segregation functions. [5,41] | Defines the disease-relevant system. |
| Published LZTR1 biology | LZTR1 helps regulate RIT1 turnover. [39] | Identifies a measurable upstream regulatory variable. |
| Published RIT1 biology | Abnormal RIT1 levels can perturb the MAD2/p31comet–checkpoint system. [40] | Connects the modifier hypothesis to mitotic fidelity. |
| Present research interpretation | A functional LZTR1 alteration might modify a BUB1B-associated cellular phenotype. | Defines a case-specific experiment; it is not yet an observed interaction. |

The proposed convergence is functional, through the shared checkpoint system. Direct LZTR1–BUBR1 binding and recognition of LZTR1 as an additional MVA disease gene are outside the evidence established here.

### 8.3 Normalization, not indiscriminate suppression

The RIT1 study also found mitotic delay and increased segregation errors after RIT1 depletion. Consequently, “less RIT1” is not a sufficient treatment objective. The proposed experiment asks whether an abnormal regulatory state can be normalized while maintaining productive, accurate division. [40]

Likewise, RIT1’s checkpoint-related behavior in the tested constructs was distinguishable from a simple basal MAPK readout. This does not prove that every mitotic effect is MEK-independent. It establishes why phospho-ERK and checkpoint function must be measured separately in the trametinib experiment.

### 8.4 What this addition changes

This literature synthesis strengthens the **reason to investigate the axis** and specifies the necessary readouts. Original scores remain unchanged; drug activity in the participant’s relevant model is the next experimental question. The supplied audit did not find the 2021 paper’s identifiers in the mechanism/modifier literature bundles it searched; that was a bounded search, not an inventory of every file. We therefore identify the bridge transparently as a post-run appraisal. [T9]

<!-- PAGEBREAK -->

## 9. From a pathway connection to a case-specific hypothesis

### 9.1 The hypothesis we propose to test

**In a qualified BUB1B-associated model, an independently demonstrated LZTR1 functional alteration may contribute to checkpoint dysfunction through abnormal RIT1 regulation, or to another measured phenotype through altered downstream signaling.**

This formulation preserves BUB1B as the disease anchor and gives LZTR1 a defined, conditional role. It also separates two possibilities that require different experiments: an additional contribution to chromosome instability, and a parallel signaling-associated phenotype. Neither possibility is assumed from a stop or splice annotation alone. [5,39–40; T3]

### 9.2 The minimum causal chain

| Link to establish | Evidence currently available | Evidence that would resolve the next step |
|---|---|---|
| Case-level LZTR1 alteration affects a productive transcript or protein | Broad genomic consequence hypotheses in the case analysis. | Orthogonal confirmation, relevant phase, allele-specific RNA/splicing and full-length protein measurements. |
| The molecular alteration changes substrate regulation or signaling | Public LZTR1/RIT1 experiments in other contexts. | Quantified RIT1 turnover and pathway-state measurements in the qualified case model. |
| The alteration contributes to a relevant functional defect | A literature-supported convergence on the checkpoint. | LZTR1 correction with the BUB1B background held constant, plus appropriate procedure controls. |
| The medicine improves that defect at suitable exposure | Candidate-specific published observations elsewhere. | Controlled, exposure-qualified drug experiments with functional and safety endpoints. |

The chain makes the proposal actionable: each result determines which subsequent experiment is warranted, and each positive molecular observation is evaluated in the context that gives it biological meaning.

### 9.3 The decisive attribution comparison

The central contrast is **the same BUB1B background with altered versus corrected LZTR1**. A reproducible improvement after LZTR1 correction would support a modifier or parallel contribution in that model. A separate BUB1B-corrected arm would estimate the primary gene’s contribution. Correcting both would reveal whether their effects are additive, interactive or largely attributable to one axis.

A modifier contribution does not require statistical synergy. It requires a defensible change attributable to LZTR1 under a controlled genetic background. Conversely, a change in a generic proliferation marker would not establish improvement of the chromosome-segregation phenotype.

### 9.4 Alternative explanations remain testable

An LZTR1 annotation may not alter the relevant transcript; a measurable molecular effect may not contribute to the functional phenotype; or the phenotype may be explained predominantly by BUB1B. These outcomes would redirect resources toward the better-supported route. They do not invalidate the executed search, because resolving competing mechanisms is the purpose of the proposed causal comparison.

Model qualification therefore precedes candidate-specific interpretation. Failure to recover a technically valid corrected line is not biological evidence against LZTR1, and failure of a drug in an unqualified model is not a decisive test of the hypothesis. Sections 14–17 specify the controls that keep these decisions distinct.

<!-- PAGEBREAK -->

## 10. Bortezomib: a gene-product intervention question

### 10.1 Why this medicine emerged

Bortezomib was retained because the discovery process found a published pharmacological observation involving **LZTR1 itself**, rather than only a nearby target or a similar chemical structure. The intervention increased detectable LZTR1 under a defined perturbation. That gives the proposal a concrete molecular starting point and a direct measurement for the first experiment. [6; T5–T6]

The scientific question is whether a relevant, recoverable LZTR1 deficit can be improved **at the level of useful function**. The nomination is conditional on demonstrating that such a deficit exists and that the model retains product capable of contributing to function.

### 10.2 The published observation and its context

Bacchiocchi et al. examined LZTR1 in melanoma-related systems. USP9X depletion reduced LZTR1 abundance, and Velcade/bortezomib restored detectable protein in the studied melanoma/melanocyte settings. The same paper reports context-dependent growth, invasion and survival roles for LZTR1 in melanoma. It therefore supports an abundance-modulation experiment while also motivating a careful proliferation and oncology assessment. [6]

Transfer to the participant’s LZTR1 context requires a specific qualification step. Pharmacologically preserving protein after increased turnover differs from recovering full-length product after transcript loss or abnormal splicing. The RNA and protein qualification stage establishes whether the published intervention logic is relevant to the actual model.

### 10.3 A traceable route from observation to experiment

| Step | Rationale | Required measurement |
|---|---|---|
| Establish a relevant deficit | A beneficial abundance intervention requires a deficit that matters to function. | RNA, isoforms, full-length LZTR1 and genetic attribution. |
| Establish recoverability | Preserving an inactive fragment would not recover the proposed activity. | Product identity, localization and substrate-regulating capacity. |
| Test the recorded pharmacological observation | Bortezomib supplies a same-gene experimental anchor. | LZTR1 abundance before and after exposure. |
| Determine the net biological effect | Proteasome inhibition changes many processes, including the route LZTR1 uses to regulate substrates. | RIT1, relevant signaling, checkpoint function and productive cell output. |

The proposal’s value is that these steps distinguish a useful functional effect from an attractive but insufficient protein-band result.

### 10.4 Preserved evidence state

The run’s final audit attributed the positive observation to **one unique primary publication**, with **zero accepted same-drug functional-rescue sources**. Eight documents and eighteen passages entered the recorded source-extraction/review context. The preserved state is **P1 partial source-bound support**, with adjusted priority **29**. Repeated appearances of the same paper do not constitute replication. [T5–T6]

The label `PRIMARY_GENE_PRODUCT_RESCUE` is retained in the original output. In the scientific interpretation it denotes **abundance support**. The next section explains the functional test needed to justify a stronger claim.

<!-- PAGEBREAK -->

## 11. Bortezomib: resolve the abundance–function trade-off

### 11.1 A mechanistically important competing effect

Castel et al. found that bortezomib prevented the reduction in RIT1 produced by LZTR1 in their experimental system. This is consistent with the proteasome being required for the substrate-regulating activity under study. Thus, increased detectable LZTR1 and impaired RIT1 degradation can coexist. [39]

This observation does not negate the abundance result reported by Bacchiocchi et al.; the studies measure different effects in different settings. It changes what would count as success. **The relevant outcome is the net functional state of the LZTR1–RIT1 system, not the LZTR1 band alone.**

The distinction is especially important when evaluating mitosis. An intervention affecting proteasomal disposal of mitotic regulators can alter exit from mitosis independently of checkpoint restoration. A longer arrest therefore needs to be interpreted alongside substrate turnover and segregation fidelity, not automatically scored as a stronger checkpoint. [39–40]

### 11.2 Proposed discriminating experiment

First establish a valid LZTR1-dependent defect through the genetic comparisons in Section 15. Then compare vehicle and bortezomib in affected and corrected cells under measured exposure conditions. Read LZTR1 quantity alongside protein identity, RIT1 abundance/turnover, relevant downstream signaling, checkpoint-complex behavior and chromosome outcomes.

A recovery or washout arm may help distinguish persistence of a useful protein effect from ongoing global proteasome inhibition. This is an experimental discriminator, not a proposed clinical dosing schedule. Any such arm needs measured recovery of proteasome function and appropriate time-matched controls before its results are interpreted.

### 11.3 Decisions defined before outcomes

| Observation | Scientific interpretation and next action |
|---|---|
| LZTR1 increases, substrate regulation normalizes and relevant function improves with preserved output | Supports further mechanistic confirmation and independent replication of this conditional drug hypothesis. |
| LZTR1 increases but RIT1 regulation or the checkpoint worsens | Supports the competing mechanism; do not advance the abundance-rescue proposition. |
| LZTR1 increases without functional improvement | Demonstrates pharmacodynamic abundance change, not useful rescue. |
| Product is absent, nonfunctional or unrelated to the model defect | Redirect the experiment; the prerequisite for this stabilization strategy is not met. |
| Apparent improvement is explained by toxicity, cytostasis or loss of abnormal cells | Classify the effect separately; it does not meet the functional-improvement criterion. |
| Exposure or assay validity is inadequate | Resolve the technical limitation before drawing a biological conclusion. |

### 11.4 Why retain a tightly qualified nomination

The medicine is retained as a **specific, experimentally anchored hypothesis with a decisive competing-mechanism test**. The evidence does not yet establish a favorable operating window. That is precisely the question the staged study would answer, and it should be answered before larger translational investment.

The original **exposure `NOT_ASSESSED`** and **safety `HOLD`** states remain unchanged, as does priority 29. The post-run RIT1 appraisal supplies the competing-mechanism test that the next experiment must resolve. [T5–T7]

<!-- PAGEBREAK -->

## 12. Trametinib: downstream compensation with a functional precedent

### 12.1 Why this medicine emerged

Trametinib offers a different research opportunity: compensation for a demonstrated downstream signaling consequence rather than preservation of residual LZTR1 protein. Its experimental anchor includes both a phenotypic response and a pathway readout in an LZTR1-mutant animal model. This makes it a substantive follow-up hypothesis while leaving the model-to-case transfer explicit. [7; T5–T6]

The proposed sequence is **qualify an LZTR1 contribution → establish the direction and relevance of signaling → test MEK inhibition → measure useful function**. An annotation alone does not establish that the relevant pathway is overactive or that inhibiting it will help.

### 12.2 The published intervention

Abe et al. studied a specific autosomal-dominant, Noonan-like LZTR1-mutant mouse model. Trametinib improved cardiac hypertrophy in the relevant mutant male group and reduced ventricular phospho-ERK. Splenic and renal enlargement did not improve. These findings establish a bounded, model-specific response rather than uniform correction of every phenotype. [7]

The dominant mutant mechanism and cardiac endpoint differ from the broad premature-stop/splice hypotheses considered in this case. The study supplies an intervention precedent whose transfer requires allele- and phenotype-specific qualification. It supports asking a specific question: whether a measured and functionally relevant signaling abnormality in a qualified case model is similarly modifiable.

### 12.3 What the candidate contributes

| Contribution | Significance for the proposed study |
|---|---|
| A named approved medicine with a defined pathway action | Provides an experimentally actionable perturbation rather than an unspecified pathway suggestion. |
| A same-gene mutant-model phenotype response | Grounds the nomination in an intervention result, not only an association. |
| A pathway-level pharmacodynamic observation | Supplies a measurement to confirm that the intended signaling intervention occurred. |
| Endpoint-specific response in the published model | Directs the new study to define its own relevant functional endpoint before treatment. |

The drug would be tested independently of bortezomib. Evidence for one branch cannot compensate for failure of the other, and no synergy or combination benefit is inferred.

### 12.4 Source reconciliation improves the scientific account

The automated passage set had omitted the ventricular phospho-ERK observation in Figure 6C. Post-run full-text review restored that observation to the report, while preserving the original receipts and classification. This is a concrete example of source review strengthening a candidate’s description without rewriting its computational history. [7; T6]

The final dossier retains **one positive primary mutant-model publication**, adjusted priority **27** and the original **P0 context-conflicted** state. The earlier priority was 32. These are historical research-priority values, not effect sizes or probabilities. Section 13 explains why the conflicting context informs the experiment without being misrepresented as a failed same-drug replication.

<!-- PAGEBREAK -->

## 13. Trametinib: separate signaling correction from checkpoint recovery

### 13.1 Interpret the counterevidence at the correct level

The source audit’s cautionary item was a secondary review describing unsuccessful **pimasertib** intervention in a different Lztr1-deficient vascular setting. It is not a direct replication of Abe et al.’s trametinib experiment. The appropriate inference is context dependence across compound, genotype, tissue and endpoint—not universal failure of MEK inhibition or of trametinib. [8; T6]

This distinction supports a more useful protocol: determine which signaling state is present in the model, and whether changing it improves the preselected functional defect. The raw conflict label remains available for reproducibility; the narrative supplies the biological interpretation it requires.

### 13.2 Two related questions, two separate readouts

The LZTR1–RIT1 literature connects protein regulation with checkpoint behavior, but does not establish that reducing phospho-ERK necessarily corrects that behavior. A pathway marker can respond while MAD2/CDC20/BUBR1 complex function remains abnormal. Conversely, a relevant non-mitotic function may improve without correcting the primary chromosome defect. [39–40]

The study should therefore distinguish **checkpoint rescue** from **downstream compensation**. Both can be research-relevant, but they support different claims. A compensatory benefit must be tied to a demonstrable case-relevant phenotype and assessed for adverse chromosome effects; it must not be renamed as correction of the initiating BUB1B defect.

### 13.3 Proposed test sequence

Before drug exposure, establish baseline RIT1/substrate regulation, pathway activation and the functional phenotype in affected and genetically informative control cells. Test whether LZTR1 correction changes those variables. This determines whether there is a causal signaling question for trametinib to address.

Then compare vehicle and trametinib across affected and corrected backgrounds. Measure pathway engagement, the predeclared functional endpoint, division completion, segregation errors and cell fate under the same exposure conditions. An orthogonal pathway perturbation can strengthen attribution, provided its own effects on viability and mitosis are measured rather than assumed.

| Result | Interpretation |
|---|---|
| Appropriate pathway correction plus functional improvement and preserved chromosome fidelity/output | Supports independent confirmation of downstream compensation; checkpoint rescue requires its own positive endpoints. |
| Lower phospho-ERK without useful functional change | Confirms pathway modulation but does not establish the proposed benefit. |
| Functional effect without evidence connecting it to the proposed signaling mechanism | Retain as a separate exploratory observation requiring attribution. |
| No relevant baseline abnormality or no LZTR1 contribution | The specific case-linked rationale does not qualify for advancement. |
| Adverse division, survival, chromosome or growth effects | The measured benefit–harm pattern does not meet the proposed advancement conditions. |

### 13.4 Present status

The nomination is supported by an identifiable intervention study and an explicit route to test transfer. **Exposure remains `NOT_ASSESSED` and safety remains `HOLD`** in the recorded output. The proposed experiments establish whether the candidate has a useful case-relevant role; they do not authorize clinical use. [T5–T7]
<!-- PAGEBREAK -->

## 14. Shared validation: establish the model before testing medicines

### 14.1 Scope and responsibility

The following program is a **proposed laboratory study**, developed from the recorded case analysis and the cited primary literature. It has not been executed or prospectively registered. The executing investigators would finalize sample access, assay qualifications, exposure conditions, analysis and stopping criteria before outcome unblinding. Existing authorized material should be prioritized; the report does not authorize new sampling or contact with the family.

The shared first stage is valuable because it qualifies both drug questions together. It establishes what the relevant alleles actually do and which functional readout can distinguish a meaningful improvement from a nonspecific cellular response.

### 14.2 Molecular and cellular qualification

| Qualification question | Proposed approach | Decision enabled |
|---|---|---|
| Are the candidate calls and their relevant relationships confirmed? | Orthogonal confirmation; phase and copy-number assessment where necessary; specimen identity controls. | Define the tested genetic configuration without assuming causality from two annotations. |
| What is expressed? | Allele-specific RNA, splice-isoform assessment, full-length versus truncated protein measurements. | Distinguish transcript loss, altered product and a potentially recoverable abundance deficit. |
| Which functions are impaired? | BUBR1 abundance/localization and checkpoint behavior; LZTR1 substrate regulation; RIT1 and relevant signaling; unperturbed division outcomes. | Select the actual functional defect and the correct direction for intervention. |
| Can the phenotype be attributed genetically? | Matched correction, reconstitution calibrated to endogenous expression, and procedure controls. | Qualify the proposed BUB1B and LZTR1 contributions. |

Assay selection is informed by the published BUBR1 and LZTR1/RIT1 systems, but must be qualified in the model being used. A source assay’s success in another cell background does not establish sensitivity or specificity here. [39–41]

### 14.3 Build a model that can support the intended inference

Patient-derived cells, where authorized and suitable, provide the relevant background. Engineered isogenic cells can separate variant effects and supply sufficient material, but should be described as mechanistic models rather than replicas of every tissue phenotype. A nontransformed background with documented checkpoint and p53 status is preferable for constitutional-cell questions; tumor models serve a separate compatibility analysis.

Use independently initiated editing or culture preparations, identity and contamination checks, sequence/phase verification and karyotype/copy-number assessment. Track passage and proliferation state. A corrected clone should be supported by more than one independently derived preparation and by procedure-matched controls; otherwise clone selection can mimic a causal effect.

### 14.4 A valid negative is different from an invalid experiment

A technically valid model with no relevant LZTR1-dependent defect would redirect the two proposed branches. An uninterpretable RNA assay, an unstable culture or an unsuccessful correction is instead a model-development result. This distinction preserves the value of the study: decisions follow the evidence actually obtained, rather than a forced positive or negative label.

<!-- PAGEBREAK -->

## 15. The decisive causal experiment: hold BUB1B constant

### 15.1 A genetic comparison with a specific answer

The strongest experiment for the modifier hypothesis compares **case-associated LZTR1 function with corrected LZTR1 function on the same BUB1B background**. It asks whether LZTR1 contributes to the measured phenotype, rather than whether the two genes appear in related pathways.

An expanded four-state design separates the primary and additional contributions. “Case-associated” below names the tested configuration; it does not prejudge pathogenicity of every allele.

| State | BUB1B configuration | LZTR1 configuration | Main comparison enabled |
|---|---|---|---|
| **A** | Case-associated | Case-associated | Baseline phenotype of the combined configuration. |
| **B** | Case-associated | Corrected/reference | LZTR1 contribution with BUB1B held constant: A versus B. |
| **C** | Corrected/reference | Case-associated | BUB1B contribution with LZTR1 held constant: A versus C. |
| **D** | Corrected/reference | Corrected/reference | Joint correction and reference for residual background effects. |

The design can be staged: the A–B contrast answers the immediate modifier question, while C–D refine attribution. Healthy and mock/procedure controls complement these comparisons rather than replacing the matched genetic contrast.

### 15.2 Interpret contribution and interaction separately

If correcting LZTR1 improves a relevant endpoint in the BUB1B-associated background, that supports a contribution in that model. The improvement need not eliminate all abnormalities, because BUB1B remains altered. Comparing the LZTR1 effect in both BUB1B backgrounds can then ask whether the effect depends on BUB1B state.

For an endpoint Y in which higher values indicate greater defect, define the LZTR1 correction contrasts as **Y(A) − Y(B)** and **Y(C) − Y(D)**. Their difference is the interaction contrast on that scale. It is not a measured result or an assumed synergy; uncertainty and scale dependence must be reported. A meaningful modifier effect can exist without a statistically established interaction.

### 15.3 Test RIT1 mediation without assuming complete removal is beneficial

After establishing a relevant LZTR1 contrast, measure whether correction normalizes RIT1 turnover and checkpoint-complex behavior. An independently controlled RIT1 perturbation can test whether changing that intermediate reproduces part of the functional effect. Because both excessive and depleted RIT1 disrupted mitotic behavior in published experiments, the intended comparison is restoration toward a qualified reference range, not indiscriminate depletion. [39–40]

Concordance among genetic correction, substrate regulation and functional outcome would strengthen attribution. Discordance would refine the model—for example, an LZTR1 effect may occur through another substrate or a non-mitotic function.

### 15.4 Link this experiment to the drug studies

The genetic comparisons establish the attainable direction and approximate scale of useful change. Drug effects can then be evaluated against those controls, without demanding complete normalization by assumption or accepting a statistically significant marker change as sufficient. This makes the proposed pharmacology both more selective and more informative.

<!-- PAGEBREAK -->

## 16. Functional endpoints: recognize genuine improvement

### 16.1 Measure what happens to the starting population

The primary outcome should be selected from the model-qualified defect before drug results are examined. For the checkpoint-modifier question, track a cohort of cells enrolled before exposure through their first attempted division and subsequent daughter-cell fate. Record chromosome errors, completion, death, arrest and productive daughter-cell output separately.

Let **N** be the number of enrolled founder cells, **C** the number completing the prespecified first division, and **E** the number of those divisions with at least one prespecified segregation error. Report **E/C**, **E/N** and **C/N**, alongside cell loss and viable daughter output. None of these ratios is sufficient on its own: fewer completed divisions can lower E/N without improving fidelity, while selective loss changes the population represented by E/C.

This proposed measurement structure is an analysis design, not a validated classifier already implemented in the submitted software.

### 16.2 Complementary endpoint families

| Endpoint family | Proposed readout | Role in interpretation |
|---|---|---|
| **Molecular mechanism** | Productive RNA, full-length protein, RIT1 turnover, relevant pathway state. | Establish that the proposed molecular intervention occurred. |
| **Checkpoint function** | Challenge-response kinetics, complex integrity and timing of mitotic progression. | Distinguish checkpoint behavior from a nonspecific arrest. |
| **Segregation fidelity** | Blinded first-division scoring of lagging chromosomes, bridges and other prespecified errors. | Test the proximal chromosome-level question. |
| **Cell fate and useful output** | Completion, death, arrest, daughter survival and subsequent division. | Identify loss of tissue-producing capacity or survivor selection. |
| **Longer-term composition** | Appropriately qualified karyotype or single-cell chromosome measurements. | Follow persistent change without substituting prevalence for error generation. |

A micronucleus endpoint needs a definition that distinguishes whole-chromosome segregation effects from other sources of damage. Baseline cell-cycle state, passage, labeling and imaging exposure should be matched and qualified. The assay itself must not create the phenotype it is intended to measure.

### 16.3 Keep attribution and benefit aligned

For bortezomib, a larger LZTR1 band must be accompanied by preserved or improved substrate regulation and relevant function. For trametinib, a lower phospho-ERK signal establishes pathway modulation, not checkpoint recovery. A distinct compensatory endpoint may be appropriate only when its relevance has been established beforehand. [6–7,39–40]

The prospective decision rule should require improvement in the primary functional endpoint **and** preservation of division/output within justified margins. It should also account for abnormal-clone persistence and toxicity. Margins must be fixed from assay performance and biological controls before outcome unblinding; this report does not invent a universally valid percentage.

### 16.4 Independent confirmation

A qualified positive should be repeated in independently prepared material and, where feasible, another appropriate background or laboratory. Agreement across endpoint families is valuable because each fails differently. Technical wells, individual cells and repeated readings do not become independent biological replicates merely by increasing their number.

<!-- PAGEBREAK -->

## 17. Exposure, safety and statistical decisions

### 17.1 Exposure is part of the experiment, not a later translation claim

Both nominees retain `NOT_ASSESSED` exposure status in the executed dossier. Before interpreting a functional response as translationally relevant, establish the tested drug identity, matrix, concentration, duration and recovery conditions. Measure free medium exposure and, where informative, intracellular exposure or a fit-for-purpose pharmacodynamic readout.

A published nominal cell concentration cannot be equated automatically with a total plasma peak. Binding, tissue distribution, intracellular access and time course can change that comparison. Clinical label information can inform the experimental boundary, but does not provide an MVA-specific effective concentration. No patient dose or culture-to-child conversion is proposed here. [T7]

### 17.2 Safety questions follow the mechanism

For bortezomib, the key mechanistic concern is that proteasome inhibition may preserve detectable LZTR1 while impairing substrate regulation or mitotic progression. For trametinib, pathway suppression may occur without restoring the relevant function or may reduce useful cell output. These concerns are tested in affected and corrected constitutional-cell comparators, not inferred solely from approval in other indications. [6–7,39–40]

The oncology context requires a separate compatibility assessment: an intervention should not be advanced on a favorable constitutional-cell marker while increasing the persistence or fitness of harmful clones. Any tumor-model or treatment-interaction study needs the appropriate tumor and regimen context. The report does not propose changes to current oncology care, and culture results cannot establish long-term pediatric safety.

### 17.3 Prespecify the unit of inference and the claim family

The executing team should use pilot information to estimate variability and feasible replication before fixing the confirmatory design. Independent preparations or editing events are higher-level units; cells and wells are nested observations. Record batch, clone, passage and imaging effects, randomize allocation and mask scoring.

Choose one primary functional claim per candidate before confirmation. Control multiplicity across those claims and any prespecified dose comparisons—for example, a justified Holm procedure for the two candidate-level claims—while labeling additional mechanisms and endpoints exploratory. This is a proposed analysis strategy, not a completed power analysis or a claim that a fixed sample size is already sufficient.

### 17.4 Four useful outcomes

| Outcome | Meaning for the research program |
|---|---|
| **Qualified positive** | Mechanism, functional improvement, exposure and preservation criteria support independent confirmation. |
| **Precise lack of useful effect** | A valid experiment excludes the prespecified meaningful improvement in the tested setting. |
| **Adverse functional effect** | A relevant harm or unfavorable direction blocks advancement of that proposition. |
| **Inconclusive/technically unqualified** | Resolve exposure, assay validity or precision; do not equate uncertainty with efficacy or falsification. |

This decision structure focuses resources on interpretable results. Laboratory costs and timelines should be determined with the executing team after model and material availability are known; neither is estimated from the computational runtime.

<!-- PAGEBREAK -->

## 18. A working system with inspectable execution

### 18.1 Completed engineering and provenance results

The recorded case execution completed the native core, downstream V3.1 research stage and separate verification with exit code zero. All **12 implemented verification checks passed**. The verifier inspected **211 direct receipt records representing 107 distinct API responses**, found no orphan summary references and no receipt failing its `store=false` consistency check. The project also records **280 passed synthetic automated tests**. [T1–T2,T8]

These results establish successful execution under the recorded environment and the behavior checked by those tests. Verification was performed by a separate project program. External reproducibility, biological attribution and predictive validity are addressed through independent runs, functional experiments and future benchmarks.

### 18.2 The recorded run

| Item | Recorded description |
|---|---|
| Environment | Existing Windows 11 workstation; a new project installation, not an externally operated virgin-machine replication. |
| Scientific execution | 22 October 2026, **16:00:05.781–17:48:20.723 UTC**. |
| Elapsed scientific time | **6,494.94 seconds**, approximately **108.25 minutes**. |
| Time accounting | Initial installation and untracked development/review effort excluded. |
| Cost accounting | No reconciled invoice or complete token-cost total available; no monetary estimate asserted. |

### 18.3 Reviewer-operated reproduction

The repository README is the operational authority. Fully extract the source ZIP and launch `JURY_START_V31.cmd` with the documented Windows permissions. The guided sequence is **1 → 2 → 3 → 6 → 7 → 4**: install, set the private API credential, check models, check readiness, run self-tests and start a new authorized-case analysis. [T8]

The reviewer supplies an authorized single-proband GRCh38 VCF, its matching structured DOCX and a billable API account. The installer creates or verifies a project-specific Ubuntu 24.04/WSL2 research environment and hash-pinned dependencies. It does not require a copy of the developer’s WSL instance.

Initial scientific-storage selection requires **at least 340 GiB free**, comprising a configured 240 GiB reference-resource budget and a 100 GiB protected reserve. This is an admission requirement, not a measurement or guarantee of download size. The package excludes patient files, credentials and bulk third-party databases, except for the attributed SIGNOR release seed.

### 18.4 What reproduction should preserve

The recorded profile includes ChEMBL 37, gnomAD 4.1, GENCODE 50 and VEP 116. Source identities, versions, hashes and receipts support traceability; “current official-source profile” does not mean every dependency automatically follows the newest release. Schema or integrity mismatches can stop a run rather than silently alter assumptions.

Future source updates and stochastic outputs may change rankings. Reproduction should preserve the method and evidence trail, not force a predetermined drug list. Exact snapshot repeatability and independent installation remain distinct evaluation questions.

<!-- PAGEBREAK -->

## 19. Scalability: what transfers and how to test it

### 19.1 A reusable unit of research

The reusable contribution is the sequence **authorized case → molecular hypotheses → directional objectives → multi-channel discovery → attributable evidence → qualified research decision**. Modular adapters and structured records allow the same operations to be applied to new compatible inputs without treating a disease name as sufficient pharmacological context. [T8]

The demonstrated scale is one completed case with thousands of catalogue entities and explicit evidence states. Gene-specific interpretation and laboratory qualification remain scientific work for each new case; the architecture supports that work rather than replacing it.

### 19.2 Transferable components and case-dependent work

| Reusable component | What must be requalified for another case |
|---|---|
| Regulatory normalization and medicine identity | Relevant jurisdiction, product, formulation and current source status. |
| Local genomic annotation and phenotype representation | Input quality, inheritance model, transcript context and causal relevance. |
| Signed-route construction and discovery adapters | Biological direction, network context and coverage for the new gene. |
| Source-bound evidence objects and critique states | Study-to-case transfer, negative evidence and missing literature. |
| Execution receipts and integrity checks | Environment, data versions, tool availability and permissions. |
| The staged experimental decision template | Model, functional endpoint, exposure and disease-specific safety. |

The three-gene downstream cap, finite retrieval budgets, source gaps, storage needs and billed reasoning are practical constraints to measure. They are not a claim of unlimited parallel throughput or established cost efficiency.

### 19.3 A prospective evaluation that would demonstrate scientific value

A useful next benchmark would freeze code, source snapshots and evaluation rules, then test multiple independently curated cases spanning different mechanism classes. Include known actionable relationships, biologically wrong-direction candidates, unsupported gene associations and cases where an appropriate result is unresolved or negative.

Compare the full pipeline with a simple gene–drug lookup and with selected channel ablations. Measure recovery of justified candidates, precision of source attribution, harmful-direction errors, missed counterevidence, appropriate uncertainty reporting, runtime, storage, API cost and reviewer effort. Avoid allowing the same evidence used to construct a test case to count as independent validation of its answer.

Expert adjudication should evaluate whether the proposed intervention follows from the case mechanism, not only whether a known drug name appears near the top. This directly tests the distinction that matters in this submission: high-quality pharmacology elsewhere is useful only after its relevance to the case is established.

### 19.4 A defined next development step

The current ranking records gene-linked routes, but does not establish a calibrated measure of each gene’s causal contribution to the disease. A future evaluation should therefore assess gene–case relevance and drug–mechanism evidence as separate dimensions. This is a prospective improvement, not a change made to the submitted software or its historical scores.

The present package already supplies the structured outputs and provenance needed to conduct that evaluation. Its scalability contribution is an executable and inspectable foundation for reuse; the proposed benchmark will measure performance across cases.

<!-- PAGEBREAK -->

## 20. Scientific contribution and conclusion

### Scientific rigor — a complete chain that can be inspected

The submission connects authorized genomic analysis to directional intervention objectives, records the discovery funnel, explains why alternatives did not become final nominations and anchors each proposed medicine to an identifiable experiment. The expanded LZTR1–RIT1 appraisal gives the modifier hypothesis a specific biological connection and defines its decisive genetic test. Competing mechanisms are converted into discriminating measurements rather than left as generic cautions. [5–7,39–43; T1–T9]

### Potential impact — a focused next research decision

If LZTR1 correction improves a qualified phenotype while BUB1B remains unchanged, that would support an additional contribution worth investigating in the model. If a candidate then improves useful function at justified exposure, it would provide a basis for independent preclinical follow-up. A precise contrary result would also focus the research effort. The immediate deliverable is this decision-ready program; improved patient outcomes, time savings and clinical utility remain prospective.

### Innovation — integrate discovery with causal interpretation

The creative contribution is the implemented combination of patient-first objective generation, several evidence channels, source-level critique and technical verification, followed by a clearly labeled post-run synthesis linking a candidate axis to the checkpoint. The individual annotation, network and fingerprint methods are established tools. The project’s contribution is their execution as an accountable research workflow and the specific experimental questions it produces, not a claim of inventing the cited biology.

### Scalability — reusable structure with measurable boundaries

Code, structured evidence, modular resources and a reviewer-operated installation route provide a transferable foundation. The executed case demonstrates operation at catalogue scale; independent multi-case performance and resource efficiency are the next validation targets. The report specifies how to test those claims rather than inferring them from software size. [T8]

### Integrated conclusion

**We nominate bortezomib and trametinib for separate, mechanism-qualified preclinical investigation.** Bortezomib tests the relationship between recoverable LZTR1 abundance and net functional regulation, including the competing RIT1-turnover effect. Trametinib tests downstream compensation where a measured signaling abnormality contributes to a relevant phenotype. The common first step is to establish the LZTR1 contribution against the BUB1B disease background.

The completed work comprises an executed case-to-candidate system, two source-linked propositions, documented alternatives and an expanded causal/functional evaluation plan. The strict `NO_VALIDATED_LEAD` outcome, exposure and safety status, finite review coverage and one-case evidence base remain visible because they define the present research stage. The work therefore supports hypothesis-driven investigation at its present stage, with a clear path toward stronger evidence.

**The request to the judging panel is to consider the quality of that scientific contribution and its suitability for independent follow-up.** The work offers a reproducible route from a difficult rare-disease case to experiments that can establish whether the proposed pharmacological opportunities are real.
<!-- PAGEBREAK -->

## Appendix A. Methods-template response and AI disclosure

### A.1 Integrated responses to the Track 2 methods questions

| Requested item | Response and location |
|---|---|
| **Team and responsibility** | ExaMedix Technologies; Matteo Caprani, research lead and developer. |
| **Detailed approach** | Sections 1–13: authorized case, mechanism hypotheses, rescue objectives, multi-channel discovery, selection and candidate-specific reasoning. |
| **Automated versus manual work** | Discovery and recorded ranking were automated. Post-run source reading, numerical reconciliation, audit interpretation and report preparation were AI-assisted and curated under author responsibility. |
| **Literature/expert curation** | Primary-source review, model/context separation and correction of retrieval omissions. No independent clinical/domain-specialist endorsement or new wet-laboratory work is claimed. |
| **Public versus proprietary evidence** | Public pharmacology, regulatory and scientific resources; organizer-controlled case inputs; paid computation, not an undisclosed proprietary drug-evidence database. |
| **Variant mechanism** | Sections 1–2, 7–9: transcript/product hypotheses and conditional functional attribution; exact controlled case details excluded. |
| **Effort** | Recorded scientific execution 6,494.94 seconds. Installation and total human effort not completely tracked; no reconciled monetary total. |
| **Abstract, strengths and boundaries** | Page 2 abstract; Sections 5–20. |
| **Reproduction and supporting material** | Section 18, Appendix D and repository README. The three-minute pitch remains a separate submission deliverable. |

The report integrates the official methods questions rather than treating the scientific narrative and methods disclosure as unrelated documents. [3–4]

### A.2 Providers, tiers and recorded settings

**Development and report assistance.** The research lead reports use of **OpenAI ChatGPT Pro**, with “Improve the model for everyone” disabled throughout relevant project use. This is an author-confirmed historical account setting, not a provider-side audit or a claim of enterprise terms. ChatGPT assisted engineering, source review, authorized artifact inspection and writing.

**Executable research.** The pipeline used commercially billed, pay-as-you-go **OpenAI Responses API** calls. Receipts record `store=false` for 107 distinct responses. Observed model strings were `gpt-5.4-mini-2026-03-17`, `gpt-5.6-terra`, `gpt-6-luna` and `gpt-6.1-sol`; these are recorded identifiers, not independently verified claims about model capabilities. Deterministic annotation, queries, counts, fingerprints and verification remain distinct local operations. [T2,T8]

**This editorial revision.** ChatGPT-assisted review used connected GitHub access for the latest report and public-web retrieval, including Firecrawl, for published sources. The supplied sanitized audit and addendum informed the rewrite. No new research run or paid pipeline analysis was initiated for the revision. No raw case genome or clinical document was provided to those public-web retrieval tools.

**Accountability.** AI output is an aid to research, not experimental evidence or external peer review. The report identifies the original studies and preserved artifacts supporting its claims; responsibility for the submitted interpretation remains with the author.

<!-- PAGEBREAK -->

## Appendix B. Privacy, controlled data and research governance

### B.1 Local case analysis and model-facing information

The implementation is designed to omit raw VCF content, exact case alleles/coordinates/HGVS, genotypes, identifiers and complete clinical documents from model-facing payloads. It uses sanitized biological abstractions and attributable public evidence. Whole-chromosome reference downloads avoid patient-specific coordinate URLs. [T2–T3,T8]

Those abstractions can still be derived information. Logged flags and receipt checks are not a comprehensive network capture or a full audit of every payload and connected service. The public report therefore describes the architecture and observed checks, without claiming a stronger privacy certification than the evidence supports.

The source package contains no participant files or API credentials. Authorized inputs and controlled derivatives must remain within the applicable access arrangements, separately from publicly distributed code and report materials.

### B.2 Training controls, storage and retention are distinct

OpenAI API training defaults and the `store=false` response setting should be distinguished from retention obligations. `store=false` is not a contractual Zero Data Retention arrangement; applicable abuse-monitoring retention and exceptions may remain. The author’s consumer-account training opt-out is a separate setting. No special enterprise agreement, provider-side deletion audit or absence of optional feedback/safety-review exceptions is asserted. [9–11]

Service and connector terms should be assessed for the actual account and material handled. For the present revision, public scientific pages and the report text were retrieved for literature/editorial work. The historical pipeline declarations and the current editorial workflow are documented separately rather than merged into a single blanket claim.

### B.3 Controlled-data obligations

The organizer’s dataset access conditions continue to govern restricted input and derivative material. The team’s CC BY license does not authorize redistribution of the genome, source clinical document or controlled intermediate data. Restrictions on reidentification, recontact and publication remain applicable. [2,30]

The recorded rules require deletion of controlled datasets and applicable intermediate/derived copies within the **30-day post-close period**, including relevant archives and backups, followed by confirmation to **cnb@sagebase.org**. This report states the obligation; it does not certify a deletion that has not occurred. Retention of sanitized submission artifacts must follow organizer guidance.

Public release for review should expose the code and permitted methods while keeping restricted material outside the repository. Report, methods, repository access and the separate pitch video must satisfy the organizers’ submission requirements. [2–4]

### B.4 Research and clinical responsibilities

No new specimen acquisition, intervention in the child or return of a research interpretation as a clinical result is authorized by this submission. Any laboratory program requires the appropriate institutional arrangements, consent scope and specialist oversight. The existing care team retains responsibility for clinical decisions.

The report’s public findings are limited to what is necessary to explain the research. The modifier hypothesis is offered for functional investigation; diagnosis, variant reclassification and clinical decisions require their respective validation processes.

<!-- PAGEBREAK -->

## Appendix C. Resource rights and reuse

### C.1 Original work and third-party terms

Original code and documentation are offered under **CC BY 4.0**, with attribution to Matteo Caprani / ExaMedix Technologies, the repository, preservation of the license and identification of modifications. The author declares a **non-commercial scientific contribution to the Challenge**. That purpose statement does not add a non-commercial restriction to the original-work license or relicense third-party resources. [2,12]

| Resource | Role in the investigation | Recorded rights treatment |
|---|---|---|
| **ChEMBL 37** | Pharmacology, assays, mechanisms, structures. | CC BY-SA 3.0 data; external bulk resource, not original CC BY project content. [14] |
| **SIGNOR October 2026** | Signed relationships and official seed. | CC BY 4.0; source/version attribution and original/packed hashes retained. [15] |
| **Enrichr/LINCS/GEO** | Public perturbational signatures. | Service and upstream-library terms; libraries not bundled or relicensed. [16–17,36–38] |
| **Ensembl VEP** | Offline consequence annotation. | Apache-2.0 software; reference data retain separate terms. [18] |
| **GENCODE/GRCh38; gnomAD** | Transcript/reference and population context. | Official source-specific policies, attribution and external acquisition. [19–21] |
| **ClinVar; HPO** | Variant assertions and phenotype concepts. | NCBI/source policies and HPO’s own license explanation. [22–23] |
| **UniProt; AlphaFold DB** | Protein and structural context. | Source CC BY 4.0 notices and attribution. [24–25] |
| **Open Targets; Reactome; DGIdb** | Target/pathway and discovery context. | Platform/database CC0 notices where applicable; software, illustrations and aggregator/source rights remain distinct. [26–28] |
| **FDA; EMA; MHRA; Health Canada** | Regulatory eligibility and normalization. | Agency/public-data, reproduction and Open Government terms; third-party material can retain separate rights. [32–35] |

### C.2 Enrichr classification remains a documented governance item

The executed adapter downloaded public signature libraries for local gene/direction matching; it did not submit the participant genome to Enrichr. The provider describes downloadable libraries and attribution, while its terms distinguish free academic/non-profit use from commercial use requiring a license. [16–17,36–38]

No written provider determination specifically covering this company-affiliated, prize-bearing entry is documented in the reviewed materials. The author’s declared research purpose does not assert academic-institution or registered-non-profit status. Provider clarification is the route to resolve that classification; later commercialization or redistribution requires its own review. Neither infringement nor unconditional authorization is inferred here.

### C.3 Dependencies and capability boundaries

RDKit retains BSD terms; NumPy, Pillow and Python retain their upstream BSD-family, MIT-CMU and PSF terms. Git, Ubuntu/WSL and other components retain their own licenses. Recorded pins include NumPy 2.5.3, Pillow 12.3.0 and RDKit 2026.3.6. Pangolin was acquired separately under GPL-3.0 but produced **zero predictions**; SpliceAI and AlphaMissense were not used in the executed default result. No qualified PK/PK-Sim model contributed. [13,29; T7–T8]

This inventory preserves the source-specific rights review. Dependency-level clearance and provider permissions remain subject to their applicable terms; the metadata supports that review rather than substituting for it.

<!-- PAGEBREAK -->

## Appendix D. Provenance, version boundaries and acknowledgment

### D.1 Preserved technical evidence

| ID | Artifact or inspected material | Evidence represented |
|---|---|---|
| **T1** | `E2E_EXECUTION_V31.json` | Timestamps, execution stages and exit codes. |
| **T2** | `CLEANROOM_REPRODUCIBILITY_ATTESTATION_V31.json` | Twelve checks, model identifiers and response-receipt accounting. |
| **T3** | `GENOMIC_REQUIRED_EVIDENCE_V31.json`; patient-state inspection | Required VEP/gnomAD execution and broad consequence hypotheses. |
| **T4** | `FINAL_RESULT_V31_NATIVE.json`; screening and pharmacology summaries | Catalogue, routes, signals, profiles and strict outcome. |
| **T5** | `TRACK2_HYPOTHESES_V31_FINAL.json`; generated dossier | Candidate states, scores, final nominees and combination result. |
| **T6** | `SOURCE_BOUND_ATTACK_SUMMARY.json` and source inspection | Accepted/rejected evidence, source identities and critique context. |
| **T7** | Splice-prediction and exposure/PK summaries | Module states and outstanding qualification. |
| **T8** | Source, README, LICENSE, requirements, notices and tests | Implemented method, installation, rights and 280 synthetic tests. |
| **T9** | `MVA_Track2_BUB1B_LZTR1_ReadOnly_Audit.json` | Supplied sanitized read-only audit of retained BUB1B and final LZTR1 decisions. |

T9 concerns the preserved run `JURY_V31_20261022_175934`. The report’s original execution evidence and this supplied audit are distinct from the current editorial work. Not every underlying artifact is a permissible public redistribution asset.

### D.2 Integrity and revision record

**Recorded final scientific-summary SHA-256**  
`411fb311d4bea3459e5c483b1116c1af94ecd439fcaac29601aa13beaf76f067`

**Recorded final hypothesis-file SHA-256**  
`d4bacb1178d4c0700d931699c7b66a90d59d7d4257c8e33ca0914e4269a4c115`

**Recorded official-source snapshot SHA-256**  
`cf5bed60a0d16afb10793c615d196b4e8bdadc2e90358bea4de1bd5c8ebdf41f`

**Recorded methods-workbook SHA-256**  
`61aab080a2868a3b724e76692b83c24812112e305cd3a8b03f8f91a6b2414441`

The four SHA-256 values above identify the recorded scientific outputs and source materials. Internal development identifiers are retained in a private audit, separate from this public submission. The LZTR1–RIT1 literature synthesis and expanded validation plan were added after the recorded run; the original software, results, evidence labels and candidate priorities remain unchanged.

### D.3 Required acknowledgment

> This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

Matteo Caprani / ExaMedix Technologies is the research lead and developer. Author responsibility and AI assistance are disclosed in Appendix A. No organizer endorsement or independent specialist approval is implied.

<!-- PAGEBREAK -->

## References · 1–23

The new mechanistic appraisal uses primary studies [39–43] alongside the preserved candidate sources [6–7]. Reference [8] is explicitly a secondary source retained to explain the original audit’s conflict state. Resource and governance references preserve the prior report’s source-specific inventory; the revision does not claim that every historical terms page was independently re-audited. Technical evidence [T1–T9] is listed in Appendix D.

**[1]** Sage Bionetworks. *MVA Hackathon 2026: Track 2 objective and weighted judging criteria.* Official About tab. https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py

**[2]** Sage Bionetworks. *Official Rules: data, licensing, acknowledgment and publication/deletion obligations.* [Official rules source](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py).

**[3]** Sage Bionetworks. *Track 2 submission instructions.* [Official submission source](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track2.py).

**[4]** Sage Bionetworks. *Methods description workbook, Track 2.* [Official methods workbook](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/resolve/main/static/templates/methods_description_form.xlsx).

**[5]** Hanks S, et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. *Nature Genetics.* 2004;36:1159–1161. DOI: [**10.1038/ng1449**](https://doi.org/10.1038/ng1449). PMID: 15475955.

**[6]** Bacchiocchi A, et al. LZTR1 is a melanoma oncogene that promotes invasion and suppresses apoptosis. *Oncogene.* 2025;44:3974–3984. DOI: [**10.1038/s41388-025-03538-2**](https://doi.org/10.1038/s41388-025-03538-2). PMID: 40885854; PMCID: PMC12500468. Relevant: USP9X/LZTR1 experiments, Figure 3.

**[7]** Abe T, Morisaki K, Niihori T, Terao M, Takada S, Aoki Y. Dysregulation of RAS proteostasis by autosomal-dominant LZTR1 mutation induces Noonan syndrome-like phenotypes in mice. *JCI Insight.* 2024;9(22):e182382. DOI: [**10.1172/jci.insight.182382**](https://doi.org/10.1172/jci.insight.182382). PMCID: PMC11601938. Relevant: Figure 6 and associated results.

**[8]** Saint-Laurent C, Mazeyrie L, Yart A, Edouard T. Novel therapeutic perspectives in Noonan syndrome and RASopathies. *European Journal of Pediatrics.* 2024;183:1011–1019. DOI: [**10.1007/s00431-023-05263-y**](https://doi.org/10.1007/s00431-023-05263-y). PMCID: PMC10951041. Secondary source for the recorded model-dependent counterevidence.

**[9]** Sage Bionetworks. *Third-party LLM API and corrected consumer-tier guidance.* Official Space discussion 2. https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2

**[10]** OpenAI. *API data controls: training, response state and retention.* https://developers.openai.com/api/docs/guides/your-data

**[11]** OpenAI. *Data Controls FAQ.* https://help.openai.com/en/articles/7730893-data-controls-in-chatgpt

**[12]** Creative Commons. *Attribution 4.0 International legal code.* https://creativecommons.org/licenses/by/4.0/legalcode.en

**[13]** RDKit contributors. *Toolkit overview and licensing.* https://www.rdkit.org/docs/Overview.html

**[14]** ChEMBL / EMBL-EBI. *General FAQ and data-reuse licensing.* https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions

**[15]** SIGNOR consortium. *Documentation and October 2026 release.* https://signor.uniroma2.it/documentation/ ; https://signor.uniroma2.it/releases/Oct2026_release.txt

**[16]** Ma’ayan Lab. *Enrichr terms.* https://maayanlab.cloud/Enrichr/templates/help/terms-submenu.html

**[17]** Chen EY, et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. *BMC Bioinformatics.* 2013;14:128. DOI: [**10.1186/1471-2105-14-128**](https://doi.org/10.1186/1471-2105-14-128).

**[18]** Ensembl. *VEP documentation and upstream software.* https://www.ensembl.org/info/docs/tools/vep/ ; https://github.com/Ensembl/ensembl-vep

**[19]** GENCODE. *Official data access and human releases.* https://www.gencodegenes.org/pages/data_access.html ; https://www.gencodegenes.org/human/

**[20]** EMBL-EBI. *Terms for data resources and services.* https://www.ebi.ac.uk/about/terms-of-use

**[21]** gnomAD / Broad Institute. *Version 4.1 and policies.* https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/ ; https://gnomad.broadinstitute.org/policies

**[22]** NCBI. *ClinVar and source policies.* https://www.ncbi.nlm.nih.gov/clinvar/ ; https://www.ncbi.nlm.nih.gov/home/about/policies/

**[23]** Human Phenotype Ontology. *License and explanation.* https://hpo.jax.org/app/license ; https://github.com/obophenotype/human-phenotype-ontology/blob/master/LICENSE.md

<!-- PAGEBREAK -->

## References · 24–43

**[24]** UniProt. *License and disclaimer.* https://www.uniprot.org/help/license

**[25]** AlphaFold Protein Structure Database / EMBL-EBI. *FAQ and terms.* https://alphafold.ebi.ac.uk/faq

**[26]** Open Targets. *Platform data licence.* https://platform-docs.opentargets.org/licence

**[27]** Reactome. *Data, software and illustration terms.* https://reactome.org/license

**[28]** DGIdb. *Resource and upstream notices.* https://www.dgidb.org/ ; https://github.com/dgidb/dgidb-v5

**[29]** Pangolin authors. *Software and upstream license.* https://github.com/tkzeng/Pangolin

**[30]** Sage Bionetworks. *Rare Disease, Real Kid: MVA Hackathon 2026 — controlled dataset.* Synapse project **syn76251147**; controlled-access distribution: https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data

**[31]** ExaMedix Technologies. *Rare Disease MVA Pharma: source, README, methods and license.* https://github.com/mattc789/rare-disease-mva-pharma-track2

**[32]** U.S. Food and Drug Administration. *Website policies and third-party rights.* https://www.fda.gov/about-fda/contact-fda/website-policies

**[33]** European Medicines Agency. *Legal notice and reproduction.* https://www.ema.europa.eu/en/about-us/about-website/legal-notice

**[34]** UK National Archives. *Open Government Licence v3.0.* https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/

**[35]** Health Canada. *Drug Product Database — All Files.* https://open.canada.ca/data/en/dataset/bf55e42a-63cb-4556-bfd8-44f26e5a36fe

**[36]** Ma’ayan Lab. *Enrichr downloadable libraries and attribution.* https://maayanlab.cloud/Enrichr/

**[37]** Kuleshov MV, et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Research.* 2016. DOI: [**10.1093/nar/gkw377**](https://doi.org/10.1093/nar/gkw377).

**[38]** Xie Z, et al. Gene set knowledge discovery with Enrichr. *Current Protocols.* 2021;1:e90. DOI: [**10.1002/cpz1.90**](https://doi.org/10.1002/cpz1.90).

**[39]** Castel P, Cheng A, Cuevas-Navarro A, et al. RIT1 oncoproteins escape LZTR1-mediated proteolysis. *Science.* 2019;363:1226–1230. DOI: [**10.1126/science.aav1444**](https://doi.org/10.1126/science.aav1444). PMID: **30872527**; PMCID: **PMC6986682**. Relevant: Figure 3 (LZTR1/CUL3-dependent turnover and bortezomib) and Figure 4 (cell models and Noonan-family fibroblasts).

**[40]** Cuevas-Navarro A, Van R, Cheng A, Urisman A, Castel P, McCormick F. The RAS GTPase RIT1 compromises mitotic fidelity through spindle assembly checkpoint suppression. *Current Biology.* 2021;31:3915–3924.e9. DOI: [**10.1016/j.cub.2021.06.030**](https://doi.org/10.1016/j.cub.2021.06.030). PMID: **34237269**; PMCID: **PMC8440430**. Relevant: Figures 3–4 and related supplementary results for RIT1/LZTR1, mitotic timing, chromosome errors and checkpoint-complex integrity.

**[41]** Suijkerbuijk SJE, van Osch MH, Bos FL, Hanks S, Rahman N, Kops GJPL. Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. *Cancer Research.* 2010;70:4891–4900. DOI: [**10.1158/0008-5472.CAN-09-4319**](https://doi.org/10.1158/0008-5472.CAN-09-4319). PMID: **20516114**; PMCID: **PMC2887387**.

**[42]** North BJ, et al. SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. *EMBO Journal.* 2014;33:1438–1453. DOI: [**10.15252/embj.201386907**](https://doi.org/10.15252/embj.201386907). PMID: **24825348**; PMCID: **PMC4194088**. Distinguish NMN exposure, genetic SIRT2 elevation and nicotinamide-inhibition experiments.

**[43]** Sieben CJ, et al. BubR1 allelic effects drive phenotypic heterogeneity in mosaic-variegated aneuploidy progeria syndrome. *Journal of Clinical Investigation.* 2020;130:171–188. DOI: [**10.1172/JCI126863**](https://doi.org/10.1172/JCI126863). PMID: **31738183**. Relevant: allele/tissue differences and the selected muscle mTORC1 observations.
