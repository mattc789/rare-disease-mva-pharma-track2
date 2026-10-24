# Rare Disease MVA Pharma — Third-Party Notices (V3.1)

This file documents third-party source/tool terms used by the research pipeline. It is an engineering compliance record, not legal advice.

Policy version: 2026-10-23-v31-submission-disclosure

## Default Challenge/public profile

- Static source-policy gate: PASS for the declared default policy. This is **not** a certification that every optional scientific service actually used in a run has unconditional reuse rights; see the actual Enrichr-use disclosure below.
- Third-party bulk datasets are not silently relicensed under the project's own licence.
- Gated/patient data and API secrets are never redistribution assets.

## Sources

### AlphaFold_DB
- Class: GREEN_ATTRIBUTION
- Enabled by default: true
- Terms: AlphaFold DB prediction data CC BY 4.0; academic and commercial use explicitly permitted by AFDB. AF2 code Apache-2.0; model parameters CC BY 4.0.
- Redistribution strategy: Runtime retrieval of selected structure/PAE files only. Preserve AFDB/DeepMind/EMBL-EBI attribution, source URL, accession and retrieval date. Do not claim clinical validation or warranty; do not bundle unrelated mirrored databases.

### AlphaMissense
- Class: YELLOW_LICENSE_CONFLICT
- Enabled by default: false
- Terms: CONFLICT: current DeepMind/Ensembl/AlphaFold DB documentation states CC BY 4.0, but the actual official AlphaMissense_hg38.tsv.gz file header states CC BY-NC-SA 4.0 and an EMBL-EBI 2023 page also describes non-commercial terms.
- Redistribution strategy: Implemented as an optional external research adapter only. Disabled in Challenge/public/commercial baseline. Do not download, redistribute, embed scores, or use for a submission artifact until rights-holder clarification resolves the conflicting official notices.

### Official ChEMBL 37 SQLite source — default V3.1 scientific provider
- Source: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/chembl_37_sqlite.tar.gz
- Release: ChEMBL 37, database curated by ChEMBL / EMBL-EBI.
- License: Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0).
- Official LICENSE, README and SHA-256 checksums downloaded alongside the database.
- Only source code and connectors belong to the ExaMedix application release.
  The original and locally derived ChEMBL data are **not bundled** or
  relicensed under the application's CC BY source license.
- If an adapted data extract is redistributed, provide source, release,
  changes, license and share-alike notices. This software's source package
  does not include patient data or the ChEMBL database.

### ChEMBL
- Class: YELLOW_SHAREALIKE_DATA
- Enabled by default: true
- Terms: ChEMBL data CC BY-SA 3.0
- Redistribution strategy: Download separately and retain ChEMBL provenance/share-alike obligations. Do not relicense or bundle an adapted ChEMBL database inside the CC-BY source archive.

### ClinVar
- Class: GREEN_RUNTIME
- Enabled by default: true
- Terms: NCBI public resource; source assertions remain third-party submissions.
- Redistribution strategy: Official download/runtime use with provenance; do not present assertions as independently validated.

### DGIdb
- Class: YELLOW_DISCOVERY_INDEX
- Enabled by default: true
- Terms: DGIdb v5 software MIT. DGIdb aggregates claims from many upstream knowledgebases/publications with source-specific terms.
- Redistribution strategy: Use DGIdb only as a discovery/index channel. Persist source names/PMIDs. A DGIdb-only interaction cannot become submission evidence until verified in a primary/public source whose reuse terms are compatible.

### Enrichr_LINCS_GEO
- Class: CONDITIONAL_THIRD_PARTY_SERVICE_TERMS
- Observed completed Challenge run: **ENABLED** (4,042 perturbational-expression signal records); the executable V3.1 screening module uses `RDMVA_ENRICHR_MODE=enabled` unless explicitly disabled. The commercial-mode gate fails closed if this channel is enabled; the supported way to exclude it is to set `RDMVA_ENRICHR_MODE=disabled`. Any commercial integration would require permission and appropriate implementation review.
- Terms: [Ma'ayan Lab Enrichr official terms](https://maayanlab.cloud/Enrichr/templates/help/terms-submenu.html) distinguish free academic/non-profit research from commercial use requiring a separate license. ExaMedix describes this Challenge entry as research, not the sale of enrichment analysis or clinical services; the provider's classification of company-associated competition research has not been independently determined.
- Use and redistribution: The executed workflow retrieved public LINCS/GEO gene-set-library text through the Enrichr service and matched it locally. Its public GitHub source release does not bundle an Enrichr library or relicense its upstream datasets. Preserve source citations and upstream terms. Public downloadability alone is not a blanket commercial license.
- Policy discrepancy disclosure: The V3 source-policy declaration lists this channel as disabled by default, but the **actual V3.1 research run enabled it**; scientific descriptions must follow the recorded execution, not the default policy label. See the integrated final report for qualification.

### Ensembl_VEP
- Class: GREEN_EXTERNAL
- Enabled by default: true
- Terms: Apache-2.0
- Redistribution strategy: External tool installed separately from official upstream.

### HPO
- Class: YELLOW_ATTRIBUTION
- Enabled by default: true
- Terms: See source-specific policy record.
- Redistribution strategy: Official runtime download; preserve attribution/version and source-file integrity.

### Illumina_SpliceAI
- Class: RED
- Enabled by default: false
- Terms: See source-specific policy record.
- Redistribution strategy: Excluded from free/open baseline because upstream software/model terms do not match the required profile.

### OpenAI_API
- Class: GREEN_SERVICE
- Observed completed AI-enabled Challenge run: true; a fresh reviewer must separately provide their own billable API key and explicitly consent to the online requests.
- Terms: External API service; no API key or service code redistributed.
- Redistribution strategy: BYO API key only; store=false; sanitized research envelope; no raw VCF/genotypes/patient identifiers.

### Open_Targets
- Class: GREEN
- Enabled by default: true
- Terms: Platform marked CC0 1.0; code is open-source under repository-specific licences.
- Redistribution strategy: Prefer runtime API/download with provenance. Do not silently relicense embedded upstream evidence.

### PK-Sim
- Class: YELLOW_GPL_EXTERNAL
- Enabled by default: false
- Terms: GPL-2.0
- Redistribution strategy: Independent external program/adapter; model inputs need their own provenance.

### Pangolin
- Class: YELLOW_GPL_EXTERNAL
- Enabled by default: true
- Terms: GPL-3.0
- Redistribution strategy: Independent external process; do not relicense as CC-BY.

### Reactome
- Class: GREEN
- Enabled by default: true
- Terms: Database data and derived data files CC0; most software Apache-2.0; illustrations/branding CC BY 4.0; some code may have dependency-specific copyleft.
- Redistribution strategy: Use pathway/interaction data through runtime API or official downloads; preserve provenance.

### gnomAD
- Class: YELLOW_POLICY
- Enabled by default: true
- Terms: See source-specific policy record.
- Redistribution strategy: Official runtime download; no bulk database redistribution in source release.

## AlphaFold DB

- Default profile: allowed with attribution.
- Prediction data are treated as theoretical research context, not clinically validated evidence.
- Selected structures/metadata should be retrieved at runtime with accession/source attribution.

## AlphaMissense

- Disabled by default.
- Current project policy records conflicting official notices for the precomputed prediction artifact.
- It is not used for the Challenge/public/commercial baseline until the rights conflict is resolved.

## ChEMBL / DGIdb

- ChEMBL remains external/source-specific data with its own share-alike obligations; the source package does not bundle an adapted ChEMBL database.
- DGIdb is discovery/index only; an aggregator-only interaction cannot become submission evidence without primary-source verification.

### Official SIGNOR October 2026 CC BY 4.0 pinned fallback (V3.1.2)

The source package includes one 6,722,506-byte public losslessly compressed
SIGNOR 4.0 October 2026 snapshot solely to withstand an official university
server outage. The file decompresses to exactly 21,425,930 bytes
(SHA-256 63db6f65e9148e81354f564c3c3cdfc9ee45dc421553bd0926e5fae511ef9d06).
Original data remains licensed under CC BY 4.0. Full attribution and original
source URL are given in reference_seed/signor/ATTRIBUTION.md.
Its use is explicitly labeled a verified pinned snapshot, not a
freshly captured live download.
