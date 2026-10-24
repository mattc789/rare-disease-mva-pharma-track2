"""Synthetic full-profile genomics gate: VEP and case-selected gnomAD are mandatory."""
import unittest
from v30.native_pipeline import verify_required_genomic_evidence

def good():
    private={"genes":{"TEST":[{"chrom":"7","pos":100,"ref":"A","alt":"G"}]}}
    vep_setup={"state":"READY"}
    vep_status={"tool_state":"COMPLETED","annotation_records_returned":1,
                "tool_provenance":{"run_command_kind":"wsl_local_perl_offline"}}
    gnomad_audit={"completion":"READY","required_jobs":2,"installed_jobs":2}
    gnomad_status={"state":"LOCAL_PARTIAL","v4_1_full_coverage":False,
                   "cohort_chromosome_completeness":{"exomes":["7"],"genomes":["7"]}}
    rows=[{"measurements":[{"status":"NOT_OBSERVED_IN_THIS_INDEXED_LOCAL_SHARD"},
                            {"status":"LOCALLY_MATCHED"}]}]
    return private,vep_setup,vep_status,gnomad_audit,gnomad_status,rows

class V31GenomeAuditTests(unittest.TestCase):
    def test_case_focused_coverage_passes_without_falsely_claiming_global(self):
        a=verify_required_genomic_evidence(*good())
        self.assertTrue(a["required_science_pass"])
        self.assertFalse(a["checks"]["gnomad_full_24_chromosome_both_cohorts"])
        self.assertEqual(a["case_chromosome_cohort_count"],2)

    def test_absent_genome_shard_blocks_the_full_profile(self):
        a=list(good())
        a[4]["cohort_chromosome_completeness"]["genomes"]=[]
        self.assertFalse(verify_required_genomic_evidence(*a)["required_science_pass"])

    def test_hash_complete_but_query_failed_is_not_valid(self):
        a=list(good())
        a[5][0]["measurements"][0]["status"]="INDEX_QUERY_FAILED"
        self.assertFalse(verify_required_genomic_evidence(*a)["required_science_pass"])

    def test_vep_installed_but_not_run_blocks_profile(self):
        a=list(good())
        a[2]["tool_state"]="UNAVAILABLE_OR_FAILED"
        self.assertFalse(verify_required_genomic_evidence(*a)["required_science_pass"])

    def test_vep_zero_annotations_blocks_profile(self):
        a=list(good())
        a[2]["annotation_records_returned"]=0
        self.assertFalse(verify_required_genomic_evidence(*a)["required_science_pass"])

    def test_gnomad_incomplete_audit_blocks_profile(self):
        a=list(good())
        a[3]["completion"]="INCOMPLETE_ERROR"
        self.assertFalse(verify_required_genomic_evidence(*a)["required_science_pass"])

if __name__=="__main__":unittest.main()
