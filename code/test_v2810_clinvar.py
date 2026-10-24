"""ClinVar no_conflicts regression: rare-disease V2.8.7, synthetic records only."""
import unittest
from unittest.mock import patch
from collections import defaultdict
from v2810 import case_science,intake

class ClinVarConcordance(unittest.TestCase):
    REV_OK="criteria_provided,_multiple_submitters,_no_conflicts"
    def test_pathogenic_consensus_regains_proper_weight(self):
        raw={"CLNSIG":"Pathogenic/Likely_pathogenic","CLNREVSTAT":self.REV_OK}
        result=case_science.clinvar_interpretation(raw)
        self.assertEqual(result["state"],"reported_pathogenic")
        self.assertEqual(result["weight"],6)
        self.assertEqual(intake._clinvar_weight(raw),6)
    def test_benign_consensus_is_vetoed_from_primary_ranking(self):
        raw={"CLNSIG":"Benign/Likely_benign","CLNREVSTAT":self.REV_OK}
        result=case_science.clinvar_interpretation(raw)
        self.assertEqual(result["state"],"reported_benign")
        self.assertEqual(result["weight"],-10)
        self.assertEqual(intake._clinvar_weight(raw),-10)
    def test_genuine_conflicts_not_mislabeled_as_concordant(self):
        for review in ("criteria_provided,_conflicting_classifications",
                       "criteria_provided,_conflicts"):
            result=case_science.clinvar_interpretation(
                {"CLNSIG":"Pathogenic/Likely_pathogenic","CLNREVSTAT":review})
            self.assertEqual(result["state"],"conflicting")
            self.assertEqual(result["weight"],0)
    def test_unreviewed_classification_is_not_overweighted(self):
        raw={"CLNSIG":"Pathogenic","CLNREVSTAT":"no_assertion_criteria_provided"}
        self.assertEqual(intake._clinvar_weight(raw),0)
    def test_fully_benign_gene_cannot_win_with_homozygous_bonus(self):
        def variant(g,gt,cv):
            return {"gene":g,"GT":gt,"filter":"PASS","GQ":99,
             "allele_balance":1.0 if gt=="1/1" else 0.5,"clinvar":cv,
             "best_annotation":{"severity":11,"consequence":"splice_donor_variant",
                                "reference_sequence_mismatch":False},
             "chrom":"1","pos":1,"ref":"A","alt":"G","phase_tags":{}}
        rows=[variant("GENE_BENIGN","1/1",
               {"CLNSIG":"Benign/Likely_benign","CLNREVSTAT":self.REV_OK}),
              variant("GENE_UNKNOWN","0/1",None)]
        with patch.object(intake,"parse_hpo_resources",
                          return_value=(defaultdict(set),defaultdict(set))):
            ranked=intake.rank_genes(rows,[])
        self.assertEqual([x["gene"] for x in ranked],["GENE_UNKNOWN"])
    def test_single_homo_call_does_not_phase_other_heterozygotes(self):
        mixed=[{"GT":"1/1"},{"GT":"0/1","phase_tags":{}}]
        state=case_science.phase_assessment(mixed)
        self.assertEqual(state["state"],"homozygous_with_other_unphased_variants")
        self.assertFalse(state["trans_alleles_established"])
        self.assertEqual(case_science.phase_assessment([{"GT":"1/1"}])["state"],"homozygous_observed")
if __name__=="__main__":unittest.main()
