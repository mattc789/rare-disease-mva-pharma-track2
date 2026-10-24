"""Adversarial research-tool checks for V2.8.9. All variants and drugs SYNTHETIC."""
import csv,hashlib,json,math,unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from v2810 import vep_local,gnomad_local,spliceai_local,quant_pk,specialist_evidence

def variant():
    return {"gene":"DEMO1","chrom":"1","pos":100,"ref":"A","alt":"G",
      "GT":"0/1","phase_tags":{},"hgvsp":"p.Ala2Gly"}
def private():
    return {"genes":{"DEMO1":[variant()]},
      "phenotype_features":[],"residual_protein_status":"unknown"}
def putcsv(p,rows):
    with Path(p).open("w",newline="",encoding="utf-8") as f:
      w=csv.DictWriter(f,fieldnames=["time_h","unbound_concentration_nM"])
      w.writeheader()
      for t,c in rows:w.writerow({"time_h":t,"unbound_concentration_nM":c})

class VEPChecks(unittest.TestCase):
 def test_vep_input_contains_no_genotype_or_patient_id(self):
    with TemporaryDirectory() as td:
      keys=vep_local.generate_vcf([variant()],Path(td)/"v.vcf")
      raw=(Path(td)/"v.vcf").read_text()
      self.assertEqual(len(keys),1)
      self.assertIn("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",raw)
      self.assertNotIn("0/1",raw)
      self.assertNotIn("DEMO1",raw)
      self.assertNotIn("SAMPLE",raw)
 def test_vep_build_and_allele_inputs_enforced(self):
    self.assertTrue(vep_local.validate_input(variant()))
    self.assertFalse(vep_local.validate_input({**variant(),"chrom":"HLA*01"}))
    self.assertFalse(vep_local.validate_input({**variant(),"ref":"ACGT*"}))
    self.assertFalse(vep_local.validate_input({**variant(),"pos":-4}))
 def test_vep_wsl_path_in_local_windows(self):
    self.assertEqual(vep_local.wsl_path(r"Z:\Data\private_file.vcf"),"/mnt/z/Data/private_file.vcf")
    with self.assertRaises(ValueError):vep_local.wsl_path("https://example.com/variant")
 def test_vep_defaults_are_derived_from_checkout_and_active_wsl_user(self):
    with patch.object(vep_local,"PROJECT_ROOT",Path(r"Z:\Portable Checkout")), patch.object(vep_local,"_wsl_home",return_value="/home/testuser"):
      d=vep_local.default_paths()
    self.assertEqual(d["vep"],"/home/testuser/rdmva_tools/ensembl-vep-116/vep")
    self.assertTrue(d["gtf"].startswith("/mnt/z/Portable Checkout/"))
    self.assertTrue(d["fasta"].startswith("/mnt/z/Portable Checkout/"))
 def test_vep_json_parser_finds_actual_transcript_terms(self):
    with TemporaryDirectory() as td:
      f=Path(td)/"out.json"
      j={"input":"1\t100\t.\tA\tG\t.\tPASS\t.","seq_region_name":"1","start":100,
       "most_severe_consequence":"splice_donor_variant",
       "transcript_consequences":[{"gene_symbol":"DEMO1","gene_id":"ENSG_TEST",
         "transcript_id":"ENST_TEST","consequence_terms":["splice_donor_variant"],
         "impact":"HIGH","protein_start":12}]}
      f.write_text(json.dumps(j)+"\n")
      rows=vep_local.parse_json_lines(f)
      self.assertEqual(rows[0]["most_severe_consequence"],"splice_donor_variant")
      self.assertEqual(rows[0]["transcript_annotations"][0]["gene_symbol"],"DEMO1")
      self.assertEqual(rows[0]["ref"],"A")
 def test_vep_indel_joins_on_original_vcf_pos_not_normalized_start(self):
    with TemporaryDirectory() as td:
      f=Path(td)/"indel.json"
      j={"input":"1\t100\t.\tAA\tA\t.\tPASS\t.",
        "seq_region_name":"1","start":101,
        "most_severe_consequence":"frameshift_variant"}
      f.write_text(json.dumps(j)+"\n",encoding="utf-8")
      row=vep_local.parse_json_lines(f)[0]
      self.assertEqual(row["pos"],100)
      self.assertEqual(row["ref"],"AA")
      self.assertEqual(row["alt"],"A")
 def test_vep_offline_stage_source_verified_by_real_output(self):
    with TemporaryDirectory() as td:
      def fake(inp,out,**kwargs):
        Path(out).write_text(json.dumps({
         "input":"1\t100\t.\tA\tG\t.\tPASS\t.",
         "seq_region_name":"1","start":100,
         "most_severe_consequence":"missense_variant"})+"\n")
        return {"exit_code":0,"run_command_kind":"wsl_local_perl_offline"}
      with patch.object(vep_local,"_run_vep",side_effect=fake):
        stat,rows=vep_local.annotation_stage(private(),td)
      self.assertEqual(stat["tool_state"],"COMPLETED")
      self.assertEqual(len(rows),1)
      self.assertFalse(stat["clinically_qualified_interpretation"])
 def test_vep_missing_dependencies_never_become_normal_biology(self):
    with TemporaryDirectory() as td,\
       patch.object(vep_local,"_run_vep",side_effect=FileNotFoundError("missing reference")):
      stat,rows=vep_local.annotation_stage(private(),td)
      self.assertEqual(stat["tool_state"],"UNAVAILABLE_OR_FAILED")
      self.assertEqual(rows,[])

class GnomADChecks(unittest.TestCase):
 def test_missing_full_gnomad_is_not_complete_or_zero(self):
    manifest,st=gnomad_local.validate_manifest(None)
    self.assertFalse(st["v4_1_full_coverage"])
    x=gnomad_local.lookup_allele(variant(),manifest)
    self.assertEqual(x["state"],"FREQUENCY_UNDETERMINED")
    self.assertIsNone(x["max_valid_global_af"])
 def test_gnomad_release_build_and_index_checked(self):
    with TemporaryDirectory() as td:
      p=Path(td)/"manifest.json"
      f=Path(td)/"1.vcf.gz"
      f.write_bytes(b"synthetic")
      Path(str(f)+".tbi").write_bytes(b"index-synthetic")
      m={"schema":gnomad_local.SCHEMA,"release":"v4.1","assembly":"GRCh38",
        "shards":[{"cohort":"exomes","chrom":"1","file":str(f),
          "official_reference":"SYNTHETIC_SOURCE"}]}
      p.write_text(json.dumps(m))
      shards,status=gnomad_local.validate_manifest(p)
      self.assertEqual(status["shards_found"],1)
      self.assertFalse(status["v4_1_full_coverage"])
      m["assembly"]="GRCh37";p.write_text(json.dumps(m))
      with self.assertRaisesRegex(ValueError,"MISMATCH"):gnomad_local.validate_manifest(p)
 def test_gnomad_multi_allele_AF_correct_index_and_AN(self):
    shards={("exomes","1"):{"file":"SOURCE.vcf.gz","reference":"gnomAD_v4.1_qualified"}}
    def query(*args):
      return ["1\t100\t.\tA\tT,G\t.\tPASS\tAF=0.2,0.0002;AC=100,1;AN=500"]
    x=gnomad_local.lookup_allele(variant(),shards,query=query)
    self.assertEqual(x["max_valid_global_af"],0.0002)
    self.assertEqual(x["measurements"][0]["allele_count"],1)
    self.assertEqual(x["measurements"][1]["status"],"DATASET_OR_CHROMOSOME_MISSING")
 def test_gnomad_no_match_not_frequency_zero(self):
    shards={("exomes","1"):{"file":"SOURCE.vcf.gz","reference":"gnomAD_v4.1_qualified"}}
    x=gnomad_local.lookup_allele(variant(),shards,query=lambda *args:[])
    self.assertIsNone(x["max_valid_global_af"])
    self.assertTrue(x["measurements"][0]["absence_is_not_proof_of_zero_frequency"])
 def test_gnomad_nonpass_cannot_produce_AF(self):
    shards={("genomes","1"):{"file":"SOURCE.vcf.gz","reference":"gnomAD_v4.1"}}
    x=gnomad_local.lookup_allele(variant(),shards,
       query=lambda *args:["1\t100\t.\tA\tG\t.\tRF\tAF=0.10;AC=100;AN=1000"])
    self.assertIsNone(x["max_valid_global_af"])
    self.assertIn("source_filter_not_pass",x["measurements"][1]["record_flags"])

class SpliceAIChecks(unittest.TestCase):
 def test_illumina_requires_explicit_commercial_license(self):
    _,st=spliceai_local.validate_manifest(None)
    self.assertFalse(st["commercial_use_automatically_permitted"])
    with TemporaryDirectory() as td:
      f=Path(td)/"s.vcf.gz";f.write_bytes(b"x")
      Path(str(f)+".tbi").write_bytes(b"x")
      p=Path(td)/"attest.json"
      m={"schema":spliceai_local.SCHEMA,"assembly":"GRCh38",
        "licensee":"ExaMedix Technologies","license_scope":"academic",
        "commercial_use_permitted":False,"license_reference":"STUB",
        "score_files":[{"file":str(f),"assembly":"GRCh38","source":"Illumina_SpliceAI"}]}
      p.write_text(json.dumps(m))
      _,s=spliceai_local.validate_manifest(p)
      self.assertEqual(s["state"],"LICENSE_NOT_ATTESTED")
 def test_exact_alternate_allele_and_gene_score_parsing(self):
    info="T|DEMO1|0.99|0|0|0|0|0|0|0,G|OTHER|0.9|0|0|0|0|0|0|0,G|DEMO1|0.12|0.51|0.03|0.01|2|-5|0|0"
    rows=spliceai_local.parse_score_text(info,"G","DEMO1")
    self.assertEqual(len(rows),1)
    self.assertEqual(rows[0]["maximum_delta_score"],0.51)
    self.assertEqual(rows[0]["prediction_band"],"recommended_threshold")
    self.assertFalse(rows[0]["splicing_effect_verified_by_RNA"])
 def test_scores_not_imputed_when_license_absent(self):
    with TemporaryDirectory() as td:
      s,rows=spliceai_local.run(private(),td)
      self.assertEqual(s["variant_scores_resolved"],0)
      self.assertIsNone(rows[0]["maximum_delta_score"])
      self.assertTrue(s["score_absence_is_not_prediction_of_normal_splicing"])

class PKChecks(unittest.TestCase):
 def test_analytical_iv_dimension_and_elimination(self):
    x=quant_pk.analytical_iv_profile(100,5,50,500,0.1,[0,10])
    self.assertAlmostEqual(x[0]["predicted_free_plasma_nM"],400,places=9)
    self.assertAlmostEqual(x[1]["predicted_free_plasma_nM"],400*math.exp(-1),places=7)
 def test_bad_fu_or_dose_rejected(self):
    with self.assertRaises(ValueError):
      quant_pk.analytical_iv_profile(100,5,50,500,1.2,[0])
    with self.assertRaises(ValueError):
      quant_pk.analytical_iv_profile(-1,5,50,500,0.5,[0])
 def test_independent_observed_profile_qualification_screen(self):
    obs=[(0,400),(5,400*math.exp(-.5)),(10,400*math.exp(-1))]
    x=quant_pk.validate_against_observed(obs,obs)
    self.assertAlmostEqual(x["absolute_average_fold_error"],1,places=6)
    self.assertTrue(x["predeclared_twofold_screen_met"])
    bad=quant_pk.validate_against_observed(obs,[(t,v/10) for t,v in obs])
    self.assertFalse(bad["predeclared_twofold_screen_met"])
 def test_qualified_pbpk_requires_project_report_hashes_and_independent_observations(self):
    with TemporaryDirectory() as td:
      base=Path(td)
      f=base/"pred.csv";obs=base/"obs.csv";pr=base/"drug.pksim5";doc=base/"evaluation.pdf"
      putcsv(f,[(0,100),(1,50),(2,25)])
      putcsv(obs,[(0,105),(1,48),(2,28)])
      pr.write_bytes(b"VALIDATED_FAKE_PROJECT")
      doc.write_bytes(b"FAKE_QUALIFICATION_REPORT")
      row={"model_kind":"OSP_PK_SIM_PREVALIDATED","compound":"DRUG",
       "source_reference":"independent study A","independent_validation_reference":"independent study B",
       "analyte":"unbound_parent_drug","population":"adult healthy volunteers",
       "dose_regimen_reference":"published dose protocol","formulation":"solution",
       "route":"intravenous","source_checked":True,"independent_evaluation_cohort":True,
       "study_reports_dose_and_sampling":True,"model_applicability_for_this_population":True,
       "clinical_equations_qualification_documented":True,
       "simulation_export_csv":str(f),"observations_csv":str(obs),
       "pksim_project_file":str(pr),"model_evaluation_report":str(doc),
       "project_sha256":hashlib.sha256(pr.read_bytes()).hexdigest(),
       "report_sha256":hashlib.sha256(doc.read_bytes()).hexdigest()}
      good=quant_pk.assess(row)
      self.assertEqual(good["state"],"SOURCE_QUALIFIED_POPULATION_PK_CONTEXT")
      self.assertFalse(good["model_is_patient_specific"])
      self.assertIsNone(good["predicted_therapeutic_benefit"])
      row["report_sha256"]="0"*64
      bad=quant_pk.assess(row)
      self.assertEqual(bad["state"],"PK_SIM_MODEL_PROVENANCE_HASH_MISMATCH")
 def test_unqualified_or_missing_models_cannot_assign_pk_efficacy(self):
    with TemporaryDirectory() as td:
      summary,models=quant_pk.run(td,manifest_path=None)
      self.assertEqual(summary["models_source_qualified"],0)
      self.assertEqual(models,[])
      self.assertEqual(summary["patients_with_validated_pk_predictions"],0)

class IntegrationChecks(unittest.TestCase):
 def test_dossier_merger_never_converts_assay_to_clinical_efficacy(self):
    v=variant()
    from v2810.expert_genomics import safe_hash
    ref="SHA256:"+safe_hash(v)
    genetics={"genes":{"DEMO1":{"variants":[{"variant_ref":ref,"gaps":[]}],
      "gene_mechanism_confidence":"GENETIC_CAUSALITY_UNCERTAIN"}}}
    vep=[{"chrom":"1","pos":100,"ref":"A","alt":"G",
      "most_severe_consequence":"splice_region_variant","transcript_annotations":[{"gene_symbol":"DEMO1"}]}]
    splice=[{"variant_ref":ref,"state":"MATCHED_SOURCE_SCORE","maximum_delta_score":.72}]
    gnomad=[{"variant_ref":ref,"state":"INDEXED_POPULATION_EVIDENCE",
      "max_valid_global_af":.11,"measurements":[{"cohort":"exomes","status":"LOCALLY_MATCHED","allele_frequency":.11}]}]
    integrated=specialist_evidence.attach(genetics,private(),vep,splice,gnomad)
    merged=genetics["genes"]["DEMO1"]["variants"][0]
    self.assertEqual(integrated["gnomad_variants_matched"],1)
    self.assertEqual(merged["independent_specialist_annotation"]["spliceai_max_delta_score"],.72)
    self.assertFalse(merged["independent_specialist_annotation"]["experimentally_confirmed_variant_impact"])
    self.assertEqual(genetics["genes"]["DEMO1"]["gene_mechanism_confidence"],"GENETIC_CAUSALITY_UNCERTAIN")
 def test_specialist_final_report_says_missed_data_explicitly(self):
    with TemporaryDirectory() as td:
      x=specialist_evidence.make_report(td,
       {"vep":{"tool_state":"UNAVAILABLE_OR_FAILED"},
        "spliceai":{"state":"LICENSE_NOT_ATTESTED"},
        "gnomad":{"state":"NOT_INSTALLED","v4_1_full_coverage":False}},
       {"state":"NO_QUALIFIED_QUANTITATIVE_DRUG_MODELS_INSTALLED","models_source_qualified":0},
       {"patient_variant_count":1})
      content=(Path(td)/"SPECIALIST_EVIDENCE_STATUS.md").read_text()
      self.assertIn("NOT_INSTALLED",content)
      self.assertIn("LICENSE_NOT_ATTESTED",content)
      self.assertIn("not prescriptions",content.lower())


class PKSimNativeBridgeChecks(unittest.TestCase):
 def test_pk_sim_native_csv_time_and_total_plasma(self):
    from v2810 import pksim_bridge
    with TemporaryDirectory() as td:
        f=Path(td)/"exp.csv"
        with f.open("w",encoding="utf-8",newline="") as stream:
            w=csv.writer(stream)
            w.writerow(["IndividualId","Time [min]","Organism|VenousBlood|Drug|Plasma [µmol/l]"])
            w.writerows([["0",0,0.0002],["0",60,0.0001],["1",0,0.0003]])
        cols=pksim_bridge.available_series(f)
        self.assertEqual(len(cols["concentration_columns"]),1)
        row=pksim_bridge.extract_curve(f,cols["concentration_columns"][0])
        self.assertEqual(len(row["observations"]),2)
        self.assertEqual(row["observations"][0]["total_plasma_nM"],0.2)
        self.assertTrue(row["unbound_fraction_source_qualification_required"])
        self.assertNotIn("calculated_free_plasma_nM",row["observations"][0])
 def test_pk_sim_unknown_series_and_bad_fu_are_rejected(self):
    from v2810 import pksim_bridge
    with TemporaryDirectory() as td:
        f=Path(td)/"x.csv"
        f.write_text('IndividualId,Time [min],A [nmol/l]\n0,0,10\n',encoding="utf-8")
        with self.assertRaises(ValueError):
            pksim_bridge.extract_curve(f,"A [nmol/l]",unbound_fraction=1.5)
        with self.assertRaises(ValueError):
            pksim_bridge.extract_curve(f,"A [arbitrary unit]")
 def test_molecule_control_is_not_patient_specific(self):
    from v2810 import pksim_bridge
    with TemporaryDirectory() as td:
        x=pksim_bridge.audit_control(td,Path(td)/"summary.json")
        self.assertEqual(x["simulations_exported"],0)
        self.assertFalse(x["independent_pk_population_qualification_from_just_csv"])
        self.assertTrue(x["control_is_not_a_candidate_for_rare_disease"])

if __name__=="__main__":unittest.main()
