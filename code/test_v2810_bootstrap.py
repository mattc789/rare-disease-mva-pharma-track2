"""Synthetic license/privacy/storage tests for the free auto-download layer."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase,main
from unittest.mock import patch
from v2810 import resource_bootstrap as r
from v2810.vep_local import wsl_path

class BootstrapTests(TestCase):
 def test_chromosomes_derived_from_alleles_not_gene_name(self):
    case={"genes":{"ALPHA":[{"chrom":"chr22"},{"chrom":"chr1"}],
                   "BETA":[{"chrom":"chr5"},{"chrom":"chr22"}]}}
    chrom,unsupported=r.needed_chromosomes(case)
    self.assertEqual(chrom,["1","5","22"])
    self.assertEqual(unsupported,[])
    plan=r.plan(case)
    self.assertEqual(len(plan["download_jobs"]),6)
    self.assertNotIn("ALPHA",json.dumps(plan))
    self.assertNotIn("BETA",json.dumps(plan))
    self.assertNotIn("pos",json.dumps(plan))
    for job in plan["download_jobs"]:
       self.assertEqual(job["assembly"],"GRCh38")
       self.assertTrue(job["source_url"].startswith("https://storage.googleapis.com/gcp-public-data--gnomad/"))
 def test_unsupported_mt_does_not_request_wrong_chr(self):
    case={"genes":{"ANOTHER":[{"chrom":"chrMT"},{"chrom":"chrX"},{"chrom":"Y"}]}}
    p=r.plan(case)
    self.assertEqual(p["needed_chromosomes"],["X","Y"])
    self.assertEqual(p["unsupported_chromosomes"],["chrMT"])
 def test_rejects_unapproved_cohorts(self):
    with self.assertRaises(ValueError):r.plan({"genes":{}},cohorts=("fees_for_predictions",))
 def test_license_profile_rejects_spliceai_and_enrichr(self):
    self.assertEqual(r.CATALOG["Illumina-SpliceAI"]["distribution"],"EXCLUDED_FROM_AUTOMATIC_BOOTSTRAP")
    self.assertEqual(r.CATALOG["Enrichr-LINCS-GEO"]["distribution"],"EXCLUDED_FROM_DEFAULT")
    self.assertEqual(r.CATALOG["VEP"]["license"],"Apache-2.0")
    self.assertEqual(r.CATALOG["Pangolin"]["license"],"GPL-3.0 standalone splice tool")
 def test_wsl_unc_paths_are_local_not_remote(self):
    ctx={"home":"/home/testuser","distribution":"Ubuntu"}
    win=r.windows_unc(ctx,"/home/testuser/rdmva_resources/test")
    self.assertEqual(win,r"\\wsl.localhost\Ubuntu\home\testuser\rdmva_resources\test")
    self.assertEqual(wsl_path(win),"/home/testuser/rdmva_resources/test")
    with self.assertRaises(ValueError):r.windows_unc(ctx,"/etc/passwd")
 def test_no_genetic_allele_urls_and_dry_run_never_downloads(self):
    case={"genes":{"FAKECASE":[{"chrom":"chr21","pos":6789,"ref":"A","alt":"T"}]}}
    with TemporaryDirectory() as td,\
         patch.object(r,"offline_external_bin_wsl",return_value=True),\
         patch.object(r,"wsl_context",return_value={"home":"/home/testuser","distribution":"Ubuntu"}),\
         patch.object(r,"_run_wsl",return_value=""),\
         patch.object(r,"_free_bytes",return_value=300*1024**3),\
         patch.object(r,"_download_meta",return_value=10*1024**3),\
         patch.object(r,"_exists",return_value=False),\
         patch.object(r,"_curl_download") as forbidden:
       manifest,audit=r.prepare(case,td,download=False)
       forbidden.assert_not_called()
       self.assertEqual(audit["completion"],"PLANNED_ONLY")
       self.assertIsNone(manifest)
       j=json.dumps(audit)
       self.assertNotIn("6789",j)
       self.assertNotIn('"ref"',j)
 def test_storage_budget_blocks_partial_jobs_and_reports_reason(self):
    case={"genes":{"FAKE":[{"chrom":"21"}]}}
    with TemporaryDirectory() as td,\
         patch.object(r,"offline_external_bin_wsl",return_value=True),\
         patch.object(r,"wsl_context",return_value={"home":"/home/testuser","distribution":"Ubuntu"}),\
         patch.object(r,"_run_wsl",return_value=""),\
         patch.object(r,"_free_bytes",return_value=18*1024**3),\
         patch.object(r,"_download_meta",return_value=10*1024**3),\
         patch.object(r,"_exists",return_value=False),\
         patch.object(r,"_curl_download") as forbidden:
       manifest,audit=r.prepare(case,td,download=True,max_gib=120)
       forbidden.assert_not_called()
       self.assertEqual(audit["completion"],"PARTIAL_STORAGE_CONSTRAINT")
       self.assertTrue(all(x["status"]=="BLOCKED_STORAGE_BUDGET" for x in audit["records"]))
       self.assertEqual(audit["installed_jobs"],0)
 def test_idempotent_installed_with_valid_sha(self):
    case={"genes":{"FAKE":[{"chrom":"21"}]}}
    def existing(path):return True
    def size(path):return 10000
    with TemporaryDirectory() as td,\
         patch.object(r,"offline_external_bin_wsl",return_value=True),\
         patch.object(r,"wsl_context",return_value={"home":"/home/testuser","distribution":"Ubuntu"}),\
         patch.object(r,"_run_wsl",return_value=""),\
         patch.object(r,"_free_bytes",return_value=100*1024**3),\
         patch.object(r,"_download_meta",return_value=10000),\
         patch.object(r,"_exists",side_effect=existing),\
         patch.object(r,"_file_size",side_effect=size),\
         patch.object(r,"_sha256",return_value="f"*64),\
         patch.object(r,"_curl_download") as forbidden:
       manifest,audit=r.prepare(case,td,download=True)
       forbidden.assert_not_called()
       self.assertEqual(audit["completion"],"READY")
       self.assertEqual(audit["installed_jobs"],2)
       self.assertTrue(Path(manifest).is_file())
       j=json.loads(Path(manifest).read_text())
       self.assertEqual(j["schema"],"V289_GNOMAD_V41_LOCAL_1")
       self.assertEqual(j["shards"][0]["source_sha256"],"f"*64)


class GenomicsCompatibilityTests(TestCase):
 def test_gnomad_official_bgz_file_extension_with_local_index(self):
    from v2810 import gnomad_local
    with TemporaryDirectory() as td:
        folder=Path(td)
        f=folder/"gnomad.exomes.v4.1.sites.chr21.vcf.bgz"
        f.write_bytes(b"fake_BGZF_only_testing_manifest")
        Path(str(f)+".tbi").write_bytes(b"synthetic_index")
        manifest=folder/"m.json"
        manifest.write_text(json.dumps({"schema":"V289_GNOMAD_V41_LOCAL_1",
          "assembly":"GRCh38","release":"v4.1",
          "shards":[{"cohort":"exomes","chrom":"21","file":str(f),
          "official_reference":"https://storage.googleapis.com/gcp-public-data--gnomad/synthetic"}]}))
        records,status=gnomad_local.validate_manifest(manifest)
        self.assertIn(("exomes","21"),records)
        self.assertFalse(status["v4_1_full_coverage"])
        self.assertEqual(status["state"],"LOCAL_PARTIAL")
 def test_free_gpl_pangolin_selected_gene_annotation_subset(self):
    import gzip
    from v2810 import open_splicing
    with TemporaryDirectory() as td:
        source=Path(td)/"public.gtf.gz"
        lines=[
          '1\\tGENCODE\\tgene\\t1\\t50\\t.\\t+\\t.\\tgene_id "g1"; gene_name "GENE_A";\\n',
          '1\\tGENCODE\\ttranscript\\t1\\t50\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1"; gene_name "GENE_A";\\n',
          '1\\tGENCODE\\texon\\t1\\t50\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1"; gene_name "GENE_A";\\n',
          '1\\tGENCODE\\tgene\\t60\\t80\\t.\\t+\\t.\\tgene_id "g2"; gene_name "GENE_B";\\n']
        # literal tab/newline characters rather than escaped source annotations
        with gzip.open(source,"wt",encoding="utf-8") as f:
            f.writelines([line.encode("utf-8").decode("unicode_escape") for line in lines])
        out=Path(td)/"selected.gtf.gz"
        v=open_splicing._write_selected_gene_gtf(source,["GENE_A"],out)
        self.assertEqual(v["public_annotation_records_selected"],3)
        with gzip.open(out,"rt",encoding="utf-8") as f:
            content=f.read()
        self.assertIn('GENE_A',content)
        self.assertNotIn('GENE_B',content)
 def test_default_online_ai_requires_explicit_opt_in(self):
    from v2810 import pipeline
    self.assertEqual(pipeline.run.__defaults__[5],False)


class FreeSplicingScientificBoundaryTests(TestCase):
 def test_open_splicing_scores_are_bound_to_hashed_alleles_only(self):
    from v2810 import open_splicing
    case={"genes":{"DEMO1":[{"gene":"DEMO1","chrom":"21","pos":100,
             "ref":"A","alt":"G","GT":"0/1"},
             {"gene":"DEMO1","chrom":"21","pos":200,
             "ref":"T","alt":"C","GT":"0/1"}]}}
    with TemporaryDirectory() as td:
       out=Path(td)
       file=out/"model_synthetic.vcf"
       file.write_text("##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n"
         "21\\t100\\t.\\tA\\tG\\t.\\tPASS\\tPangolin=DEMO1|+3:0.7|-4:-0.2\\n"
         "21\\t200\\t.\\tT\\tC\\t.\\tPASS\\t.\\n".encode("utf8").decode("unicode_escape"),
         encoding="utf-8")
       rows=open_splicing.parse_private_pangolin_scores(case,file,out)
       self.assertEqual(len(rows),2)
       self.assertEqual(sum(x["state"]=="COMPUTATIONAL_SPLICE_SCORE" for x in rows),1)
       scored=next(x for x in rows if x["state"]=="COMPUTATIONAL_SPLICE_SCORE")
       self.assertAlmostEqual(scored["max_absolute_splice_score_change"],0.7)
       self.assertFalse(scored["RNA_splicing_validated"])
       outjson=(out/"PANGOLIN_PRIVATE_HASHED_SCORES.jsonl").read_text()
       self.assertNotIn('"pos"',outjson)
       self.assertNotIn('"ref"',outjson)

if __name__=="__main__":main()
