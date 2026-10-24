"""Verify GRCh38 reference cache recovery without any patient data."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import tempfile,unittest
from v2810 import resource_bootstrap as r,vep_local

URI="https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

class VEPReferenceCacheRecovery(unittest.TestCase):
 def test_missing_index_uses_existing_fasta_without_any_network(self):
  with tempfile.TemporaryDirectory() as td:
   fa=Path(td)/"Homo_sapiens.GRCh38.dna.primary_assembly.fa"
   with patch.object(r,"_vep_fasta_index_verified",side_effect=[False,True]),\
        patch.object(r,"_cached_large_vep_fasta",return_value=True),\
        patch.object(vep_local,"wsl_path",return_value="/mnt/c/fa"),\
        patch.object(r,"_run_wsl") as run,\
        patch.object(r,"_download_meta",side_effect=AssertionError("unexpected HEAD")),\
        patch.object(r,"_curl_download",side_effect=AssertionError("unexpected download")):
    state=r._restore_vep_reference_from_cache(fa,uri=URI)
    self.assertEqual(state,"REBUILT_INDEX_FOR_CACHED_GRCH38_FASTA")
    run.assert_called_once_with(["samtools","faidx","/mnt/c/fa"],timeout=1800)

 def test_valid_existing_fasta_and_index_never_downloaded(self):
  with tempfile.TemporaryDirectory() as td:
   fa=Path(td)/"cache.fa"
   with patch.object(r,"_vep_fasta_index_verified",return_value=True),\
        patch.object(r,"_run_wsl",side_effect=AssertionError("no commands")) ,\
        patch.object(r,"_download_meta",side_effect=AssertionError("no web")):
    self.assertEqual(r._restore_vep_reference_from_cache(fa,uri=URI),
                     "VERIFIED_EXISTING_FASTA_AND_INDEX")

 def test_complete_cached_gzip_reused_when_no_extracted_fasta(self):
  with tempfile.TemporaryDirectory() as td:
   fa=Path(td)/"cache.fa"
   archive=Path(str(fa)+".gz")
   archive.write_bytes(b"synthetic test only")
   def fake_wsl(p):return "/mnt/z/"+Path(p).name
   with patch.object(r,"_vep_fasta_index_verified",side_effect=[False,True,True]),\
        patch.object(r,"_cached_large_vep_fasta",return_value=False),\
        patch.object(r,"_cached_official_fasta_archive",return_value=True),\
        patch.object(r,"_download_meta",return_value=881964081),\
        patch.object(r.shutil,"disk_usage",return_value=SimpleNamespace(free=50*1024**3)),\
        patch.object(vep_local,"wsl_path",side_effect=fake_wsl),\
        patch.object(r,"_run_wsl") as run,\
        patch.object(r,"_curl_download",side_effect=AssertionError("unnecessary download")):
    state=r._restore_vep_reference_from_cache(fa,uri=URI)
    self.assertEqual(state,"RESTORED_FROM_VERIFIED_CACHED_ARCHIVE")
    cmds=[call.args[0] for call in run.call_args_list]
    self.assertTrue(any(c[0:2]==["gzip","-t"] for c in cmds))
    self.assertTrue(any(c[0:2]==["samtools","faidx"] for c in cmds))

 def test_invalid_small_or_bad_index_fails_without_false_ready(self):
  with tempfile.TemporaryDirectory() as td:
   fa=Path(td)/"cache.fa"
   fa.write_bytes(b"TRUNCATED")
   Path(str(fa)+".fai").write_text(
     "1\t248956422\t56\t60\t61\n7\t159345973\t100\t60\t61\n",encoding="ascii")
   self.assertFalse(r._vep_fasta_index_verified(fa))

 def test_gencode_gtf_is_coordinate_sorted_before_tabix(self):
  # GENCODE contains transcripts with decreasing starts in source order.
  # Merely re-bgzipping it makes tabix abort on the first inverted position.
  source=Path(r.__file__).read_text(encoding="utf-8")
  self.assertIn("LC_ALL=C sort --stable",source)
  self.assertIn(" -k1,1 -k4,4n",source)
  self.assertIn('["bash","-o","pipefail","-c",cmd]',source)
  self.assertIn('["tabix","-f","-p","gff",wsl_path(stage)]',source)
  self.assertIn('os.replace(stage_index,Path(str(gtf)+".tbi"))',source)

 def test_vep_gtf_fasta_uses_local_mode_without_cache_or_database(self):
  with tempfile.TemporaryDirectory() as td:
   infile=Path(td)/"synthetic.vcf"
   outfile=Path(td)/"synthetic.json"
   infile.write_text("##fileformat=VCFv4.2\n",encoding="ascii")
   outfile.write_text("{}\n",encoding="ascii")
   paths={"vep":"/home/test/vep","gtf":"/mnt/t/synthetic.gtf.gz","fasta":"/mnt/t/synthetic.fa"}
   with patch.object(vep_local,"default_paths",return_value=paths),\
        patch.object(vep_local,"wsl_path",return_value="/mnt/t/synthetic"),\
        patch.object(vep_local.subprocess,"run",return_value=SimpleNamespace(
          returncode=0,stderr="",stdout="")) as invoke:
    vep_local._run_vep(infile,outfile,vep=paths["vep"],
                       gtf=paths["gtf"],fasta=paths["fasta"])
    args=invoke.call_args.args[0]
    self.assertIn("--gtf",args)
    self.assertIn("--fasta",args)
    self.assertNotIn("--offline",args)
    self.assertNotIn("--cache",args)
    self.assertNotIn("--database",args)
    vep_local._run_vep(infile,outfile,vep=paths["vep"],
                       gtf=paths["gtf"],fasta=paths["fasta"],cache_dir="/cache")
    args=invoke.call_args.args[0]
    self.assertIn("--offline",args)
    self.assertIn("--cache",args)
    self.assertNotIn("--database",args)

 def test_vep_executable_dependency_preflight(self):
  source=Path(r.__file__).read_text(encoding="utf-8")
  for package in ("libtry-tiny-perl","liburi-perl","libwww-perl"):
   self.assertIn(package,source)
  self.assertIn('vep,"--help"',source)

 def test_source_pinned_to_prior_real_grch38_reference_hash(self):
  source=Path(r.__file__).read_text(encoding="utf-8")
  self.assertIn("1e74081a49ceb9739cc14c812fbb8b3db978eb80ba8e5350beb80d8ad8dfef3b",source)
  self.assertIn("V31_ENSEMBL_GRCH38_FASTA_PROVENANCE_SHA256_MISMATCH",source)

if __name__=="__main__":unittest.main()
