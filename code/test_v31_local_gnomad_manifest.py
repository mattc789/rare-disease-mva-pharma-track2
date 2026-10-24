"""Regression for supplied local gnomAD: no network, no patient fixtures."""
from __future__ import annotations
import hashlib, json, tempfile, unittest
from pathlib import Path
from v30.native_pipeline import audit_supplied_gnomad_manifest

def digest(b): return hashlib.sha256(b).hexdigest()

class SuppliedManifestAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
    def make(self,chroms=("7",),cohorts=("exomes","genomes"),bad=None):
        items=[]
        for co in cohorts:
            for chrom in chroms:
                path=self.root/co/f"gnomad.{co}.v4.1.sites.chr{chrom}.vcf.bgz"
                path.parent.mkdir(parents=True,exist_ok=True)
                source=b"synthetic_public_variant_frequency_only\n"
                index=b"synthetic_index_fixture"
                path.write_bytes(source)
                Path(str(path)+".tbi").write_bytes(index)
                items.append({
                    "cohort":co,"chrom":chrom,"file":str(path),
                    "official_reference":"https://example.invalid/public-synthetic-fixture",
                    "source_sha256":digest(source),
                    "index_sha256":digest(index),
                })
        if bad:bad(items)
        manifest=self.root/"manifest.json"
        manifest.write_text(json.dumps({
            "schema":"V289_GNOMAD_V41_LOCAL_1",
            "assembly":"GRCh38","release":"v4.1","shards":items,
        }),encoding="utf-8")
        return manifest
    def private(self,chroms=("7",)):
        return {"genes":{"SYNTHETIC":[
            {"chrom":c,"pos":101,"ref":"A","alt":"G"} for c in chroms]}}
    def test_uses_official_hash_pinned_local_cohorts_and_never_calls_download(self):
        m=self.make()
        a=audit_supplied_gnomad_manifest(self.private(),m)
        self.assertEqual(a["completion"],"READY")
        self.assertEqual(a["mode"],"EXPLICIT_LOCAL_MANIFEST_NO_DOWNLOAD")
        self.assertEqual(a["required_jobs"],2)
        self.assertEqual(a["installed_jobs"],2)
        self.assertEqual(a["downloaded_gib_this_run"],0)
        self.assertFalse(a["v4_1_full_coverage"])
        self.assertFalse(a["patient_variant_coordinates_sent_to_public_hosts"])
    def test_multi_chromosome_coverage_has_precise_scope(self):
        m=self.make(chroms=("7","15","22"))
        a=audit_supplied_gnomad_manifest(self.private(("7","15","22")),m)
        self.assertEqual((a["required_jobs"],a["installed_jobs"]),(6,6))
        self.assertEqual(a["selected_case_chromosome_count"],3)
    def test_missing_genomes_fails_no_network(self):
        m=self.make(cohorts=("exomes",))
        with self.assertRaisesRegex(RuntimeError,"MISSING_CASE_SHARDS"):
            audit_supplied_gnomad_manifest(self.private(),m)
    def test_missing_required_chromosome_fails_no_network(self):
        m=self.make()
        with self.assertRaisesRegex(RuntimeError,"MISSING_CASE_SHARDS"):
            audit_supplied_gnomad_manifest(self.private(("7","15")),m)
    def test_unaudited_source_checksum_fails_closed(self):
        m=self.make(bad=lambda rows: rows[0].pop("source_sha256"))
        with self.assertRaisesRegex(RuntimeError,"NEEDS_SOURCE_AND_INDEX_SHA256_PINS"):
            audit_supplied_gnomad_manifest(self.private(),m)
    def test_unaudited_index_checksum_fails_closed(self):
        m=self.make(bad=lambda rows: rows[0].pop("index_sha256"))
        with self.assertRaisesRegex(RuntimeError,"NEEDS_SOURCE_AND_INDEX_SHA256_PINS"):
            audit_supplied_gnomad_manifest(self.private(),m)
    def test_changed_existing_source_rejected(self):
        m=self.make()
        path=Path(json.loads(m.read_text())["shards"][0]["file"])
        path.write_bytes(b"modified")
        with self.assertRaisesRegex(RuntimeError,"SOURCE_SHA256_MISMATCH"):
            audit_supplied_gnomad_manifest(self.private(),m)
    def test_changed_existing_index_rejected(self):
        m=self.make()
        path=json.loads(m.read_text())["shards"][0]["file"]+".tbi"
        Path(path).write_bytes(b"changed index")
        with self.assertRaisesRegex(RuntimeError,"INDEX_SHA256_MISMATCH"):
            audit_supplied_gnomad_manifest(self.private(),m)
    def test_invalid_manifest_missing_index_fails(self):
        m=self.make()
        path=json.loads(m.read_text())["shards"][0]["file"]+".tbi"
        Path(path).unlink()
        with self.assertRaisesRegex(RuntimeError,"MISSING_CASE_SHARDS"):
            audit_supplied_gnomad_manifest(self.private(),m)
    def test_unsupported_chromosome_fails(self):
        m=self.make()
        with self.assertRaisesRegex(RuntimeError,"CASE_CHROMOSOMES_INVALID"):
            audit_supplied_gnomad_manifest(self.private(("MT",)),m)
    def test_empty_case_fails(self):
        m=self.make()
        with self.assertRaisesRegex(RuntimeError,"CASE_CHROMOSOMES_INVALID"):
            audit_supplied_gnomad_manifest(self.private(()),m)

if __name__=="__main__":unittest.main()
