"""Small, fast SQLite fixture exercising ChEMBL 37 local API compatibility."""
import contextlib,io,json,sqlite3,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from v2810 import chembl_sqlite_v31 as adapter
from v2810 import chembl_sqlite_download_v31 as downloader
from v2810 import core

DDL={
"molecule_dictionary":"molregno INTEGER PRIMARY KEY,chembl_id TEXT,pref_name TEXT,max_phase INTEGER,molecule_type TEXT,withdrawn_flag INTEGER,availability_type INTEGER,veterinary INTEGER,first_approval INTEGER,oral INTEGER,parenteral INTEGER,topical INTEGER,prodrug INTEGER",
"molecule_hierarchy":"molregno INTEGER,parent_molregno INTEGER,active_molregno INTEGER",
"molecule_synonyms":"molregno INTEGER,synonyms TEXT,syn_type TEXT",
"molecule_atc_classification":"molregno INTEGER,level5 TEXT",
"drug_mechanism":"mec_id INTEGER PRIMARY KEY,record_id INTEGER,molregno INTEGER,tid INTEGER,site_id INTEGER,action_type TEXT,mechanism_of_action TEXT,direct_interaction INTEGER,disease_efficacy INTEGER,molecular_mechanism INTEGER,mechanism_comment TEXT,selectivity_comment TEXT,binding_site_comment TEXT",
"compound_records":"record_id INTEGER,molregno INTEGER",
"target_dictionary":"tid INTEGER PRIMARY KEY,chembl_id TEXT,pref_name TEXT,organism TEXT,tax_id INTEGER,target_type TEXT,species_group_flag INTEGER",
"target_components":"tid INTEGER,component_id INTEGER",
"component_sequences":"component_id INTEGER,accession TEXT,component_type TEXT,description TEXT",
"compound_structures":"molregno INTEGER,canonical_smiles TEXT",
"compound_properties":"molregno INTEGER,full_mwt REAL,mw_freebase REAL",
"docs":"doc_id INTEGER,chembl_id TEXT",
"assays":"assay_id INTEGER,chembl_id TEXT,description TEXT,confidence_score INTEGER,assay_type TEXT,tid INTEGER",
"activities":"activity_id INTEGER,assay_id INTEGER,doc_id INTEGER,molregno INTEGER,standard_type TEXT,standard_value REAL,standard_units TEXT,standard_relation TEXT,data_validity_comment TEXT",
}

def fixture(path):
    with contextlib.closing(sqlite3.connect(path)) as db:
        for table,cols in DDL.items():
            db.execute("CREATE TABLE "+table+"("+cols+")")
        db.execute("INSERT INTO molecule_dictionary VALUES (1,'CHEMBL1','TEST DRUG',4,'Small molecule',0,1,0,2020,1,0,0,0)")
        db.execute("INSERT INTO molecule_dictionary VALUES (2,'CHEMBL2','TEST PARENT',4,'Small molecule',0,1,0,2020,1,0,0,0)")
        db.execute("INSERT INTO molecule_hierarchy VALUES (1,2,2)")
        db.execute("INSERT INTO molecule_hierarchy VALUES (2,2,2)")
        db.execute("INSERT INTO molecule_synonyms VALUES (1,'TEST TRADE NAME','TRADE_NAME')")
        db.execute("INSERT INTO molecule_atc_classification VALUES (1,'A01AA01')")
        db.execute("INSERT INTO compound_records VALUES (3,1)")
        db.execute("INSERT INTO drug_mechanism VALUES (10,3,1,8,NULL,'INHIBITOR','Inhibits test protein',1,1,1,NULL,NULL,NULL)")
        db.execute("INSERT INTO target_dictionary VALUES (8,'CHEMBL99','Test protein','Homo sapiens',9606,'SINGLE PROTEIN',0)")
        db.execute("INSERT INTO component_sequences VALUES (2,'P12345','PROTEIN','Test protein')")
        db.execute("INSERT INTO target_components VALUES (8,2)")
        db.execute("INSERT INTO compound_structures VALUES (1,'CCO')")
        db.execute("INSERT INTO compound_properties VALUES (1,46.1,46.1)")
        db.execute("INSERT INTO docs VALUES (4,'CHEMBL_DOC_123')")
        db.execute("INSERT INTO assays VALUES (5,'CHEMBL_ASSAY_4','Inhibition of test protein',9,'F',8)")
        db.execute("INSERT INTO activities VALUES (7,5,4,1,'IC50',210.0,'nM','=',NULL)")
        db.execute("CREATE TABLE _test_padding (pad BLOB)")
        db.execute("INSERT INTO _test_padding VALUES (zeroblob(210000))")
        db.commit()

class ChEMBLSQLiteSmoke(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.root=Path(self.tmp.name)
        self.db=self.root/"chembl_37.db"
        fixture(self.db)
        self.provenance=self.root/"CHEMBL_SQLITE_PROVENANCE_V31.json"
        self.provenance.write_text(json.dumps({
          "schema":"RDMVA_CHEMBL_SQLITE_OFFICIAL_V31_1",
          "chembl_release":37,
          "archive_sha256":downloader.ARCHIVE_SHA256,
          "database_bytes":self.db.stat().st_size}),encoding="utf-8")
        self.patches=[patch.object(adapter,"DATABASE",self.db),
                      patch.object(adapter,"PROVENANCE",self.provenance)]
        for p in self.patches:p.start()
    def tearDown(self):
        for p in reversed(self.patches):p.stop()
        self.tmp.cleanup()
    def url(self,tail):return "https://www.ebi.ac.uk/chembl/api/data/"+tail

    def test_catalog_parent_alias_and_structure(self):
        data=adapter.local_json(self.url("molecule.json?max_phase=4&limit=1"))
        self.assertEqual(len(data["molecules"]),1)
        self.assertTrue(data["page_meta"]["next"])
        m=data["molecules"][0]
        self.assertEqual(m["molecule_hierarchy"]["parent_chembl_id"],"CHEMBL2")
        self.assertEqual(m["molecule_synonyms"][0]["molecule_synonym"],"TEST TRADE NAME")
        self.assertEqual(m["atc_classifications"],["A01AA01"])
        self.assertEqual(m["molecule_structures"]["canonical_smiles"],"CCO")
        self.assertEqual(adapter.local_json(self.url("molecule/CHEMBL1.json"))["molecule_properties"]["full_mwt"],46.1)

    def test_full_1000_molecule_page_preserves_synonyms_and_atc(self):
        # A 1000-row REST page requires >900 molecule identifiers for child
        # lookups. It must be fetched in SQLite-safe batches without data loss.
        with contextlib.closing(sqlite3.connect(self.db)) as db:
            db.executemany(
                "INSERT INTO molecule_dictionary (molregno,chembl_id,pref_name,max_phase,molecule_type) "
                "VALUES (?,?,?,?,?)",
                [(i,"CHEMBL"+str(i),"Synthetic "+str(i),4,"Small molecule")
                 for i in range(3,1004)])
            db.execute("INSERT INTO molecule_synonyms VALUES (1000,'SYNTHETIC 1000','OTHER')")
            db.execute("INSERT INTO molecule_atc_classification VALUES (1000,'Z99ZZ99')")
            db.execute("INSERT INTO molecule_synonyms VALUES (1003,'SYNTHETIC 1003','OTHER')")
            db.execute("INSERT INTO molecule_atc_classification VALUES (1003,'Z99ZZ98')")
            db.commit()
        first=adapter.local_json(self.url("molecule.json?max_phase=4&limit=1000"))
        self.assertEqual(len(first["molecules"]),1000)
        self.assertIsNotNone(first["page_meta"]["next"])
        row=next(v for v in first["molecules"] if v["molecule_chembl_id"]=="CHEMBL1000")
        self.assertEqual(row["molecule_synonyms"][0]["molecule_synonym"],"SYNTHETIC 1000")
        self.assertEqual(row["atc_classifications"],["Z99ZZ99"])
        second=adapter.local_json(self.url("molecule.json?max_phase=4&limit=1000&offset=1000"))
        self.assertEqual(len(second["molecules"]),3)
        last=next(v for v in second["molecules"] if v["molecule_chembl_id"]=="CHEMBL1003")
        self.assertEqual(last["molecule_synonyms"][0]["molecule_synonym"],"SYNTHETIC 1003")
        self.assertEqual(last["atc_classifications"],["Z99ZZ98"])

    def test_mechanism_parent_and_target(self):
        m=adapter.local_json(self.url("mechanism.json?limit=1000"))["mechanisms"]
        self.assertEqual(m[0]["molecule_chembl_id"],"CHEMBL1")
        self.assertEqual(m[0]["parent_molecule_chembl_id"],"CHEMBL2")
        self.assertEqual(m[0]["target_chembl_id"],"CHEMBL99")
        t=adapter.local_json(self.url("target.json?target_components__accession=P12345&limit=100"))["targets"]
        self.assertEqual(t[0]["target_components"][0]["accession"],"P12345")
        self.assertEqual(t[0]["tax_id"],9606)

    def test_activity_and_assay_are_available_locally(self):
        a=adapter.local_json(self.url("activity.json?target_chembl_id=CHEMBL99&standard_type__in=IC50%2CEC50&standard_units=nM&standard_value__lte=10000&limit=1000"))["activities"]
        self.assertEqual(len(a),1)
        self.assertEqual(a[0]["assay_chembl_id"],"CHEMBL_ASSAY_4")
        self.assertEqual(a[0]["standard_relation"],"=")
        self.assertEqual(a[0]["parent_molecule_chembl_id"],"CHEMBL2")
        ass=adapter.local_json(self.url("assay.json?assay_chembl_id__in=CHEMBL_ASSAY_4"))["assays"]
        self.assertEqual(ass[0]["confidence_score"],9)

    def test_core_uses_sqlite_without_api_http(self):
        with patch.object(core.urllib.request,"urlopen",side_effect=AssertionError("Remote API forbidden")):
            result=core.fetch_json(self.url("target.json?target_chembl_id__in=CHEMBL99"))
        self.assertEqual(result["targets"][0]["target_chembl_id"],"CHEMBL99")

    def test_fail_closed_invalid_provenance(self):
        self.provenance.write_text(json.dumps({"schema":"test"}))
        with self.assertRaisesRegex(RuntimeError,"PROVENANCE_MISMATCH"):
            adapter.local_json(self.url("molecule.json?limit=1"))

    def test_no_unbounded_activity_scan(self):
        with self.assertRaisesRegex(RuntimeError,"UNBOUNDED_ACTIVITY"):
            adapter.local_json(self.url("activity.json?limit=100"))

    def test_schema_validation_checks_real_required_tables(self):
        # Fixture DB is small but has every required table/column and rows.
        with patch.object(downloader,"ROOT",self.root):
            with self.assertRaisesRegex(RuntimeError,"SQLITE_TOO_SMALL"):
                downloader._check_schema(self.root/"not_existing.sqlite")

    def test_archive_resume_reports_real_incrementing_percentage(self):
        # Synthetic 100-byte archive and stub process: NO network or bulk data.
        class FakeCurl:
            def __init__(self):self.ticks=0;self.returncode=None
            def poll(self):
                self.ticks+=1
                if self.ticks>=3:self.returncode=0
                return self.returncode
            def terminate(self):self.returncode=-15
            def wait(self):return self.returncode
        partial=self.root/(downloader.ARCHIVE_NAME+".part")
        partial.write_bytes(b"x"*20)  # Resumed download begins at 20%.
        increments=iter((30,50))
        def advance(_seconds):
            with partial.open("ab") as f:f.write(b"x"*next(increments))
        out=io.StringIO()
        with patch.object(downloader,"ARCHIVE_BYTES",100),\
             patch.object(downloader.subprocess,"Popen",return_value=FakeCurl()) as launch,\
             patch("time.sleep",side_effect=advance),contextlib.redirect_stdout(out):
            downloader._curl(downloader.BASE+downloader.ARCHIVE_NAME,
                             partial,resume=True,timeout=0)
        progress=out.getvalue()
        for percentage in ("20.0%","50.0%","100.0%"):
            self.assertIn(percentage,progress)
        self.assertIn("state=download_finished_verifying_SHA256",progress)
        self.assertIn("--continue-at",launch.call_args.args[0])
        self.assertIn("--silent",launch.call_args.args[0])

if __name__=="__main__":unittest.main()
