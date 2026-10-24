from __future__ import annotations
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v2810.download_resilience_v31 import (
    DownloadWarnings, resilient_download, RETRIES, RETRY_DELAY_SECONDS,
)

class PublicDataRecoveryTests(unittest.TestCase):
    def test_exactly_three_retries_then_alternative(self):
        with tempfile.TemporaryDirectory() as d:
            ledger = DownloadWarnings(d)
            calls=[]
            waited=[]
            def failing():
                calls.append("official")
                raise RuntimeError("HTTP 503 Service Unavailable")
            def alternative():
                calls.append("alternative")
                return {"verified": True, "source":"government_zip"}
            output=io.StringIO()
            with redirect_stdout(output):
                result=resilient_download("CA_APPROVED","https://health-products.canada.ca/api/drug/",
                    failing,fallback=alternative,fallback_name="Canada_official_ZIP",
                    warnings=ledger,sleeper=waited.append)
            self.assertEqual(result["source"],"government_zip")
            self.assertEqual(calls,["official"]*4+["alternative"])
            self.assertEqual(waited,[120,120,120])
            stored=json.loads((Path(d)/"INSTALL_WARNINGS.json").read_text())
            self.assertEqual(stored["skipped_sources"],[])
            self.assertEqual(len(stored["recovered_sources"]),1)
            self.assertIn("HTTP 503",output.getvalue())

    def test_skip_proceeds_to_next_download_and_menu_reminder_persists(self):
        with tempfile.TemporaryDirectory() as d:
            ledger=DownloadWarnings(d)
            waited=[]
            out=io.StringIO()
            def failure():raise OSError("upstream timeout")
            with redirect_stdout(out):
                result=resilient_download("UNAVAILABLE_REFERENCE",
                    "https://ftp.ebi.ac.uk/missing",failure,
                    fallback=lambda:failure(),fallback_name="alternate_official",
                    warnings=ledger,sleeper=waited.append)
                next_download=resilient_download("NEXT_REFERENCE",
                    "https://www.ema.europa.eu/official",
                    lambda:{"verified":True},warnings=ledger,sleeper=waited.append)
                ledger.print_summary()
            self.assertIsNone(result)
            self.assertTrue(next_download["verified"])
            self.assertEqual(waited,[120,120,120])
            self.assertIn("INSTALLATION",out.getvalue().upper())
            self.assertIn("UNAVAILABLE_REFERENCE",out.getvalue())
            self.assertIn("MISSING",out.getvalue())
            persisted=DownloadWarnings(d)
            self.assertEqual(len(persisted.data["skipped_sources"]),1)
            with redirect_stdout(io.StringIO()):
                resilient_download("UNAVAILABLE_REFERENCE",
                    "https://ftp.ebi.ac.uk/missing",
                    lambda:{"verified":True},warnings=persisted,
                    sleeper=lambda _:self.fail("must not retry if healthy"))
            self.assertFalse(DownloadWarnings(d).data["skipped_sources"])

    def test_transient_success_after_one_120_second_wait(self):
        calls=[]
        sleeps=[]
        def first_fails():
            calls.append(1)
            if len(calls)==1:raise TimeoutError("offline")
            return True
        with redirect_stdout(io.StringIO()):
            self.assertTrue(resilient_download("TRANSIENT","https://valid.com",
                first_fails,sleeper=sleeps.append))
        self.assertEqual(len(calls),2)
        self.assertEqual(sleeps,[120])

    def test_policy_constants_and_menu_integrity(self):
        self.assertEqual(RETRIES,3)
        self.assertEqual(RETRY_DELAY_SECONDS,120)
        root=Path(__file__).resolve().parents[1]
        ps=(root/"JURY_START_V31.ps1").read_text(encoding="utf-8")
        for option in range(8):
            self.assertIn("'"+str(option)+"' {",ps)
        self.assertIn("Show-DownloadWarnings",ps)
        self.assertNotIn("Type 6 to review",ps)

if __name__=="__main__":
    unittest.main()
