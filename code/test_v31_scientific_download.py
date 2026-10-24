"""Offline tests for accelerated PUBLIC data only (no patient input or network)."""
from __future__ import annotations
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from unittest.mock import patch
import hashlib
import os
import tempfile
import threading
import time
import unittest

from v2810 import scientific_download_v31 as sd

DATA = bytes(range(256)) * 4096  # 1 MiB of synthetic public-like bytes
DIGEST = hashlib.sha256(DATA).hexdigest()


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    def log_message(self, *_args):
        pass

    def do_GET(self):
        server = self.server
        with server.counter_lock:
            server.gets += 1
        original_range = self.headers.get("Range")
        ranged = self.path != "/single" and original_range is not None
        start, end = 0, len(DATA) - 1
        if ranged:
            import re
            match = re.fullmatch(r"bytes=(\d+)-(\d+)", original_range)
            if not match:
                self.send_error(416)
                return
            start, end = map(int, match.groups())
            if start > end or end >= len(DATA):
                self.send_error(416)
                return
        payload = DATA[start:end+1] if ranged else DATA
        self.send_response(206 if ranged else 200)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("ETag", '"known-synthetic-source-1"')
        if ranged:
            # Probe is valid. An erroneous later range must be rejected.
            claimed = start + 1 if self.path == "/bad" and start else start
            self.send_header("Content-Range", f"bytes {claimed}-{end}/{len(DATA)}")
        self.end_headers()
        if ranged and start:
            with server.counter_lock:
                server.active += 1
                server.max_active = max(server.max_active, server.active)
            time.sleep(.012)
        try:
            self.wfile.write(payload)
            self.wfile.flush()
        except (ConnectionResetError, BrokenPipeError):
            pass
        finally:
            if ranged and start:
                with server.counter_lock:
                    server.active -= 1


class AcceleratedPublicDatasetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        cls.server.counter_lock = threading.Lock()
        cls.server.gets = cls.server.active = cls.server.max_active = 0
        cls.worker = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.worker.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.worker.join(timeout=2)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "synthetic.vcf.bgz"
        self.url = f"http://127.0.0.1:{self.server.server_port}"
        with self.server.counter_lock:
            self.server.gets = 0
            self.server.max_active = 0

    def get(self, path, **kw):
        return sd.download_public_data(
            self.url + path, self.out, expected_bytes=len(DATA),
            max_bytes=len(DATA), expected_sha256=DIGEST,
            parallel_threshold=32_768, segment_size=65_536,
            connections=4, allow_test_http=True, **kw
        )

    def test_parallel_range_download_is_bit_identical_and_reused(self):
        result = self.get("/range")
        self.assertEqual(result["mode"], "segmented_https")
        self.assertEqual(self.out.read_bytes(), DATA)
        self.assertGreater(self.server.max_active, 1)
        before = self.server.gets
        reused = self.get("/range")
        self.assertEqual(reused["mode"], "already_present")
        self.assertEqual(before, self.server.gets)

    def test_partial_checkpoints_are_verified_and_resumed(self):
        original = sd._one_range
        def fail_last(url, path, start, end, total, validator):
            if end == len(DATA)-1:
                time.sleep(.15)
                raise sd.TransferError("SYNTHETIC_INTERRUPTION")
            return original(url, path, start, end, total, validator)
        with patch.object(sd, "_one_range", side_effect=fail_last):
            with self.assertRaisesRegex(sd.TransferError, "SYNTHETIC_INTERRUPTION"):
                self.get("/range")
        self.assertFalse(self.out.exists())
        manifest = Path(str(self.out)+".part.rdmva.json")
        self.assertTrue(manifest.is_file())
        self.assertTrue(sd.json.loads(manifest.read_text())["done"])
        resumed = self.get("/range")
        self.assertEqual(resumed["mode"], "segmented_https")
        self.assertEqual(self.out.read_bytes(), DATA)
        self.assertFalse(manifest.exists())

    def test_no_range_host_streams_single_file(self):
        res = self.get("/single")
        self.assertEqual(res["mode"], "single_https")
        self.assertEqual(self.out.read_bytes(), DATA)

    def test_broken_content_range_never_promotes_output(self):
        with self.assertRaisesRegex(sd.TransferError, "RANGE_IGNORED"):
            self.get("/bad")
        self.assertFalse(self.out.exists())

    def test_wrong_sha_never_promotes_single_stream(self):
        with self.assertRaisesRegex(sd.TransferError, "SHA256_MISMATCH"):
            sd.download_public_data(self.url+"/single",self.out,
                expected_bytes=len(DATA),max_bytes=len(DATA),
                expected_sha256="0"*64,allow_test_http=True)
        self.assertFalse(self.out.exists())

    def test_caps_source_identity_and_https_guard(self):
        with self.assertRaisesRegex(sd.TransferError,"HTTPS_SOURCE_REQUIRED"):
            sd.download_public_data("http://example.net/file",self.out,max_bytes=100)
        with self.assertRaisesRegex(sd.TransferError,"EXCEEDS_ALLOWED_CAP"):
            sd.download_public_data(self.url+"/range",self.out,
                max_bytes=len(DATA)//2,allow_test_http=True)

    def test_all_scientific_source_bootstraps_use_shared_transport(self):
        root = Path(__file__).resolve().parents[1]
        paths = ("bootstrap_current_public_v31.py","bootstrap_v28.py",
                 "bootstrap_regulatory_v283.py",
                 "code/v2810/resource_bootstrap.py",
                 "code/v2810/chembl_sqlite_download_v31.py")
        for name in paths:
            content = (root/name).read_text(encoding="utf-8")
            self.assertIn("download_public_data",content,name)


if __name__ == "__main__":
    unittest.main()
