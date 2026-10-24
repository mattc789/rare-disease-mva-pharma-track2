"""Accelerated, integrity-preserving public SCIENTIFIC data transfers (V3.1).

Standard-library only: nothing new to install on a judge's Windows PC.
HTTPS sources are selected/pinned by each scientific bootstrapper. No patient
genome, case metadata, credentials or API request is handled by this module.
Range-capable, immutable responses are fetched in bounded parallel segments;
other servers use a safe single HTTPS stream. Exact caller SHA-256, size, ZIP/
VCF validation remain mandatory. Checkpoint files are excluded from releases.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, build_opener, HTTPRedirectHandler
import hashlib
import json
import os
import re
import threading
import time

MIB = 1024 * 1024
DEFAULT_THRESHOLD = 32 * MIB
DEFAULT_SEGMENT = 64 * MIB
RETRYABLE = {408, 425, 429, 500, 502, 503, 504}
UA = "RareDiseaseMVAPharma-V31-PublicDatasets/1.0"
SCHEMA = "RDMVA_SCIENTIFIC_DATA_SEGMENTS_V31_1"


class TransferError(RuntimeError):
    pass


class _SecureRedirects(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urlsplit(newurl).scheme != "https" and not (
            _LOOPBACK_TEST_MODE and urlsplit(newurl).scheme == "http"
            and urlsplit(newurl).hostname in ("127.0.0.1", "localhost")
        ):
            raise TransferError("RDMVA_DATA_REDIRECT_TO_NON_TLS_REFUSED")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


_LOOPBACK_TEST_MODE = False
_lock = threading.RLock()
_OPENER = build_opener(_SecureRedirects())


def _allowed_url(url, allow_test_http=False):
    parsed = urlsplit(url)
    if parsed.scheme == "https" and parsed.hostname and not parsed.username:
        return
    if (allow_test_http and parsed.scheme == "http" and
            parsed.hostname in ("localhost", "127.0.0.1")):
        return
    raise TransferError("RDMVA_DATA_OFFICIAL_HTTPS_SOURCE_REQUIRED")


def _open(url, headers, timeout=55):
    request = Request(url, headers={"User-Agent": UA, "Accept-Encoding": "identity", **headers})
    return _OPENER.open(request, timeout=timeout)


def _digest(path, begin=0, length=None):
    h = hashlib.sha256()
    with Path(path).open("rb") as src:
        src.seek(begin)
        remaining = length
        while True:
            block = src.read(MIB if remaining is None else min(MIB, remaining))
            if not block:
                break
            h.update(block)
            if remaining is not None:
                remaining -= len(block)
                if remaining == 0:
                    break
    if remaining not in (None, 0):
        raise TransferError("RDMVA_DATA_CHECKPOINT_TRUNCATED")
    return h.hexdigest()


def _atomic_json(path, obj):
    temp = Path(str(path) + ".tmp")
    temp.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temp, path)


def _probe(url, *, expected_bytes=None):
    """One-byte Range GET: verifies server support; never reads the whole object."""
    try:
        with _open(url, {"Range": "bytes=0-0"}, timeout=45) as response:
            status = response.status
            content_range = response.headers.get("Content-Range", "")
            match = re.fullmatch(r"bytes 0-0/(\d+)", content_range.strip())
            if status == 206 and match and response.headers.get("Content-Encoding", "identity").lower() == "identity":
                if len(response.read(2)) != 1:
                    raise TransferError("RDMVA_DATA_RANGE_PROBE_INVALID_BODY")
                size = int(match.group(1))
                if expected_bytes is not None and size != expected_bytes:
                    raise TransferError("RDMVA_DATA_OFFICIAL_LENGTH_CHANGED")
                etag = response.headers.get("ETag", "").strip()
                modified = response.headers.get("Last-Modified", "").strip()
                # Weak ETags are not safe If-Range validators.
                validator = etag if etag.startswith('"') and etag.endswith('"') else modified
                return size, validator, bool(validator)
            # Range-ignoring hosts are handled by the single-stream fallback.
            if status == 200:
                size = response.headers.get("Content-Length")
                if expected_bytes is not None and size and int(size) != expected_bytes:
                    raise TransferError("RDMVA_DATA_OFFICIAL_LENGTH_CHANGED")
                return int(size) if size else expected_bytes, "", False
            raise TransferError("RDMVA_DATA_INVALID_RANGE_PROBE_RESPONSE")
    except HTTPError as exc:
        if exc.code in (400, 405, 416, 501):
            return expected_bytes, "", False
        raise


def _one_range(url, path, start, end, total, validator):
    headers = {"Range": f"bytes={start}-{end}", "If-Range": validator}
    with _open(url, headers, timeout=120) as response:
        header = response.headers.get("Content-Range", "").strip()
        if response.status != 206 or header != f"bytes {start}-{end}/{total}":
            raise TransferError("RDMVA_DATA_SOURCE_CHANGED_OR_RANGE_IGNORED")
        if response.headers.get("Content-Encoding", "identity").lower() != "identity":
            raise TransferError("RDMVA_DATA_UNEXPECTED_WIRE_COMPRESSION")
        declared = response.headers.get("Content-Length")
        if declared is not None and int(declared) != end - start + 1:
            raise TransferError("RDMVA_DATA_SEGMENT_HTTP_LENGTH_MISMATCH")
        h = hashlib.sha256()
        left = end - start + 1
        with Path(path).open("r+b", buffering=0) as output:
            output.seek(start)
            while left:
                block = response.read(min(MIB, left))
                if not block:
                    raise TransferError("RDMVA_DATA_SEGMENT_TRUNCATED")
                output.write(block)
                h.update(block)
                left -= len(block)
            if response.read(1):
                raise TransferError("RDMVA_DATA_SEGMENT_OVERSIZED")
            os.fsync(output.fileno())
        return h.hexdigest()


def _retry_transient(operation, label):
    for attempt in range(3):
        try:
            return operation()
        except (URLError, OSError, TimeoutError) as exc:
            if isinstance(exc, HTTPError) and exc.code not in RETRYABLE:
                raise
            if attempt == 2:
                raise
            print("RDMVA_DATA_TRANSIENT_RETRY",label,"attempt="+str(attempt+1),flush=True)
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def _single_stream(url, dest, *, expected_bytes, max_bytes, label, expected_sha256=None):
    partial = Path(str(dest) + ".part")
    # Legacy unproven .part content is not mixed with a new remote version.
    # The old incomplete bytes are overwritten only when a new transfer starts.
    received = 0
    with _open(url, {}, timeout=120) as remote, partial.open("wb") as output:
        if remote.status != 200 or remote.headers.get("Content-Encoding", "identity").lower() != "identity":
            raise TransferError("RDMVA_DATA_SINGLE_STREAM_UNEXPECTED_RESPONSE")
        header = remote.headers.get("Content-Length")
        if header and int(header) > max_bytes:
            raise TransferError("RDMVA_DATA_SOURCE_EXCEEDS_ALLOWED_CAP")
        if expected_bytes is not None and header and int(header) != expected_bytes:
            raise TransferError("RDMVA_DATA_OFFICIAL_LENGTH_CHANGED")
        while True:
            block = remote.read(MIB)
            if not block:
                break
            received += len(block)
            if received > max_bytes:
                raise TransferError("RDMVA_DATA_SOURCE_EXCEEDS_ALLOWED_CAP")
            output.write(block)
        output.flush()
        os.fsync(output.fileno())
    if expected_bytes is not None and received != expected_bytes:
        raise TransferError("RDMVA_DATA_SINGLE_STREAM_LENGTH_MISMATCH")
    if expected_sha256 and _digest(partial) != expected_sha256.lower():
        raise TransferError("RDMVA_DATA_OFFICIAL_SHA256_MISMATCH")
    os.replace(partial, dest)
    print("RDMVA_DATA_TRANSFER_COMPLETE", label, "mode=single_https", "bytes="+str(received), flush=True)
    return {"bytes": received, "mode": "single_https"}


def download_public_data(url, destination, *, expected_bytes=None, max_bytes=None,
                         expected_sha256=None, label=None, connections=None,
                         parallel_threshold=DEFAULT_THRESHOLD,
                         segment_size=DEFAULT_SEGMENT, allow_test_http=False):
    """Fetch official public bytes; callers retain their own data-format audit.

    No direct patient data inputs. Does NOT bypass missing upstream SHA pins.
    max_bytes is required for size-changing live scientific sources.
    """
    global _LOOPBACK_TEST_MODE
    _allowed_url(url, allow_test_http=allow_test_http)
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    name = label or target.name
    if expected_bytes is not None and (not isinstance(expected_bytes, int) or expected_bytes <= 0):
        raise ValueError("invalid expected_bytes")
    limit = max_bytes if max_bytes is not None else expected_bytes
    if not isinstance(limit, int) or limit <= 0:
        raise ValueError("RDMVA_DATA_DOWNLOAD_REQUIRES_SOURCE_SIZE_OR_CAP")
    if expected_bytes is not None and expected_bytes > limit:
        raise ValueError("expected_bytes exceeds cap")
    if target.is_file() and target.stat().st_size == expected_bytes and (
        expected_sha256 is None or _digest(target) == expected_sha256.lower()
    ):
        print("RDMVA_DATA_VERIFIED_EXISTING", name, flush=True)
        return {"bytes": target.stat().st_size, "mode": "already_present"}
    count = connections
    if count is None:
        try:
            count = int(os.environ.get("RDMVA_DATA_CONNECTIONS", "4"))
        except ValueError:
            count = 4
    count = max(1, min(int(count), 8))
    if segment_size < 4096:
        raise ValueError("segment size too small")
    if allow_test_http:
        with _lock:
            _LOOPBACK_TEST_MODE = True
    try:
        total, validator, ranged = _retry_transient(
            lambda: _probe(url, expected_bytes=expected_bytes), name)
        if total is not None and total > limit:
            raise TransferError("RDMVA_DATA_SOURCE_EXCEEDS_ALLOWED_CAP")
        # Multi-range only with explicit source identity and verified 206.
        parallel = bool(ranged and validator and total and total >= parallel_threshold and count > 1)
        if not parallel:
            state = Path(str(target) + ".part.rdmva.json")
            if state.exists():
                raise TransferError("RDMVA_DATA_EXISTING_SEGMENT_CHECKPOINT_NEEDS_SAME_RANGE_SOURCE")
            print("RDMVA_DATA_TRANSFER_MODE", name, "single_https", flush=True)
            result = _retry_transient(
                lambda: _single_stream(url,target,expected_bytes=expected_bytes,
                                       max_bytes=limit,label=name,
                                       expected_sha256=expected_sha256), name)
        else:
            part = Path(str(target) + ".part")
            state = Path(str(part) + ".rdmva.json")
            fingerprint = hashlib.sha256(url.encode("utf-8")).hexdigest()
            # A pre-existing legacy partial without provenance cannot be accepted.
            if part.exists() and not state.exists():
                raise TransferError("RDMVA_DATA_LEGACY_PARTIAL_REQUIRES_VERIFIED_SOURCE_RESTART")
            pieces = (total + segment_size - 1) // segment_size
            if state.exists():
                receipt = json.loads(state.read_text(encoding="utf-8"))
                if (receipt.get("schema") != SCHEMA or receipt.get("url_sha256") != fingerprint
                        or receipt.get("bytes") != total or receipt.get("validator") != validator
                        or receipt.get("segment_size") != segment_size
                        or not part.is_file() or part.stat().st_size != total):
                    raise TransferError("RDMVA_DATA_SEGMENT_CHECKPOINT_SOURCE_OR_LENGTH_CHANGED")
            else:
                receipt = {"schema": SCHEMA, "url_sha256": fingerprint, "bytes": total,
                           "validator": validator, "segment_size": segment_size, "done": {}}
                with part.open("wb") as stream:
                    stream.truncate(total)
                _atomic_json(state, receipt)
            for ix, digest in list(receipt["done"].items()):
                start = int(ix) * segment_size
                if (int(ix) >= pieces or start < 0 or not isinstance(digest, str)
                        or _digest(part, start, min(segment_size, total-start)) != digest):
                    receipt["done"].pop(ix, None)
            _atomic_json(state, receipt)
            remaining = [i for i in range(pieces) if str(i) not in receipt["done"]]
            print("RDMVA_DATA_TRANSFER_MODE", name, "segmented_https",
                  "segments="+str(pieces), "connections="+str(count),
                  "resuming="+str(pieces-len(remaining)), flush=True)
            transferred = 0
            last = time.monotonic()
            with ThreadPoolExecutor(max_workers=count) as pool:
                futures = {}
                for ix in remaining:
                    start = ix*segment_size
                    end = min(start + segment_size, total) - 1
                    def attempt(a=start, b=end):
                        for retry in range(3):
                            try:
                                return _one_range(url, part, a, b, total, validator)
                            except (HTTPError, URLError, TimeoutError, ConnectionError, OSError) as exc:
                                if isinstance(exc, HTTPError) and exc.code not in RETRYABLE:
                                    raise
                                if retry == 2:
                                    raise TransferError("RDMVA_DATA_SEGMENT_RETRIES_EXHAUSTED") from exc
                                time.sleep(1 << retry)
                        raise AssertionError("unreachable")
                    futures[pool.submit(attempt)] = (ix, end-start+1)
                for future in as_completed(futures):
                    ix, size = futures[future]
                    digest = future.result()
                    receipt["done"][str(ix)] = digest
                    _atomic_json(state, receipt)
                    transferred += size
                    if time.monotonic() - last >= 5 or len(receipt["done"]) == pieces:
                        print("RDMVA_DATA_DOWNLOAD_PROGRESS", name,
                              f"{100*len(receipt['done'])/pieces:.1f}%",
                              "completed_segments="+str(len(receipt["done"]))+"/"+str(pieces),
                              "bytes_new="+str(transferred), flush=True)
                        last = time.monotonic()
            if len(receipt["done"]) != pieces:
                raise TransferError("RDMVA_DATA_SEGMENTS_INCOMPLETE")
            # Caller verifies full official hashes; pinned archives verify here too.
            if expected_sha256 and _digest(part) != expected_sha256.lower():
                raise TransferError("RDMVA_DATA_OFFICIAL_SHA256_MISMATCH")
            os.replace(part, target)
            state.unlink(missing_ok=True)
            print("RDMVA_DATA_TRANSFER_COMPLETE", name, "mode=segmented_https",
                  "bytes="+str(total), flush=True)
            result = {"bytes": total, "mode": "segmented_https"}
        if expected_sha256 and _digest(target) != expected_sha256.lower():
            raise TransferError("RDMVA_DATA_OFFICIAL_SHA256_MISMATCH")
        return result
    finally:
        if allow_test_http:
            with _lock:
                _LOOPBACK_TEST_MODE = False
