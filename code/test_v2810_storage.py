"""Data-storage safety contracts for auto-downloaded free public resources.

This file is synthetic and contains no actual patient's genetic or medical data.
"""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase,main
from unittest.mock import patch
from subprocess import CompletedProcess
from v2810 import resource_bootstrap as bootstrap

GiB=1024**3

class ResourceStoragePolicyTests(TestCase):
 def test_sparse_wsl_vhd_uses_host_capacity_not_virtual(self):
    context={"home":"/home/test","distribution":"Ubuntu"}
    with TemporaryDirectory() as td,\
         patch.object(bootstrap,"_run_wsl",return_value=""),\
         patch.object(bootstrap,"_free_bytes",return_value=500*GiB),\
         patch.object(bootstrap,"_physical_wsl_host_free_bytes",return_value=18*GiB):
        selected=bootstrap._public_storage_location(context,td)
        self.assertEqual(selected["free_bytes"],18*GiB)
        self.assertEqual(selected["reserve_bytes"],12*GiB)
        self.assertEqual(selected["location"],"WSL_LINUX_VHDX_WITH_PHYSICAL_HOST_CAPACITY_GATE")
 def test_external_drive_is_verified_before_download(self):
    context={"home":"/home/test","distribution":"Ubuntu"}
    usb={"Name":"TEST USB DEVICE","Bus":"USB","FileSystem":"NTFS",
         "DiskBytes":2000398934016,"VolumeBytes":2000396832768,"FreeBytes":1700*GiB}
    with TemporaryDirectory() as td:
        path=Path(td)/"config";path.mkdir()
        (path/"local_public_storage_v2810.json").write_text(json.dumps({
           "schema":"V2810_LOCAL_PUBLIC_STORAGE_1",
           "windows_root":"Z:\\Test_Research_References",
           "expected_disk_name":"TEST USB DEVICE",
           "reserve_free_gib":1000,
           "require_external_usb":True
        }),encoding="utf-8")
        process=CompletedProcess([],0,json.dumps(usb),"")
        with patch.object(bootstrap.subprocess,"run",return_value=process),\
             patch.object(bootstrap,"_run_wsl",return_value="Z:"),\
             patch.object(bootstrap,"_free_bytes",return_value=1700*GiB):
            selected=bootstrap._public_storage_location(context,td)
        self.assertEqual(selected["location"],"EXPLICIT_VERIFIED_WINDOWS_REMOVABLE_STORAGE")
        self.assertEqual(selected["reserve_bytes"],1000*GiB)
        self.assertTrue(selected["windows_root"].startswith("Z:\\"))
        self.assertTrue(selected["linux_root"].startswith("/mnt/z/"))
 def test_external_disk_change_refuses_download(self):
    context={"home":"/home/test","distribution":"Ubuntu"}
    changed={"Name":"ANOTHER USB DRIVE","Bus":"USB","FileSystem":"NTFS",
             "VolumeBytes":2000*GiB,"FreeBytes":1700*GiB}
    with TemporaryDirectory() as td:
        path=Path(td)/"config";path.mkdir()
        (path/"local_public_storage_v2810.json").write_text(json.dumps({
           "schema":"V2810_LOCAL_PUBLIC_STORAGE_1",
           "windows_root":"Z:\\Test_Research_References",
           "expected_disk_name":"TEST USB DEVICE",
           "reserve_free_gib":1000,"require_external_usb":True
        }),encoding="utf-8")
        with patch.object(bootstrap.subprocess,"run",return_value=CompletedProcess([],0,json.dumps(changed),"")):
            with self.assertRaisesRegex(RuntimeError,"DEVICE_IDENTITY_CHANGED"):
                bootstrap._public_storage_location(context,td)
 def test_physical_low_capacity_cannot_download_public_chromosome(self):
    case={"genes":{"GENE_A":[{"chrom":"chr21"}]}}
    free=18*GiB;reserve=12*GiB
    fake_dest={"linux_root":"/home/test/resources/gnomad/v4.1",
       "windows_root":None,"free_bytes":free,
       "reserve_bytes":reserve,
       "location":"WSL_LINUX_VHDX_WITH_PHYSICAL_HOST_CAPACITY_GATE"}
    with TemporaryDirectory() as td,\
         patch.object(bootstrap,"offline_external_bin_wsl",return_value=True),\
         patch.object(bootstrap,"wsl_context",return_value={"home":"/home/test","distribution":"Ubuntu"}),\
         patch.object(bootstrap,"_public_storage_location",return_value=fake_dest),\
         patch.object(bootstrap,"_run_wsl",return_value=""),\
         patch.object(bootstrap,"_free_bytes",return_value=500*GiB),\
         patch.object(bootstrap,"_exists",return_value=False),\
         patch.object(bootstrap,"_download_meta",return_value=30*GiB),\
         patch.object(bootstrap,"_curl_download") as forbidden:
         manifest,report=bootstrap.prepare(case,td,download=True,max_gib=100)
         forbidden.assert_not_called()
         self.assertIsNone(manifest)
         self.assertEqual(report["completion"],"PARTIAL_STORAGE_CONSTRAINT")
         self.assertTrue(all(job["status"]=="BLOCKED_STORAGE_BUDGET" for job in report["records"]))

if __name__=="__main__":main()
