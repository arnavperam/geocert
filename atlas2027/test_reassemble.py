"""Archive integrity and safe extraction checks for the delivery helper."""
import hashlib
import importlib.util
import json
from pathlib import Path
import stat
import tempfile
import unittest
import zipfile

spec = importlib.util.spec_from_file_location("reassemble", Path(__file__).with_name("reassemble.py"))
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)


class ArchiveTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def make(self, member="package/file.txt", symlink=False):
        archive = self.root / "original.zip"
        with zipfile.ZipFile(archive, "w") as zipped:
            info = zipfile.ZipInfo(member)
            if symlink:
                info.create_system = 3
                info.external_attr = (stat.S_IFLNK | 0o777) << 16
            zipped.writestr(info, b"exact scientific result\n")
        data = archive.read_bytes()
        parts = self.root / "artifact_parts"
        parts.mkdir()
        manifest = {"archive": "merged.zip", "root": "package", "sha256": hashlib.sha256(data).hexdigest(), "size": len(data), "parts": []}
        for i, chunk in enumerate([data[:40], data[40:]]):
            name = f"part{i:03d}"
            (parts / name).write_bytes(chunk)
            manifest["parts"].append({"name": name, "sha256": hashlib.sha256(chunk).hexdigest(), "size": len(chunk)})
        path = parts / "manifest.json"
        path.write_text(json.dumps(manifest))
        return path

    def test_valid_parts_and_extract(self):
        archive, root = helper.reassemble(self.make())
        self.assertEqual(archive.read_bytes(), (self.root / "original.zip").read_bytes())
        dest = self.root / "extracted"
        helper.extract(archive, dest, root)
        self.assertEqual((dest / "package/file.txt").read_bytes(), b"exact scientific result\n")

    def test_modified_part_rejected_without_output(self):
        manifest = self.make()
        (manifest.parent / "part000").write_bytes(b"malicious")
        with self.assertRaisesRegex(ValueError, "part checksum"):
            helper.reassemble(manifest)
        self.assertFalse((self.root / "merged.zip").exists())
        self.assertEqual(list(self.root.glob("*.partial")), [])

    def test_reordered_parts_rejected(self):
        manifest = self.make()
        value = json.loads(manifest.read_text())
        value["parts"].reverse()
        manifest.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, "complete archive"):
            helper.reassemble(manifest)

    def test_traversal_part_rejected(self):
        manifest = self.make()
        value = json.loads(manifest.read_text())
        value["parts"][0]["name"] = "../original.zip"
        manifest.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, "invalid part"):
            helper.reassemble(manifest)

    def test_traversal_archive_rejected_before_writing(self):
        archive, root = helper.reassemble(self.make("package/../../outside.txt"))
        dest = self.root / "extracted"
        with self.assertRaisesRegex(ValueError, "unsafe archive"):
            helper.extract(archive, dest, root)
        self.assertFalse(dest.exists())

    def test_absolute_archive_rejected(self):
        archive, root = helper.reassemble(self.make("/package/file.txt"))
        with self.assertRaisesRegex(ValueError, "unsafe archive"):
            helper.extract(archive, self.root / "extracted", root)

    def test_wrong_root_rejected(self):
        archive, root = helper.reassemble(self.make("other/file.txt"))
        with self.assertRaisesRegex(ValueError, "unsafe archive"):
            helper.extract(archive, self.root / "extracted", root)

    def test_symlink_rejected(self):
        archive, root = helper.reassemble(self.make(symlink=True))
        with self.assertRaisesRegex(ValueError, "unsafe archive"):
            helper.extract(archive, self.root / "extracted", root)

    def test_existing_content_preserved(self):
        archive, root = helper.reassemble(self.make())
        dest = self.root / "extracted"
        dest.mkdir()
        (dest / "existing").write_text("keep")
        with self.assertRaisesRegex(ValueError, "empty directory"):
            helper.extract(archive, dest, root)
        self.assertEqual((dest / "existing").read_text(), "keep")

    def test_duplicate_member_rejected(self):
        path = self.make()
        archive, root = helper.reassemble(path)
        with zipfile.ZipFile(archive, "a") as zipped:
            zipped.writestr("package/file.txt", "duplicate")
        with self.assertRaisesRegex(ValueError, "duplicate archive"):
            helper.extract(archive, self.root / "extracted", root)


if __name__ == "__main__":
    unittest.main()
