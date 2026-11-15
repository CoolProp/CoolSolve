#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Write the ZIP fixtures of the [zip] tests (tests/fixtures/zip/).

The GUI opens a model bundle that any zip tool made (docs/gui.md §7.2). Each fixture is the same small
model (model.eescode, model.initials) archived by one tool, the way a user would: Info-ZIP `zip` (Linux,
Deflate with extra fields, streamed with data descriptors, a folder, a password), 7-Zip (Deflate, BZIP2),
bsdtar (libarchive), Python's zipfile (ZIP64, LZMA), the layout of macOS' Compress (__MACOSX/, .DS_Store)
and of Windows' Compressed folder (names in code page 437). tests/test_zip_bundle.cpp reads them.

    python3 -I tools/zip_fixtures.py [--out tests/fixtures/zip]

Needs `zip`, `7z` and `bsdtar` on the PATH; a fixture whose tool is missing is skipped with a message.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))

# Long enough for Deflate to compress it (a tool stores a file that would not shrink)
MODEL = ('"The model of the ZIP fixtures (tests/test_zip_bundle.cpp): T in °C"\n'
         + "".join(f"y_{i} = x^{i} + {i}\n" for i in range(1, 21))
         + "x = 2\n")
INITIALS = "".join(f"y_{i}=1\n" for i in range(1, 21))


def run(cmd, cwd, out_path=None):
    """Run a tool in `cwd`; its standard output goes to `out_path` when given. False when it is missing."""
    if shutil.which(cmd[0]) is None:
        print(f"  skipped ({cmd[0]} not found): {os.path.basename(out_path or cmd[-1])}")
        return False
    if out_path:  # through a pipe, which the tool cannot seek back on
        with open(out_path, "wb") as f:
            f.write(subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, check=True).stdout)
    else:
        subprocess.run(cmd, cwd=cwd, stdout=subprocess.DEVNULL, check=True)
    return True


def write_model(folder, stem="model"):
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, stem + ".eescode"), "w", encoding="utf-8") as f:
        f.write(MODEL)
    with open(os.path.join(folder, stem + ".initials"), "w", encoding="utf-8") as f:
        f.write(INITIALS)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=os.path.join(HERE, "..", "tests", "fixtures", "zip"))
    out = os.path.abspath(parser.parse_args().out)
    os.makedirs(out, exist_ok=True)
    files = ["model.eescode", "model.initials"]

    with tempfile.TemporaryDirectory() as tmp:
        write_model(tmp)
        write_model(os.path.join(tmp, "mymodel"))
        o = lambda name: os.path.join(out, name)
        for name in os.listdir(out):
            if name.endswith(".zip"):
                os.remove(o(name))

        # Info-ZIP (the zip command of Linux): Deflate, extra fields; to a pipe: data descriptors
        run(["zip", "-q", o("infozip_deflate.zip"), *files], tmp)
        run(["zip", "-q", "-", *files], tmp, o("infozip_stream.zip"))
        run(["zip", "-q", "-r", o("infozip_folder.zip"), "mymodel"], tmp)
        run(["zip", "-q", "-P", "secret", o("infozip_encrypted.zip"), *files], tmp)
        # 7-Zip: Deflate by default, BZIP2 on request
        run(["7z", "a", "-tzip", o("7zip_deflate.zip"), *files], tmp)
        run(["7z", "a", "-tzip", "-mm=BZip2", o("7zip_bzip2.zip"), *files], tmp)
        # bsdtar (libarchive, macOS' and Windows' tar): Deflate with data descriptors
        run(["bsdtar", "--format", "zip", "-cf", o("bsdtar.zip"), *files], tmp)

        # Python's zipfile: ZIP64 records, LZMA
        with zipfile.ZipFile(o("python_zip64.zip"), "w", zipfile.ZIP_DEFLATED) as zf:
            for name in files:
                with open(os.path.join(tmp, name), "rb") as src, zf.open(name, "w", force_zip64=True) as dst:
                    dst.write(src.read())
        with zipfile.ZipFile(o("python_lzma.zip"), "w", zipfile.ZIP_LZMA) as zf:
            for name in files:
                zf.write(os.path.join(tmp, name), name)

        # macOS' Compress of a folder: the folder, __MACOSX/ with the AppleDouble files, .DS_Store
        with zipfile.ZipFile(o("macos_folder.zip"), "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("mymodel/", b"")
            for name in files:
                zf.write(os.path.join(tmp, name), "mymodel/" + name)
            zf.writestr("mymodel/.DS_Store", b"\x00\x00\x00\x01Bud1" + bytes(64))
            zf.writestr("__MACOSX/", b"")
            zf.writestr("__MACOSX/mymodel/", b"")
            zf.writestr("__MACOSX/mymodel/._model.eescode", b"\x00\x05\x16\x07\x00\x02\x00\x00" + bytes(74))

        # Windows' Compressed folder: an accented name in code page 437, without the UTF-8 flag. zipfile
        # writes ASCII names as they are, so the names are written with a placeholder and patched.
        placeholder, stem = b"#changeur", "Échangeur"
        with zipfile.ZipFile(o("windows_cp437.zip"), "w", zipfile.ZIP_DEFLATED) as zf:
            for name in files:
                zf.write(os.path.join(tmp, name), name.replace("model", placeholder.decode()))
        with open(o("windows_cp437.zip"), "rb") as f:
            data = f.read()
        with open(o("windows_cp437.zip"), "wb") as f:
            f.write(data.replace(placeholder, stem.encode("cp437")))

    for name in sorted(os.listdir(out)):
        if name.endswith(".zip"):
            print(f"  {name}: {os.path.getsize(os.path.join(out, name))} bytes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
