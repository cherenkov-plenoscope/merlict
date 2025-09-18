"""
Almagamate the merlict_c89 sources and set the version
This is only ever executed when the merlict_c89 sources change.
This is not expected to be executed with every pip install.

```bash
you@com: merlict/merlict/c89$ python almagamate_merlict_c89_and_set_version.py
```
"""

import os
import glob
import shutil
import subprocess


def rmtree(path):
    try:
        shutil.rmtree(path)
        print(f"rm -r {path:s}")
    except FileNotFoundError as e:
        pass


def rm(path):
    try:
        os.remove(path)
        print(f"rm {path:s}")
    except FileNotFoundError as e:
        pass


print(
    "Almagamate the merlict_c89 sources\n"
    "==================================\n"
)

# list the libs from within merlict_c89 which will be almagamated.
merlict_c89_module_paths = glob.glob(
    os.path.join(".", "merlict_c89", "src", "*")
)
merlict_c89_header_path = "mli.h"
merlict_c89_source_path = "mli.c"

# Remove all old installations and caches
# ---------------------------------------

print("1. ) uninstalling merlict python package ...")
subprocess.call(["pip", "uninstall", "--yes", "merlict"])

print(
    "2.1 ) removing old builds, dists, test caches, and cython artifacts ..."
)
merlict_dir = os.path.join("..", "..")
rmtree(os.path.join(merlict_dir, "build"))
rmtree(os.path.join(merlict_dir, "dist"))
rmtree(os.path.join(merlict_dir, "merlict.egg-info"))
rmtree(os.path.join(merlict_dir, ".pytest_cache"))
rmtree(os.path.join(merlict_dir, "merlict", "__pycache__"))
rmtree(os.path.join(merlict_dir, "merlict", "c89", "__pycache__"))

print("2.2 ) remove old almagamated sources ...")
rm(os.path.join(merlict_dir, "merlict", "c89", merlict_c89_header_path))
rm(os.path.join(merlict_dir, "merlict", "c89", merlict_c89_source_path))

print("2.3 ) remove old cython code ...")
rm(os.path.join(merlict_dir, "merlict", "c89", "wrapper.c"))


print("3. ) almagamate the sources from merlict_c89 ...")
_outdir = "."
subprocess.call(
    [
        "python",
        os.path.join(".", "merlict_c89", "tools", "almagamate.py"),
        _outdir,
    ]
    + merlict_c89_module_paths
)

print("4. ) gathering merlict_c89 version from almagamated sources ...")
MERLICT_C89_VERSION = {
    "MLI_VERSION_MAYOR": -1,
    "MLI_VERSION_MINOR": -1,
    "MLI_VERSION_PATCH": -1,
}
MERLICT_C89_VERSION_DIGIT_POS = len("#define MLI_VERSION_MAYOR ")

with open(merlict_c89_header_path, "rt") as f:
    txt = f.read()
    keys = list(MERLICT_C89_VERSION.keys())
    for line in str.splitlines(txt):
        for key in keys:
            if key in line:
                MERLICT_C89_VERSION[key] = int(
                    line[MERLICT_C89_VERSION_DIGIT_POS:]
                )

MERLICT_C89_VERSION_STR = "{:d}.{:d}.{:d}".format(
    MERLICT_C89_VERSION["MLI_VERSION_MAYOR"],
    MERLICT_C89_VERSION["MLI_VERSION_MINOR"],
    MERLICT_C89_VERSION["MLI_VERSION_PATCH"],
)

print(
    f"5. ) writing merlict_c89 version '{MERLICT_C89_VERSION_STR:s}' "
    "to be parsed by python package ..."
)
with open(os.path.join("..", "version_merlict_c89.py"), "wt") as f:
    f.write("# I was written by: ")
    f.write("merlict/c89/almagamate_merlict_c89_and_set_version.py. ")
    f.write("Do not modify me manually.\n")
    f.write('__version__ = "' + MERLICT_C89_VERSION_STR + '"')
    f.write("\n")

print("SUCCESS")
print("Now you can 'pip install ./merlict'.")
