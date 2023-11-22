"""
Almagamate the merlict_c89 sources and set the version
This is only ever executed once by the developers.

```bash
you@com: merlict/merlict/c89$ python almagamate_merlict_c89_and_set_version.py
```
"""
import os
import shutil
import subprocess


def rmtree(path):
    try:
        shutil.rmtree(path)
    except FileNotFoundError as e:
        print(e)


def rm(path):
    try:
        os.remove(path)
    except FileNotFoundError as e:
        print(e)


# list the libs from within merlict_c89 which will be almagamated.
merlict_c89_libs = ["chk_debug", "mli_core", "mli_viewer"]
merlict_c89_header_path = "merlict_c89.h"
merlict_c89_source_path = "merlict_c89.c"

# Remove all old installations and caches
# ---------------------------------------

# uninstall python-package
subprocess.call(["pip", "uninstall", "merlict"])

# remove old builds, dists, test-caches, and cython-artifacts
merlict_dir = os.path.join("..", "..")
rmtree(os.path.join(merlict_dir, "build"))
rmtree(os.path.join(merlict_dir, "dist"))
rmtree(os.path.join(merlict_dir, "merlict.egg-info"))
rmtree(os.path.join(merlict_dir, ".pytest_cache"))
rmtree(os.path.join(merlict_dir, "merlict", "__pycache__"))
rmtree(os.path.join(merlict_dir, "merlict", "c89", "__pycache__"))

# remove the old almagamated sources
rm(os.path.join(merlict_dir, "merlict", "c89", merlict_c89_header_path))
rm(os.path.join(merlict_dir, "merlict", "c89", merlict_c89_source_path))

# remove the cython-code
rm(os.path.join(merlict_dir, "merlict", "c89", "wrapper.c"))

# Set merlict python package's version
# ------------------------------------
MERLICT_PYTHON_VERSION_STR = "0.1.5"  # <- set version here.

# Almagamate the sources from merlict
# -----------------------------------
subprocess.call(
    [
        "python",
        os.path.join(".", "merlict_c89", "tools", "almagamate.py"),
        "--header",
        "./merlict_c89.h",
        "--source",
        "./merlict_c89.c",
    ]
    + [os.path.join(".", "merlict_c89", l) for l in merlict_c89_libs]
)

# automatically gather merlict_c89 version
# ----------------------------------------
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

# combine python-version with c89-version
# ---------------------------------------
VERSION_STR = "{:s}.{:s}".format(
    MERLICT_PYTHON_VERSION_STR,
    MERLICT_C89_VERSION_STR,
)

# write version
# -------------
with open(os.path.join("..", "version.py"), "wt") as f:
    f.write("# I was written by: ")
    f.write("merlict/c89/almagamate_merlict_c89_and_set_version.py\n")
    f.write('__version__ = "' + VERSION_STR + '"')
    f.write("\n")
