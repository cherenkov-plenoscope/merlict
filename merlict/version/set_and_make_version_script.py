"""
Set and create the version
"""
import os

# merlict python package's version
# --------------------------------
MERLICT_PYTHON_VERSION_STR = "0.0.2"

c89_src_dir = os.path.join("..", "c89", "merlict_c89", "merlict")

# gather c89 merlict version
# --------------------------
MERLICT_C89_VERSION = {
    "MAYOR": -1,
    "MINOR": -1,
    "PATCH": -1,
}
MERLICT_C89_VERSION_DIGIT_POS = len("#define MLI_VERSION_MAYOR ")

with open(os.path.join(c89_src_dir, "mli_version.h"), "rt") as f:
    txt = f.read()
    keys = list(MERLICT_C89_VERSION.keys())
    for line in str.splitlines(txt):
        for key in keys:
            if key in line:
                MERLICT_C89_VERSION[key] = int(
                    line[MERLICT_C89_VERSION_DIGIT_POS:]
                )

MERLICT_C89_VERSION_STR = "{:d}.{:d}.{:d}".format(
    MERLICT_C89_VERSION["MAYOR"],
    MERLICT_C89_VERSION["MINOR"],
    MERLICT_C89_VERSION["PATCH"],
)

# combine versions
# ----------------
VERSION_STR = "{:s}.{:s}".format(
    MERLICT_PYTHON_VERSION_STR, MERLICT_C89_VERSION_STR,
)

# export version
# --------------
version_path = os.path.join("__init__.py")

with open(os.path.join(version_path), "wt") as f:
    f.write("# I was written by: ")
    f.write("merlict/version/set_and_make_version_script.py\n")
    f.write('__version__ = "' + VERSION_STR + '"')
    f.write("\n")