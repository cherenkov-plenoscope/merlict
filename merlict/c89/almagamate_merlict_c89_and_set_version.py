"""
Almagamate the merlict_c89 sources and set the version
This is only ever executed once by the developers.
"""
import os
import subprocess

# merlict python package's version
# --------------------------------
MERLICT_PYTHON_VERSION_STR = "0.0.8"

merlict_c89_libs = ["chk_debug", "mli_core", "mli_viewer"]

subprocess.call(
    [
        "python",
        os.path.join(".", "merlict_c89", "tools", "mli_almagamate.py"),
        ".",
    ] + [os.path.join(".", "merlict_c89", l) for l in merlict_c89_libs]
)

src_name = str.join("_", merlict_c89_libs)
src_h = src_name + ".h"
src_c = src_name + ".c"

# gather c89 merlict version
# --------------------------
MERLICT_C89_VERSION = {
    "MLI_VERSION_MAYOR": -1,
    "MLI_VERSION_MINOR": -1,
    "MLI_VERSION_PATCH": -1,
}
MERLICT_C89_VERSION_DIGIT_POS = len("#define MLI_VERSION_MAYOR ")

with open(src_h, "rt") as f:
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

# combine versions
# ----------------
VERSION_STR = "{:s}.{:s}".format(
    MERLICT_PYTHON_VERSION_STR,
    MERLICT_C89_VERSION_STR,
)

# export version
# --------------
with open(os.path.join(os.path.join("..", "version.py")), "wt") as f:
    f.write("# I was written by: ")
    f.write("merlict/version/set_and_make_version_script.py\n")
    f.write('__version__ = "' + VERSION_STR + '"')
    f.write("\n")
