import merlict
import numpy as np
import os
import tempfile
import importlib
from importlib import resources


SCENERY_PATH = os.path.join(
    str(importlib.resources.files("merlict")),
    "tests",
    "resources",
    "segmented_reflector.tar",
)


def test_dump_and_open_dump_again():
    with tempfile.TemporaryDirectory(prefix="merlict_") as tmpdir:
        DUMP_PATH = os.path.join(tmpdir, "dump.bin")

        mli = merlict.open(SCENERY_PATH)
        mli.dump(DUMP_PATH)

        mli2 = merlict.open(DUMP_PATH)
