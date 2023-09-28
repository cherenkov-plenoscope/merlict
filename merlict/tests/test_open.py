import merlict
import numpy as np
import pkg_resources
import os
import tempfile

SCENERY_PATH = pkg_resources.resource_filename(
    package_or_requirement="merlict",
    resource_name=os.path.join(
        "tests", "resources", "segmented_reflector.tar"
    ),
)


def test_dump_and_open_dump_again():
    with tempfile.TemporaryDirectory(prefix="merlict_") as tmpdir:
        DUMP_PATH = os.path.join(tmpdir, "dump.bin")

        scn = merlict.open(SCENERY_PATH)
        scn.dump(DUMP_PATH)

        scn2 = merlict.open(DUMP_PATH)
