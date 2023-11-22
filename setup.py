import setuptools
import os
import Cython
from Cython import Build as _
import numpy


with open("README.rst", "r", encoding="utf-8") as f:
    long_description = f.read()


with open(os.path.join("merlict", "version.py")) as f:
    txt = f.read()
    last_line = txt.splitlines()[-1]
    version_string = last_line.split()[-1]
    version = version_string.strip("\"'")


extensions = [
    setuptools.Extension(
        name="merlict.c89.wrapper",
        sources=[
            os.path.join("merlict", "c89", "wrapper.pyx"),
            os.path.join("merlict", "c89", "bridge.c"),
            os.path.join("merlict", "c89", "merlict_c89.c"),
        ],
        language="c",
        include_dirs=[numpy.get_include()],
    )
]

setuptools.setup(
    name="merlict",
    version=version,
    description="Ray tracing in python",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/cherenkov-plenoscope/merlict",
    author="Sebastian Achim Mueller",
    author_email="sebastian-achim.mueller@mpi-hd.mpg.de",
    packages=[
        "merlict",
        "merlict.materials",
        "merlict.materials.surfaces",
        "merlict.materials.media",
        "merlict.photon",
        "merlict.c89",
        "merlict.utils",
        "merlict.intersectionSurfaceNormal",
        "merlict.intersection",
        "merlict.ray",
        "merlict.scenery",
        "merlict.scenery.string_format",
    ],
    package_data={
        "merlict": [
            os.path.join("tests", "resources", "*"),
        ],
        "merlict.materials.media": [
            os.path.join("resources", "*"),
        ],
        "merlict.materials.surfaces": [
            os.path.join("resources", "*"),
        ],
        "merlict.c89": [
            os.path.join("*.pyx"),
            os.path.join("*.pxd"),
            "bridge.h",
            "bridge.c",
            "merlict_c89.h",
            "merlict_c89.c",
        ],
    },
    install_requires=[
        "json_numpy_sebastian-achim-mueller",
        "triangle_mesh_io>=0.0.2",
    ],
    ext_modules=Cython.Build.cythonize(extensions, language_level=3),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
