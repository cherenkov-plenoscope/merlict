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
        name="merlict.c89",
        sources=[
            os.path.join("merlict", "c89", "__init__.pyx"),
            os.path.join("merlict", "c89", "bridge.c"),
        ],
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
        "merlict.photon",
        "merlict.c89",
        "merlict.viewer",
        "merlict.utils",
        "merlict.intersectionSurfaceNormal",
        "merlict.intersection",
        "merlict.ray",
        "merlict.scenery",
    ],
    package_data={
        "merlict": [
            os.path.join("tests", "resources", "*"),
        ],
        "merlict.materials": [
            os.path.join("media", "*"),
            os.path.join("surfaces", "*"),
        ],
        "merlict.c89": [
            os.path.join("merlict_c89", "merlict", "*"),
        ]
    },
    install_requires=[
        "json_numpy_sebastian-achim-mueller",
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


"""
"merlict.c89": [
    os.path.join("*"),
    os.path.join("merlict_c89", "*"),
    os.path.join("merlict_c89", "merlict", "*"),
],
"""