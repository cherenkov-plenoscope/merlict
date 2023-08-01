import setuptools
import numpy
import os


with open("README.md", "r") as f:
    long_description = f.read()


_version_path = os.path.join("merlict", "version", "__init__.py")
with open(_version_path, "r") as f:
    _txt = f.read()
    _version_line = _txt.splitlines()[-1]
    _version_expression = _version_line.split()[-1]
    version = _version_expression.strip("\"'")


setuptools.setup(
    name="merlict_sebastian-achim-mueller",
    version=version,
    description="Ray tracing in python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cherenkov-plenoscope/merlict.git",
    author="Sebastian Achim Mueller",
    author_email="sebastian-achim.mueller@mpi-hd.mpg.de",
    packages=["merlict"],
    package_data={
        "merlict": [
            os.path.join("test", "resources", "*"),
            os.path.join("materials", "media", "*"),
            os.path.join("materials", "surfaces", "*"),
        ]
    },
    install_requires=[
        "setuptools>=18.0",
        "cython",
        "json_numpy_sebastian-achim-mueller",
    ],
    zip_safe=False,
    ext_modules=[
        setuptools.Extension(
            "merlict.c89.wrapper",
            sources=[
                os.path.join("merlict", "c89", "wrapper.pyx"),
                os.path.join(
                    "merlict", "c89", "merlict_c89", "merlict", "mli.c"
                ),
                os.path.join(
                    "merlict", "c89", "merlict_c89", "merlict", "mli_viewer.c"
                ),
                os.path.join("merlict", "c89", "bridge.c"),
            ],
            include_dirs=[
                numpy.get_include(),
                os.path.join("merlict", "c89"),
            ],
            language="c",
        ),
    ],
)
