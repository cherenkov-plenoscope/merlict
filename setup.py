import setuptools
import numpy
import os


with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="merlict",
    version="0.0.0",
    description="Ray tracing in python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cherenkov-plenoscope/merlict.git",
    author="Sebastian Achim Mueller",
    author_email="sebastian-achim.mueller@mpi-hd.mpg.de",
    packages=["merlict"],
    package_data={"merlict": []},
    install_requires=["setuptools>=18.0", "cython",],
    zip_safe=False,
    ext_modules=[
        setuptools.Extension(
            "merlict.c89.wrapper",
            sources=[
                os.path.join(
                    "merlict", "c89", "wrapper.pyx"
                ),
                os.path.join(
                    "merlict", "c89", "merlict_c89", "merlict", "mli.c"
                ),
                os.path.join(
                    "merlict", "c89", "merlict_c89", "merlict", "mli_viewer.c"
                ),
                os.path.join(
                    "merlict", "c89", "bridge.c"
                ),
            ],
            include_dirs=[
                numpy.get_include(),
                os.path.join("merlict", "c89"),
            ],
            language="c",
        ),
    ],
)
