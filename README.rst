|ImgMerlictPythonLogo|

|TestStatus| |PyPiStatus| |BlackStyle| |PackStyleBlack| |GPLv3Logo|


More light than you can handle! Also: This is in beta state. Don't judge me!


*******
Install
*******

.. code-block:: bash

    pip install merlict


***************
Minimal example
***************

.. code-block:: python

    import merlict

    scenery = merlict.c89.wrapper.Server(path="merlict/tests/resources/segmented_reflector.tar")
    scenery.view()


.. |BlackStyle| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

.. |TestStatus| image:: https://github.com/cherenkov-plenoscope/merlict/actions/workflows/test.yml/badge.svg?branch=main
    :target: https://github.com/cherenkov-plenoscope/merlict/actions/workflows/test.yml

.. |PyPiStatus| image:: https://img.shields.io/pypi/v/merlict
    :target: https://pypi.org/project/merlict

.. |PackStyleBlack| image:: https://img.shields.io/badge/pack%20style-black-000000.svg
    :target: https://github.com/cherenkov-plenoscope/black_pack

.. |GPLv3Logo| image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg
    :target: https://www.gnu.org/licenses/gpl-3.0

.. |ImgMerlictPythonLogo| image:: https://github.com/cherenkov-plenoscope/merlict/blob/main/readme/merlict-python-logo-inkscape.png?raw=True

