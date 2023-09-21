|ImgMerlictPythonLogo|

|TestStatus| |PyPiStatus| |BlackStyle| |PackStyleBlack| |GPLv3Logo|


More light than you can handle! Also: This is in beta state. Don't judge me!


*******
Install
*******

.. code-block::

    pip install merlict


****************************************
Highly techincal talking about sceneries
****************************************

``sceneryPy``
---------------
A ``python`` ``dict``. Objects are ``dict``s with their vertices and faces in ``numpy.arrays`` or ``lists``. Geometric relations between objects and materials are also ``dicts``.

``sceneryDs``
---------------
A ``py``-dictionary where both the keys and the payload are ``str``-strings.
The objects are represented in wavefront-files (``.obj``). Materials are represented in ``JSON``-strings.
Lossless: [``sceneryTr``, ``sceneryUs``, ``sceneryDr``]

``sceneryTr``
---------------
A tape-archive (``.tar``). The filenames inside are the keys, and the files payloads are the values used in ``sceneryDs``.
Lossless: [``sceneryDs``, ``sceneryUs``, ``sceneryDr``]

``sceneryDr``
---------------
The same as ``sceneryTr`` but dumped into the filesystem, i.e. into directories and files according to the tape-archive.
Lossless: [``sceneryDs``, ``sceneryUs``, ``sceneryTr``]

``sceneryUs``
---------------
Merlict's in-memory-represenation used to build the oc-trees for acceleration.
It still has the hirachy defined in ``geometry/relations``.
Lossless: [``sceneryDs``, ``sceneryTr``]
Lossy: [``sceneryCa``]

``sceneryCa``
---------------
Merlict's in-memory-represenation used for the ray tracing. It is based on oc-trees. The object's faces are `flattened` into the scenery and then bundled in global oc-trees. This is good for cache-aware access. A seperate look-up-table is created to resolve the user's ``geometry/relations`` but is not needed when ray-face-intersections are queried.
Lossless: [``sceneryCs``]
Lossy: []

``sceneryCs``
---------------
Merlict can serialize and dump its ``sceneryCa`` into a binary blob.
Since the population of oc-trees in the conversion from ``sceneryUs`` to ``sceneryCa`` might take a significant time, ``sceneryCs`` is used to dump and load the scenery into caches.
This is only expected to work on a single machine. This is not meant to share sceneries and is expected to fail when read by machines with different architectures.
Lossless: [``sceneryCa``]


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

