|ImgMerlictPythonLogo|

|TestStatus| |PyPiStatus| |BlackStyle| |BlackPackStyle| |GPLv3LicenseBadge|


More light than you can handle! Also: This is in beta state. Don't judge me!
Merlict would not exist without the author's past and present affiliations:

- Max-Planck-Institute for Nuclear Physics,
  Saupfercheckweg 1, 69117 Heidelberg, Germany

- Institute for Particle Physics and Astrophysics,
  ETH-Zurich, Otto-Stern-Weg 5, 8093 Zurich, Switzerland

- Experimental Physics Vb, Astroparticle Physics,
  TU-Dortmund, Otto-Hahn-Str. 4a, 44227 Dortmund, Germany


*******
Install
*******

.. code-block:: bash

    pip install merlict


***************
Minimal example
***************


Load an existing scenery
========================

.. code-block:: python

    import merlict
    import pkg_resources
    import os

    sceneryPy = merlict.scenery.read_tar(
        path=pkg_resources.resource_filename(
            "merlict",
            os.path.join("tests", "resources", "segmented_reflector.tar"),
        ),
    )

    mli = merlict.compile(sceneryPy)
    mli.view()


Query the intersection of rays with the scenery
===============================================

.. code-block:: python

    rays = merlict.ray.init(size=1)  # only one ray for demonstration

    rays["support.x"] = 0.3
    rays["support.y"] = 0.1
    rays["support.z"] = 2.3

    rays["direction.x"] = 0
    rays["direction.y"] = 0
    rays["direction.z"] = -1

    hits, intersections = mli.query_intersectionSurfaceNormal(rays)

    if hits[0]:
        print("The ray intersects with the surface")
        print("of (object, face): ({:d}, {:d}),".format(
                intersections[0]["geometry_id.robj"],
                intersections[0]["geometry_id.face"],
            )
        )
        print("in distance: {:f},".format(intersections[0]["distance_of_ray"]))
        print("at position: ({:f}, {:f}, {:f}),".format(
                intersections[0]["position.x"],
                intersections[0]["position.y"],
                intersections[0]["position.z"],
            )
        )
        print("and with surface-normal: ({:f}, {:f}, {:f}).".format(
                intersections[0]["surface_normal.x"],
                intersections[0]["surface_normal.y"],
                intersections[0]["surface_normal.z"],
            )
        )
        print(
            "With respect to the frame "
            "of the intersected object the ray intersected"
        )
        print("at position: ({:f}, {:f}, {:f}),".format(
                intersections[0]["position_local.x"],
                intersections[0]["position_local.y"],
                intersections[0]["position_local.z"],
            )
        )
        print("and with surface-normal: ({:f}, {:f}, {:f}).".format(
                intersections[0]["surface_normal_local.x"],
                intersections[0]["surface_normal_local.y"],
                intersections[0]["surface_normal_local.z"],
            )
        )
    else:
        print("The ray does not intersect with any surface.")


will yield


.. code-block::

    The ray intersects with the surface
    of (object, face): (28, 35),
    in distance: 2.287463,
    at position: (0.300000, 0.100000, 0.012537),
    and with surface-normal: (-0.075000, -0.025000, 0.996870).
    With respect to the frame of the intersected object the ray intersected
    at position: (0.009082, 0.044013, 0.000270),
    and with surface-normal: (-0.002270, -0.011003, 0.999937).



.. |BlackStyle| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

.. |TestStatus| image:: https://github.com/cherenkov-plenoscope/merlict/actions/workflows/test.yml/badge.svg?branch=main
    :target: https://github.com/cherenkov-plenoscope/merlict/actions/workflows/test.yml

.. |PyPiStatus| image:: https://img.shields.io/pypi/v/merlict
    :target: https://pypi.org/project/merlict

.. |BlackPackStyle| image:: https://img.shields.io/badge/pack%20style-black-000000.svg
    :target: https://github.com/cherenkov-plenoscope/black_pack

.. |GPLv3LicenseBadge| image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg
    :target: https://www.gnu.org/licenses/gpl-3.0

.. |ImgMerlictPythonLogo| image:: https://github.com/cherenkov-plenoscope/merlict/blob/main/readme/merlict-python-logo-inkscape.png?raw=True

