from .version import __version__
from . import c89
from . import materials
from . import scenery
from . import ray
from . import photon
from . import intersection
from . import intersectionSurfaceNormal


def open(path):
    """
    Parameters
    ----------
    path : str
        Path to a scenery. This can be either the path to a tape-archive
        `.tar` containing a scenery, or it can be the path to merlict's
        binary memory-dump of a scenery.
    """
    return c89.wrapper.Merlict(path=path)


def compile(sceneryPy):
    """
    Parameters
    ----------
    sceneryPy : dict
        The scenery represented only with python builtin dicts, lists and
        optional numpy.arrays.
    """
    sceneryStr = scenery.string_format.convert.sceneryPy_to_sceneryStr(
        sceneryPy=sceneryPy, indent=4, relations_indent=0
    )

    return c89.wrapper.Merlict(sceneryStr=sceneryStr)
