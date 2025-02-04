from . import recarray
from . import resources
import numpy as np


def is_all_greater_zero(x):
    return bool(np.all(x > 0))


def is_monotonically_increasing(x):
    return is_all_greater_zero(np.gradient(x))


def is_free_of_nan(x):
    return not bool(np.any(np.isnan(x)))
