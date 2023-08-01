import numpy as np


def isdtype(r, dtype):
    """
    Returns True if recarray r has the expected dtype.

    Parameters
    ----------
    r : np.core.records.recarray
        Recarray to check
    dtype : list of tuples(str, dtype)
        Expected dtype
    """
    if not isinstance(r, np.core.records.recarray):
        return False

    if len(dtype) != len(r.dtype.names):
        return False

    for i in range(len(dtype)):
        expected_dtype_name_i = dtype[i][0]
        expected_dtype_i = dtype[i][1]
        if expected_dtype_name_i != r.dtype.names[i]:
            return False
        if expected_dtype_i != r.dtype[i]:
            return False

    return True


def zeros(size, dtype):
    out = np.core.records.recarray(shape=size, dtype=dtype)
    for key in out.dtype.names:
        out[key] = np.zeros(size, dtype=out.dtype[key])
    return out
