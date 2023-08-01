import numpy as np


def isdtype(recarray, dtype):
    if not isinstance(recarray, np.core.records.recarray):
        return False
    for expected_name, expected_dtype in dtype:
        if expected_name not in recarray.dtype.names:
            return False
        if not recarray.dtype[expected_name] == expected_dtype:
            return False
    return True
