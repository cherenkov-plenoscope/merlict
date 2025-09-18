from merlict.scenery.string_format.function_csv import loads
from merlict.scenery.string_format.function_csv import dumps
import numpy as np


def test_function_csv_empty():
    text = dumps(x=[], y=[], x_label="", y_label="")
    f = loads(text=text)
    assert len(f["x"]) == 0
    assert len(f["x"]) == len(f["y"])
    assert f["x_label"] == ""
    assert f["y_label"] == ""


def test_function_csv():
    func = np.asarray([[1, 2.5], [2, 3.5]])
    header = ("a/m", "b/m")

    text = dumps(
        x=func[:, 0], y=func[:, 1], x_label=header[0], y_label=header[1]
    )
    f = loads(text=text)
    header_back = f["x_label"], f["y_label"]
    func_back = np.asarray([f["x"], f["y"]]).T

    assert header == header_back
    np.testing.assert_array_equal(func, func_back)


def test_function_csv_ignores_trailing_newlines():
    ex = [1]
    ey = [2]
    exl = "a/m"
    eyl = "b/m"
    for data_lines in ["1,2", "1,2\n", "1,2\n\n"]:
        header_line = f"{exl:s},{eyl:s}\n"
        text = header_line + data_lines
        f = loads(text=text)
        assert (exl, eyl) == (f["x_label"], f["y_label"])
        np.testing.assert_array_equal(
            np.asarray([f["x"], f["y"]]), np.asarray([ex, ey])
        )
