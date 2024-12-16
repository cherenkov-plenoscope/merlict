from merlict.scenery.string_format.function_csv import loads
from merlict.scenery.string_format.function_csv import dumps
import numpy as np


def test_function_csv_empty():
    text = dumps(func=np.zeros(shape=(0, 2)), header=None)
    func_back = loads(text=text, first_line_is_haeder=False)
    assert func_back.shape == (0, 2)


def test_function_csv():
    func = np.asarray([[1, 2.5], [2, 3.5]])
    header = ("a/m", "b/m")

    text = dumps(func=func, header=header)
    func_back, header_back = loads(text=text, first_line_is_haeder=True)

    assert header == header_back
    np.testing.assert_array_equal(func, func_back)


def test_function_csv_ignores_trailing_newlines():
    expected = np.asarray([[1, 2]])
    for text in ["1,2", "1,2\n", "1,2\n\n"]:
        func = loads(text=text, first_line_is_haeder=False)
        np.testing.assert_array_equal(func, expected)

    header = ("a/m", "b/m")
    for data_lines in ["1,2", "1,2\n", "1,2\n\n"]:
        header_line = f"{header[0]:s},{header[1]}\n"
        text = header_line + data_lines
        func, head = loads(text=text, first_line_is_haeder=True)
        assert head == header
        np.testing.assert_array_equal(func, expected)
