import numpy as np
import io


def dumps(x, y, x_label, y_label):
    # RFC 4180
    assert len(x) == len(
        y
    ), f"len(x) {len(x):d} == {len(y):d} x:{x_label:s}, y:{y_label:s}."
    out = io.StringIO()
    assert "\n" not in x_label, "Expected x_label to be a single line."
    assert "\n" not in y_label, "Expected y_label to be a single line."
    out.write(f"{x_label:s},{y_label:s}\n")

    for i in range(len(x)):
        line = f"{x[i]:e},{y[i]:e}\n"
        out.write(line)
    out.seek(0)
    return out.read()


def loads(text):
    # RFC 4180
    x_label = ""
    y_label = ""
    x = []
    y = []
    is_header = True
    for line in text.splitlines():
        if len(line) == 0:
            break
        _x, _y = line.split(",")
        if is_header:
            x_label = str.strip(_x)
            y_label = str.strip(_y)
            is_header = False
        else:
            x.append(float(_x))
            y.append(float(_y))

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    return x, y, x_label, y_label
