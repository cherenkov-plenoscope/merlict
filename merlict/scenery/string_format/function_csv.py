import numpy as np
import io


def dumps(func, header):
    # RFC 4180
    func = np.asarray(func)
    assert func.shape[1] == 2

    out = io.StringIO()
    if header is not None:
        assert len(header) == 2, "Expected header string for x and y."
        for head in header:
            assert "\n" not in head, "Expected header to be a single line."
        out.write(f"{header[0]:s},{header[1]:s}\n")

    for i in range(func.shape[0]):
        x = func[i, 0]
        y = func[i, 1]
        line = f"{x:e},{y:e}\n"
        out.write(line)
    out.seek(0)
    return out.read()


def loads(text, first_line_is_haeder=False):
    # RFC 4180
    header = None
    func = []
    found_numbers = False
    for line in text.splitlines():
        if len(line) == 0:
            break
        _x, _y = line.split(",")
        if first_line_is_haeder and header is None:
            header = (_x, _y)
        else:
            found_numbers = True
            func.append([float(_x), float(_y)])

    if found_numbers:
        func = np.asarray(func)
    else:
        func = np.zeros(shape=(0, 2), dtype=float)

    if header is None:
        return func
    else:
        return func, header
