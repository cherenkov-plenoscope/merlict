from ... import utils


def assert_spectrum_is_valid(spectrum, ymin=None, ymax=None):
    assert len(spectrum["x"]) >= 2
    assert len(spectrum["x"]) == len(spectrum["y"])

    assert utils.is_free_of_nan(spectrum["x"])
    assert utils.is_all_greater_zero(spectrum["x"])
    assert utils.is_monotonically_increasing(spectrum["x"])

    assert utils.is_free_of_nan(spectrum["y"])

    if ymin is not None:
        assert np.all(spectrum["y"] >= ymin)

    if ymax is not None:
        assert np.all(spectrum["y"] <= ymax)

    assert "\n" not in spectrum["x_label"]
    assert "\n" not in spectrum["y_label"]
