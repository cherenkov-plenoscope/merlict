import merlict


def test_version():
    v = merlict.__version__
    assert len(v.split(".")) == 3 + 3
