from . import fileorder
import tarfile
import io
import fnmatch
import posixpath
import warnings


def write(sceneryStr, path):
    with tarfile.open(name=path, mode="w|") as tar:
        for item in sceneryStr:
            filepath, payload = item
            _tar_append_str(tar=tar, filepath=filepath, payload=payload)


def read(path):
    sceneryDS = []
    with tarfile.open(name=path, mode="r") as tar:
        while True:
            tar_info = tar.next()

            if tar_info is None:
                break

            if tar_info.type != tarfile.REGTYPE:
                continue

            if _is_name_in_expected_name_patterns(name=tar_info.name):
                payload = _tar_read_str(tar=tar, tar_info=tar_info)
                item = (posixpath.normpath(tar_info.name), payload)
                sceneryDS.append(item)
            else:
                warnings.warn(
                    f"Ignoring '{tar_info.name:s}' while reading tarfile."
                )

    return sceneryDS


def _is_name_in_expected_name_patterns(name):
    norm_name = posixpath.normpath(name)
    for filepath in fileorder.list():
        out = fnmatch.filter(names=[norm_name], pat=filepath)
        if len(out) > 0:
            return True
    return False


def _tar_append_str(tar, filepath, payload):
    with io.BytesIO() as buff:
        tar_info = tarfile.TarInfo(filepath)
        tar_info.size = buff.write(str.encode(payload))
        buff.seek(0)
        tar.addfile(tar_info, buff)


def _tar_read_str(tar, tar_info):
    return bytes.decode(tar.extractfile(tar_info).read())
