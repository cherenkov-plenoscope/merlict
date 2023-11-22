from . import fileorder
import tarfile
import io
import fnmatch


def write(sceneryStr, path):
    with tarfile.open(name=path, mode="w|") as tar:
        for item in sceneryStr:
            filepath, payload = item
            _tar_append_str(tar=tar, filepath=filepath, payload=payload)


def read(path):
    sceneryDS = []
    with tarfile.open(name=path, mode="r") as tar:
        filepaths = [member.name for member in tar.getmembers()]
        for filepath in fileorder.list():
            if "*" in filepath:
                for ifilepath in fnmatch.filter(filepath):
                    payload = _tar_read_str(tar=tar, filepath=ifilepath)
                    item = (ifilepath, payload)
                    sceneryDS.append(item)
            else:
                payload = _tar_read_str(tar=tar, filepath=filepath)
                item = (filepath, payload)
                sceneryDS.append(item)
    return sceneryDS


def _tar_append_str(tar, filepath, payload):
    with io.BytesIO() as buff:
        info = tarfile.TarInfo(filepath)
        info.size = buff.write(str.encode(payload))
        buff.seek(0)
        tar.addfile(info, buff)


def _tar_read_str(tar, filepath):
    return bytes.decode(tar.extractfile(filepath).read())
