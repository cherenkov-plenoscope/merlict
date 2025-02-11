from .wrapper cimport *

cimport numpy as cnumpy
cnumpy.import_array()

from libc cimport stdint
from libc cimport stdlib
from libc cimport string

import termios
import sys
import numpy as np
from .. import ray as _ray
from .. import photon as _photon
from .. import intersection as _intersection
from .. import intersectionSurfaceNormal as _intersectionSurfaceNormal


cdef _chk_msg(rc, msg):
    assert rc == CHK_SUCCESS, msg


cdef _mli_Vec2py(mli_Vec mliv):
    return np.array([mliv.x, mliv.y, mliv.z], dtype=np.float64)


cdef _mli_Vec(v):
    cdef mli_Vec mliv
    mliv.x = v[0]
    mliv.y = v[1]
    mliv.z = v[2]
    return mliv


cdef _mli_Image2py(mli_Image mliimg):
    out = np.zeros(
        shape=(mli_Image_num_cols(&mliimg), mli_Image_num_rows(&mliimg), 3),
        dtype=np.float32)
    cdef mli_Color c
    for ix in range(mli_Image_num_cols(&mliimg)):
        for iy in range(mli_Image_num_rows(&mliimg)):
            c = mli_Image_get_by_col_row(&mliimg, ix, iy)
            out[ix, iy, 0] = c.r
            out[ix, iy, 1] = c.g
            out[ix, iy, 2] = c.b
    return out


cdef _mli_View(position, rotation, field_of_view):
    cdef mli_View view
    view.position = _mli_Vec(position)
    view.rotation = _mli_Vec(rotation)
    view.field_of_view = field_of_view
    return view


cdef _mli_viewer_Config_init(
    step_length,
    preview_num_cols,
    preview_num_rows,
    export_num_cols,
    export_num_rows,
    view_position,
    view_rotation_tait_bryan_xyz,
    field_of_view,
    aperture_camera_f_stop_ratio,
    aperture_camera_image_sensor_width,
    random_seed,
    gain,
    gamma
):
    assert step_length > 0
    assert preview_num_cols > 0
    assert preview_num_rows > 0
    assert export_num_cols > 0
    assert export_num_rows > 0
    assert 0 < field_of_view <= np.pi
    assert aperture_camera_f_stop_ratio > 0
    assert aperture_camera_image_sensor_width > 0
    assert gamma > 0.0

    cdef mli_viewer_Config _c
    _c.random_seed = int(random_seed)
    _c.preview_num_cols = int(preview_num_cols)
    _c.preview_num_rows = int(preview_num_rows)
    _c.export_num_cols = int(export_num_cols)
    _c.export_num_rows = int(export_num_rows)
    _c.step_length = float(step_length)
    _c.view = _mli_View(
        position=view_position,
        rotation=view_rotation_tait_bryan_xyz,
        field_of_view=float(field_of_view),
    )
    _c.gain = float(gain)
    _c.gamma = float(gamma)
    _c.aperture_camera_f_stop_ratio = float(aperture_camera_f_stop_ratio)
    _c.aperture_camera_image_sensor_width = float(
        aperture_camera_image_sensor_width
    )
    return _c


cdef _mli_Archive_push_back_path_and_payload(
    mli_Archive *archive,
    path,
    payload,
):
    cdef int rc

    _path = str(path)
    _payload = str(payload)

    cdef bytes _py_path = _path.encode()
    cdef bytes _py_payload = _payload.encode()

    cdef stdint.uint64_t path_length = np.uint64(len(_py_path))
    cdef stdint.uint64_t payload_length = len(_py_payload)

    cdef char* _cpath = _py_path
    cdef char* _cpayload = _py_payload

    _chk_msg(mli_Archive_push_back_cstr(
        archive,
        _cpath,
        path_length,
        _cpayload,
        payload_length
    ), "Failed to push back cstr into mli_Archive.")
    return


cdef bytes _make_sure_bytes(b):
    if type(b) is bytes:
        # Fast path for most common case(s).
        return <bytes>b
    else:
        raise TypeError("Input must be bytes.")


cdef class Merlict:
    """
    A scenery of objects inside a tree structure to accelerate ray-intersection
    queries, ray tracing and path tracing.

    functions
    ---------
    init_from_sceneryStr
        Init from a scenery represented in sceneryStr (list of srings).
        A new tree structure for acceleration will be build from scratch.
    view
        An interactive viewer that renders the scenery in real-time and
        prints the images to stdout.
    query_intersection
        The intersection of a ray (or many rays) with the scenery.
    query_intersectionSurfaceNormal
        Like 'query_intersection' but with additional information about the
        surface-normal.
    dump
        Serialize to a path. Dumps include the tree structure for acceleration.
    init_from_dump
        Init from a dump without setting up the tree structure for
        acceleration again from just a sceneryStr.
    """
    cdef mli_Scenery scenery

    def __cinit__(self):
        self.scenery = mli_Scenery_init()

    def __dealloc__(self):
        mli_Scenery_free(&self.scenery)

    def __init__(self, path=None, sceneryStr=None):
        if path and not sceneryStr:
            try:
                self.init_from_tar(path)
            except AssertionError:
                try:
                    self.init_from_dump(path)
                except AssertionError:
                    raise AssertionError(
                        "Can not read scenery from path {:s}".format(path)
                    )

        elif sceneryStr and not path:
            self.init_from_sceneryStr(sceneryStr)
        else:
            raise ValueError("Either 'path' or 'sceneryStr', but not both.")

    def view(
        self,
        step_length=0.1,
        preview_num_cols=160,
        preview_num_rows=45,
        export_num_cols=1280,
        export_num_rows=720,
        view_position=np.array([0.0, 0.0, 0.0]),
        view_rotation_tait_bryan_xyz=np.array([np.deg2rad(90.0), 0.0, 0.0]),
        field_of_view=np.deg2rad(80.0),
        aperture_camera_f_stop_ratio=2.0,
        aperture_camera_image_sensor_width=24e-3,
        random_seed=0,
        gain=0.1,
        gamma=1.0,
    ):
        """
        An interactive view which displays images in stdout.
        The viewer waits for user-input (the user pressing keys) to
        manipulate the view. At each moment, a high resolution image can
        be rendered which will be written to the current working directory.

        Press [h] to print the help.

        parameters
        ----------
        step_length : float
            How far the view-port moves in each step.
        preview_num_cols : int
            Number of columns in the image printed to stdout.
        preview_num_rows : int
            Number of rows in the image printed to stdout. The pixels on stdout
            are ssumed to be twice as high as they are wide. For this, to
            obtain a quadratic image, num-rows must be approx. 1/2 num-cols.
        view_position : [float, float, float]
            Initial position of the view-port.
        view_rotation_tait_bryan_xyz : [float, float, float]
            Initial orientation of the view-port. In units of rad.
        field_of_view : float
            Initial field-of-view of the view-port. In units of rad.
        random_seed : int
            For the path-tracer.
        export_num_cols : int
            Number of columns in high-res image.
        export_num_rows : int
            Number of rows in high-res image. The high-res image's pixel are
            square.
        aperture_camera_f_stop_ratio : float
            F-stop of the camera rendering the high-res image.
        aperture_camera_image_sensor_width : float
            Physical width (along the columns) of the image sensor in the
            camera rendering the high-res image.
        gain : float
            Color and image gain.
        gamma : float
            Color values in the output image are raised to the power of gamma
            in order to compensate for different illumination levels.
        """
        cdef mli_viewer_Config cconfig

        cconfig = _mli_viewer_Config_init(
            step_length=step_length,
            preview_num_cols=preview_num_cols,
            preview_num_rows=preview_num_rows,
            export_num_cols=export_num_cols,
            export_num_rows=export_num_rows,
            view_position=view_position,
            view_rotation_tait_bryan_xyz=view_rotation_tait_bryan_xyz,
            field_of_view=field_of_view,
            aperture_camera_f_stop_ratio=aperture_camera_f_stop_ratio,
            aperture_camera_image_sensor_width=aperture_camera_image_sensor_width,
            random_seed=random_seed,
            gain=gain,
            gamma=gamma)

        fd = sys.stdin.fileno()
        old_attr = termios.tcgetattr(fd)
        new_attr = termios.tcgetattr(fd)
        C_FLAG = 3

        new_attr[C_FLAG] = new_attr[C_FLAG] & ~termios.ICANON
        try:
            termios.tcsetattr(fd, termios.TCSADRAIN, new_attr)
            _rc = mli_viewer_run_interactive_viewer(&self.scenery, cconfig)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_attr)

    def init_from_tar(self, path):
        cdef int rc
        _path = str(path)
        cdef bytes _py_path = _path.encode()
        cdef char* _cpath = _py_path  # Who is responsible for this memory?

        cdef mli_Archive archive = mli_Archive_init()
        try:
            _chk_msg(mli_Archive__from_path_cstr(&archive, _cpath),
                "Failed to read mli_Archive from path.")
            _chk_msg(mli_Scenery_malloc_from_Archive(&self.scenery, &archive),
                "Failed to make mli_Scenery from mli_Archive.")
        finally:
            mli_Archive_free(&archive)

    def init_from_dump(self, path):
        """
        Loads Merlict from a previous dump.

        Warning
        -------
        A dump is not meant to exchange sceneries with others!
        Only read dumps written on the same platform and by the same version
        of merlict. See also dump().

        Parameters
        ----------
        path : str
            Path to read the dump from.
        """
        cdef int rc
        _path = str(path)
        cdef bytes _py_path = _path.encode()
        cdef char* _cpath = _py_path
        _chk_msg(mli_Scenery_malloc_from_path(&self.scenery, _cpath),
            "Failed to malloc scenery from dump in path.")

    def dump(self, path):
        """
        Dumps the compiled Merlict to path.

        Warning
        -------
        A dump is not meant to exchange sceneries with others!
        A dump is specific to your platform's architecture and version of
        merlict and will probably crash on other machines.
        Think of it as a `pickle` thingy.
        The dump is meant to export and import a compiled scenery so that
        merlict does not need to compile the user's scenery again.
        So the objective of the dump is for local caching only!

        Parameters
        ----------
        path : str
            Path to write the dump to.
        """
        cdef int rc
        _path = str(path)
        cdef bytes _py_path = _path.encode()
        cdef char* _cpath = _py_path
        _chk_msg(mli_Scenery_write_to_path(&self.scenery, _cpath),
            "Failed to dump mli_Scenery into path.")

    def dumps(self):
        """
        Returns a binary dump of the merlict scenery (the c89 binary dump).
        This can be loaded again with loads(dump).

        Do not use merlict_c89 binary dumps to share sceneries with others
        as you would also never share pickle.dumps with others.
        Only use your own merlict_c89 dumps.
        """
        cdef int rc
        cdef mli_IO buff = mli_IO_init()
        cdef bytes pyout

        try:
            _chk_msg(mli_IO_open_memory(&buff), "Failed to open buffer.")
            assert buff.type == MLI_IO_TYPE_MEMORY
            assert buff.data.memory.size == 0

            _chk_msg(mli_Scenery_to_io(&self.scenery, &buff),
                "Failed to serialize scenery to buffer.")

            pyout = buff.data.memory.cstr[:buff.data.memory.size]

        finally:
            _chk_msg(mli_IO_close(&buff),
                "Failed to close buffer.")

        return pyout

    def loads(self, dump):
        """
        Loads a merlict_c89 binary dump of the scenery.
        The binary dump can be dumped using dumps().

        Do not use merlict_c89 binary dumps to share sceneries with others
        as you would also never share pickle.dumps with others.
        Only use your own merlict_c89 dumps.
        """
        cdef bytes _dump = _make_sure_bytes(dump)
        cdef mli_IO buff = mli_IO_init()
        cdef stdint.uint64_t size = len(_dump)

        try:
            _chk_msg(mli_IO_open_memory(&buff), "Failed to open buffer.")
            assert buff.type == MLI_IO_TYPE_MEMORY
            assert buff.data.memory.size == 0

            _chk_msg(mli_IoMemory__malloc_capacity(&buff.data.memory, size),
                "Failed to malloc buffer from dump.")

            assert buff.data.memory.capacity == size
            assert buff.data.memory.pos == 0
            buff.data.memory.size = size

            for i in range(size):
                buff.data.memory.cstr[i] = _dump[i]

            # The copy operations below fail. No idea what is going on there.
            # string.memcpy(buff.data.memory.cstr, &_dump[0], size)
            # buff.data.memory.cstr[0:size] = _dump[0:size]

            _chk_msg(mli_Scenery_from_io(&self.scenery, &buff),
                "Failed to load scenery from buffer.")
        finally:
            _chk_msg(mli_IO_close(&buff), "Failed to close buffer.")

    def init_from_sceneryStr(self, sceneryStr):
        cdef mli_Archive tmp_archive = mli_Archive_init()
        try:
            _chk_msg(mli_Archive_malloc(&tmp_archive),
                "Failed to malloc mli_Archive")

            for item in sceneryStr:
                filename, payload = item
                _mli_Archive_push_back_path_and_payload(
                    &tmp_archive,
                    filename,
                    payload)

            _chk_msg(mli_Scenery_malloc_from_Archive(
                    &self.scenery,
                    &tmp_archive),
                "Failed to malloc mli_Scenery from mli_Archive.")
        finally:
            mli_Archive_free(&tmp_archive)

    def query_intersection(self, rays):
        assert _ray.israys(rays)
        cdef stdint.uint64_t num_ray = rays.shape[0]
        isecs = _intersection.init(size=num_ray)

        cdef cnumpy.ndarray[mli_Ray, mode="c"] crays = np.ascontiguousarray(
            rays
        )

        cdef cnumpy.ndarray[
            mli_Intersection, mode="c"
        ] cisecs = np.ascontiguousarray(
            isecs
        )

        cdef cnumpy.ndarray[
            stdint.int64_t, mode="c"
        ] cis_valid_isecs = np.ascontiguousarray(
            np.zeros(rays.shape[0], dtype=np.int64)
        )

        if num_ray:
            mli_Bridge_query_many_intersection(
                &self.scenery,
                num_ray,
                &crays[0],
                &cisecs[0],
                &cis_valid_isecs[0])

        isecs_mask = cis_valid_isecs.astype(dtype=np.bool_)

        return isecs_mask, isecs

    def query_intersectionSurfaceNormal(self, rays):
        assert _ray.israys(rays)
        cdef stdint.uint64_t num_ray = rays.shape[0]
        isecs = _intersectionSurfaceNormal.init(size=num_ray)

        cdef cnumpy.ndarray[mli_Ray, mode="c"] crays = np.ascontiguousarray(
            rays
        )

        cdef cnumpy.ndarray[
            mli_IntersectionSurfaceNormal, mode="c"
        ] cisecs = np.ascontiguousarray(
            isecs
        )

        cdef cnumpy.ndarray[
            stdint.int64_t, mode="c"
        ] cis_valid_isecs = np.ascontiguousarray(
            np.zeros(rays.shape[0], dtype=np.int64)
        )

        if num_ray:
            mli_Bridge_query_many_intersectionSurfaceNormal(
                &self.scenery,
                num_ray,
                &crays[0],
                &cisecs[0],
                &cis_valid_isecs[0])

        isecs_mask = cis_valid_isecs.astype(dtype=np.bool_)

        return isecs_mask, isecs

    def __repr__(self):
        out = "{:s}()".format(self.__class__.__name__)
        return out
