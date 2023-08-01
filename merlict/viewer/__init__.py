import numpy as np


def init_config():
    return {
        "random_seed": 42,
        "preview_num_cols": 160,
        "preview_num_rows": 90 // 2,
        "export_num_cols": 640,
        "export_num_rows": 360,
        "step_length": 0.1,
        "view": {
            "position": [0, 0, 0],
            "rotation": [np.deg2rad(90), 0, 0],
            "field_of_view": np.deg2rad(80.0),
        },
        "aperture_camera_f_stop_ratio": 2.0,
        "aperture_camera_image_sensor_width": 24e-3,
    }
