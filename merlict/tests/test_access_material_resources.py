import merlict


def test_list_material_surfaces():
    keys = merlict.materials.surfaces.list_resources()
    assert len(keys) > 0


def test_list_material_media():
    keys = merlict.materials.media.list_resources()
    assert len(keys) > 0


def test_read_material_surfaces():
    surf = merlict.materials.surfaces.init(key="perfect_absorber")
    assert "material" in surf
    assert "specular_reflection" in surf
    assert "diffuse_reflection" in surf
    assert "color" in surf


def test_read_material_media():
    medi = merlict.materials.media.init(key="water_T293K")
    assert "refraction" in medi
    assert "absorption" in medi
