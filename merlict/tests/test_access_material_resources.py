import merlict


def test_list_material_spectra():
    keys = merlict.materials.spectra.list_resources()
    assert len(keys) > 0


def test_list_material_surfaces():
    keys = merlict.materials.surfaces.list_resources()
    assert len(keys) > 0


def test_list_material_media():
    keys = merlict.materials.media.list_resources()
    assert len(keys) > 0


def test_read_material_surfaces():
    surf = merlict.materials.surfaces.init_from_resources(
        key="perfect_absorber"
    )
    assert "type" in surf
    assert surf["type"] == "cook-torrance"
    assert "reflection_spectrum" in surf
    assert "diffuse_weight" in surf
    assert "specular_weight" in surf
    assert "roughness" in surf


def test_read_material_media():
    medi = merlict.materials.media.init_from_resources(key="water_T293K")
    assert "refraction_spectrum" in medi
    assert "absorption_spectrum" in medi
