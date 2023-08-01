import json_numpy


def sceneryDs_from_scenery(scenery, indent=4, relations_indent=0, readme=None):
    sceneryDs = {}

    default_readme = "Scenery\n"
    default_readme += "=======\n"
    sceneryDs["README.md"] = readme if readme else default_readme

    sceneryDs["geometry"] = {}

    sceneryDs["geometry"]["objects"] = {}
    for okey in uScenery["geometry"]["objects"]:
        okey_filename = "{:s}.obj".format(okey)
        wfr = obj.init_from_mesh(mesh=scenery["geometry"]["objects"][okey])
        sceneryDs["geometry"]["objects"][okey_filename] = obj.dumps(obj=wfr)

    sceneryDs["geometry"]["relations.json"] = json_numpy.dumps(
        scenery["geometry"]["relations"], indent=relations_indent
    )

    sceneryDs["materials"] = {}

    sceneryDs["materials"]["media"] = {}
    for mkey in scenery["materials"]["media"]:
        mkey_filename = "{:s}.json".format(mkey)
        mjsn = json_numpy.dumps(
            scenery["materials"]["media"][mkey], indent=indent
        )
        sceneryDs["materials"]["media"][mkey_filename] = mjsn

    sceneryDs["materials"]["surfaces"] = {}
    for skey in scenery["materials"]["surfaces"]:
        skey_filename = "{:s}.json".format(skey)
        sjsn = json_numpy.dumps(
            scenery["materials"]["surfaces"][skey], indent=indent
        )
        sceneryDs["materials"]["media"][mkey_filename] = sjsn

    sceneryDs["materials"]["boundary_layers.json"] = json_numpy.dumps(
        scenery["materials"]["boundary_layers"], indent=indent,
    )

    sceneryDs["materials"]["default_medium.txt"] = str(
        scenery["materials"]["default_medium"]
    )

    return sceneryDs
