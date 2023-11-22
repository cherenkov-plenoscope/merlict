import posixpath


def list():
    return [
        "README.md",
        posixpath.join("geometry", "objects", "*.obj"),
        posixpath.join("geometry", "relations.json"),
        posixpath.join("materials", "media", "*.json"),
        posixpath.join("materials", "surfaces", "*.json"),
        posixpath.join("materials", "boundary_layers.json"),
        posixpath.join("materials", "default_medium.txt"),
    ]
