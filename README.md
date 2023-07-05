![merlict python logo](/readme/merlict-python-logo-inkscape.png)

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

More light than you can handle! Also: This is in beta state. Don't judge me!


Highly techincal talking about sceneries
----------------------------------------

```sceneryPy```
---------------
A ```python``` ```dict```. Objects are ```dict```s with their vertices and faces in ```numpy.arrays``` or ```lists```. Geometric relations between objects and materials are also ```dicts```.

```sceneryDs```
---------------
A ```py```-dictionary where both the keys and the payload are ```str```-strings.
The objects are represented in wavefront-files (```.obj```). Materials are represented in ```JSON```-strings.
Lossless: [```sceneryTr```, ```sceneryUs```, ```sceneryDr```]

```sceneryTr```
---------------
A tape-archive (```.tar```). The filenames inside are the keys, and the files payloads are the values used in ```sceneryDs```.
Lossless: [```sceneryDs```, ```sceneryUs```, ```sceneryDr```]

```sceneryDr```
---------------
The same as ```sceneryTr``` but dumped into the filesystem, i.e. into directories and files according to the tape-archive.
Lossless: [```sceneryDs```, ```sceneryUs```, ```sceneryTr```]

```sceneryUs```
---------------
Merlict's in-memory-represenation used to build the oc-trees for acceleration.
It still has the hirachy defined in ```geometry/relations```.
Lossless: [```sceneryDs```, ```sceneryTr```]
Lossy: [```sceneryCa```]

```sceneryCa```
---------------
Merlict's in-memory-represenation used for the ray tracing. It is based on oc-trees. The object's faces are `flattened` into the scenery and then bundled in global oc-trees. This is good for cache-aware access. A seperate look-up-table is created to resolve the user's ```geometry/relations``` but is not needed when ray-face-intersections are queried.
Lossless: [```sceneryCs```]
Lossy: []

```sceneryCs```
---------------
Merlict can serialize and dump its ```sceneryCa``` into a binary blob.
Since the population of oc-trees in the conversion from ```sceneryUs``` to ```sceneryCa``` might take a significant time, ```sceneryCs``` is used to dump and load the scenery into caches.
This is only expected to work on a single machine. This is not meant to share sceneries and is expected to fail when read by machines with different architectures.
Lossless: [```sceneryCa```]
