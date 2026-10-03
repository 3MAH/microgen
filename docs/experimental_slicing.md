Experimental meshless slicing
=============================

`microgen.slicing` supplies shared layer data, a contour-source interface, streaming implicit raster sampling, and adapters for UVtools and a caller-configured PySLM hatcher. All coordinates use millimeters in one explicit build frame. This interface is experimental and may change.

Install the optional dependencies:

```console
pip install -e ".[slicing]"
```

A field-only solid requires no TPMS grid or triangle mesh:

```python
import numpy as np
from microgen import surface_functions
from microgen.slicing import BoundedSolid, ImplicitContourSource, LayerPlan, slice_contours

k = 2 * np.pi / 2.0
def sheet(x, y, z):
    return np.abs(surface_functions.gyroid(k*x, k*y, k*z)) - 0.3

solid = BoundedSolid(sheet, bounds=(0, 8, 0, 8, 0, 8))
plan = LayerPlan.uniform(z_min=0, z_max=8, thickness=0.05)
source = ImplicitContourSource(solid)

for layer in slice_contours(source, plan, tolerance=0.025):
    print(layer.position.z_top, len(layer.regions))
```

The sheet's `0.3` is a level-set offset, not a wall thickness in millimeters. Gradient-normalized fields are distance estimates. Slicing itself needs a negative-inside field, not an exact SDF.

`BoundedSolid` intersects the field with its local clipping box. Its optional `build_transform` is a 4x4 affine matrix mapping local coordinates into build coordinates. `BoundedSolid.from_shape(shape)` uses the shape's field and bounds. It does not infer the renderer's center/orientation. Choose the solid phase and build placement explicitly.

`LayerPlan` records bottom, top, and sample Z separately. Uniform plans sample at each slab's midpoint and clip the final slab to `z_max`. Arbitrary callables remain supported. No expression graph or OEM kernel integration is required.

Contours use a uniform XY grid, marching squares, and bracketed edge-root refinement. `tolerance` limits grid pitch, not certified boundary error. Equal-sign corners can hide an interior feature, and saddle cases follow the algorithm's low-valued connectivity convention. Every result reports that topology is uncertified. The default point cap rejects oversized planes; increase tolerance or explicitly raise the cap. Adaptive contour extraction is future work.

Each `ContourRegion` has a closed counterclockwise exterior and explicit clockwise holes, validated as a polygon. Empty layers are retained. A custom source can implement `contours_at(z, *, tolerance)` and return the same `ContourLayer` data. Mesh and BREP sources can be added at that seam when required.

To feed a real PySLM hatcher:

```console
pip install PythonSLM==0.6.1
```

```python
from pyslm.hatching import Hatcher
from microgen.slicing import hatch_with_pyslm

hatcher = Hatcher()
# Configure your process settings on hatcher before use.
for position, native_layer in hatch_with_pyslm(
    slice_contours(source, plan, tolerance=0.025), hatcher,
):
    # Empty layers yield None. Assign machine/build-style/model IDs downstream.
    print(position.z_top, native_layer)
```

The adapter passes writable copies of ordinary polygon paths to the existing hatcher. It does not require PySLM to understand Microgen geometry. Use the accompanying `LayerSpec` to assign the native layer's Z and index in your receiver's units, then assign model and build-style IDs. Mixed-part layers are rejected; hatch each part separately to preserve process identity. Build files and machine control remain outside this module.

For raster output, `PixelGrid` specifies dimensions, X/Y pitch, and a lower-left origin. Columns increase along build X. Image row zero corresponds to highest build Y. Samples lie at pixel centers. Supersampling averages binary occupancy into uint8 coverage values. Grayscale coverage is not an exposure calibration.

Export a UVJ from UVtools using a known printer/material profile, then use it as the template:

```python
from microgen.slicing import PixelGrid, read_uvj_profile, slice_rasters, write_uvj

profile = read_uvj_profile("printer-profile.uvj")
size = profile["Properties"]["Size"]
grid = PixelGrid(
    size["X"], size["Y"],
    (size["Millimeter"]["X"] / size["X"],
     size["Millimeter"]["Y"] / size["Y"]),
)
height = size["LayerHeight"]
plan = LayerPlan.uniform(z_min=0, z_max=np.ceil(8 / height) * height, thickness=height)
layers = slice_rasters(solid, plan, grid=grid, supersampling=2)
write_uvj(layers, "gyroid.uvj", profile=profile)
```

The writer requires full-display images at origin zero, matching profile dimensions and constant layer height, with contiguous slabs starting at zero. A shorter final slab is rejected. Choose a total plan height that is a whole number of printer layers. The bounded solid can end earlier than that plan.

Global normal and bottom exposure/motion settings are copied from the profile. Source per-layer overrides are not copied. Original pixel images and geometry are replaced. Output is written atomically after validation, with zero-based `slice/00000000.png` entries, `config.json`, and top-of-slab Z positions. Image stacks stream one layer at a time; only per-layer metadata accumulates in memory.

UVJ follows the existing [UVtools implementation](https://github.com/sn4k3/UVtools/blob/e892d3f0c466ce795f0ffffd7365fae9c64d6b78/UVtools.Core/FileFormats/UVJFile.cs). Open the result in UVtools and convert using the selected printer's existing format. Confirm orientation with an asymmetric calibration coupon. UVtools conversion does not supply conventional 3D supports; include any necessary support geometry before slicing.

Run the example with `python examples/ImplicitSlicing/gyroid_layers.py`. Pass `--profile printer-profile.uvj --output gyroid.uvj` for a raster stack. Large display/layer counts can take time; use the shared computing queue for heavy runs.

Validation commands:

```console
pytest tests/test_slicing.py -q
```

Set `UVTOOLS_CMD` to the UVtools command-line executable to enable `tests/test_uvtools_integration.py`. That test asks real UVtools to normalize a test UVJ, converts generated layers to GOO, then decodes them back into UVJ and compares pixels and metadata. The focused suite passed with UVtools 7.0.1 and PythonSLM 0.6.1 on Windows, including the native GOO round-trip and hatching an annulus without crossing its hole. Test process values are synthetic and are not printer recommendations. Physical printing and commercial LPBF interoperability remain unverified.

Validation on 2026-10-03 produced 25 passing focused tests. The full suite produced 432 passes, 62 skips, and 22 failures. Rerunning all 22 failing cases on the parent commit `38e1912` in the same environment reproduced the same test names and failure messages. These failures concern existing TPMS fixtures and CAD operations. The full suite is therefore not green in this environment, although the new slicing tests pass.
