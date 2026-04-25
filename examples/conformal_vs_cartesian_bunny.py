"""Side-by-side: Cartesian :class:`Infill` vs :class:`Conformal` on a bunny.

Both fill the same bunny envelope with the same gyroid TPMS at the same
offset and ``repeat_cell``.  Two visualizations:

- **Left** — ``Infill``: the gyroid is in the world Cartesian frame, the bunny
  surface just *clips* it.  Cells are aligned to global axes; near the
  envelope you see them sliced on arbitrary planes.
- **Right** — ``Conformal``: the gyroid's *radial* coord is the signed
  distance to the bunny surface.  Cells stack as concentric "shells" along
  the surface normal — an onion-skin scaffold that always presents
  perpendicular to the envelope.

Run with ``python examples/conformal_vs_cartesian_bunny.py``.
"""

from __future__ import annotations

import numpy as np
import pyvista as pv
from pyvista import examples

from microgen.shape.surface_functions import gyroid
from microgen.shape.tpms import Conformal, Infill


# -----------------------------------------------------------------------------
# Build the envelope
# -----------------------------------------------------------------------------

bunny = examples.download_bunny()
# Original bunny is ~6cm; scale to a reasonable working size and recenter.
bunny.transform(np.diag([40.0, 40.0, 40.0, 1.0]), inplace=True)
bunny.translate(-np.array(bunny.center_of_mass()), inplace=True)
print(f"bunny: volume={bunny.volume:.3f}  bounds={tuple(np.round(bunny.bounds, 2))}")


# -----------------------------------------------------------------------------
# Two TPMS infills with identical parameters
# -----------------------------------------------------------------------------

OFFSET = 0.4
REPEAT_CELL = 5

cartesian_infill = Infill(
    obj=bunny,
    surface_function=gyroid,
    offset=OFFSET,
    repeat_cell=REPEAT_CELL,
)
cartesian_sheet = cartesian_infill.generate_vtk(type_part="sheet")
print(
    f"Cartesian Infill : sheet_volume={abs(cartesian_sheet.volume):.3f}  "
    f"density_vs_bunny={abs(cartesian_sheet.volume) / bunny.volume:.2%}",
)

conformal = Conformal(
    envelope=bunny,
    surface_function=gyroid,
    offset=OFFSET,
    repeat_cell=REPEAT_CELL,
    default_tangent_axis=(1.0, 0.0, 0.0),
)
conformal_sheet = conformal.generate_vtk(type_part="sheet")
print(
    f"Conformal        : sheet_volume={abs(conformal_sheet.volume):.3f}  "
    f"density_vs_bunny={abs(conformal_sheet.volume) / bunny.volume:.2%}",
)


# -----------------------------------------------------------------------------
# Cut both meshes through y for a clean cross-section view
# -----------------------------------------------------------------------------

cut_origin = (0.0, -0.3 * float(bunny.bounds[3] - bunny.bounds[2]), 0.0)
cartesian_view = cartesian_sheet.clip("y", origin=cut_origin, invert=False)
conformal_view = conformal_sheet.clip("y", origin=cut_origin, invert=False)


# -----------------------------------------------------------------------------
# Side-by-side plot
# -----------------------------------------------------------------------------

plotter = pv.Plotter(shape=(1, 2), window_size=(1600, 800))

plotter.subplot(0, 0)
plotter.add_text("Cartesian Infill\n(gyroid in world frame, bunny clips)", font_size=10)
plotter.add_mesh(bunny, color="white", opacity=0.15)
plotter.add_mesh(cartesian_view, color="orange")
plotter.view_isometric()

plotter.subplot(0, 1)
plotter.add_text(
    "Conformal\n(gyroid in surface-normal frame — onion shells)",
    font_size=10,
)
plotter.add_mesh(bunny, color="white", opacity=0.15)
plotter.add_mesh(conformal_view, color="cyan")
plotter.view_isometric()

plotter.link_views()
plotter.show()
