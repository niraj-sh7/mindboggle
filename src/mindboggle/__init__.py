import os

import numpy as np

from .version import __version__

# from .info import (LONG_DESCRIPTION as __doc__,
#                  __version__)
# __doc__ += """
# """

# Set up package information function
try:
    from .pkg_info import get_pkg_info as _get_pkg_info
except:
    get_info = lambda: ""
else:
    get_info = lambda: _get_pkg_info(os.path.dirname(__file__))

# module imports
# from . import blah as blah
# object imports
# from .blah import blah, blah

INIT_MSG = "Running {packname} version {version} (latest: {latest})".format
latest = {"version": "Unknown"}

# NumPy 2 removed several legacy aliases that this codebase and its doctests
# still exercise. Restore the aliases in place so older doctests keep working
# without rewriting hundreds of examples.
for _name, _value in (
    ("float", float),
    ("int", int),
):
    if _name not in np.__dict__:
        setattr(np, _name, _value)

try:
    from nibabel.spatialimages import SpatialImage as _SpatialImage

    def _get_data(self):
        return np.asanyarray(self.dataobj)

    def _get_affine(self):
        return self.affine

    def _get_header(self):
        return self.header

    _SpatialImage.get_data = _get_data
    _SpatialImage.get_affine = _get_affine
    _SpatialImage.get_header = _get_header
except Exception:
    pass

try:
    import etelemetry

    latest = etelemetry.get_project("nipy/mindboggle")
except Exception as e:
    print("Could not check for version updates: ", e)
finally:
    print(
        INIT_MSG(packname="mindboggle", version=__version__, latest=latest["version"])
    )


if not hasattr(np, "asscalar"):
    np.asscalar = lambda a: a.item()
