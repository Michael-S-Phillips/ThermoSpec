"""DEM -> triangular mesh for real-topography thermal modeling (terrain sub-project 2).

`DEMMesh` subclasses `crater.CraterMesh` so it exposes the identical attribute set
(`normals, areas, centroids, vertices, faces, sub_*, tangent1/2`) and drops straight into
`ShadowTester`, `CraterRadiativeTransfer`, and `view_factors.compute_view_factors`. Only the
mesh construction is new: a heightfield E[nrows,ncols] on a regular (dx,dy) grid becomes two
triangles per cell in a local Cartesian frame (x=col*dx, y=row*dy, z=elevation), wound so facet
normals point up (+z) -- valid for any single-valued heightfield (no overhangs).

`load_dem` reads a DEM from GeoTIFF (rasterio, lazy), ESRI ASCII (.asc), or NumPy (.npy, needs an
explicit scale) and returns (elevation, dx, dy).
See the analysis handoff spec (psr_topography_implementation_spec.md, section 2.1).
"""
import numpy as np

from crater import CraterMesh


class DEMMesh(CraterMesh):
    """A CraterMesh built from a DEM heightfield instead of a mesh file."""

    def __init__(self, elevation, dx, dy=None, origin='centroid'):
        E = np.asarray(elevation, dtype=float)
        if E.ndim != 2:
            raise ValueError("elevation must be a 2D array [nrows, ncols]")
        self.nrows, self.ncols = E.shape
        self.dx = float(dx)
        self.dy = float(dx if dy is None else dy)
        self.elevation = E

        self.vertices, self.faces = self._build_from_dem(E)
        if origin == 'centroid':
            self.vertices = self.vertices - self.vertices.mean(axis=0)

        # reuse CraterMesh geometry machinery unchanged
        self.sub_vertices, self.sub_faces, self.sub_face_index = self.subdivide()
        (self.normals, self.areas, self.centroids,
         self.sub_normals, self.sub_areas, self.sub_centroids) = self.compute_geometry()
        self.tangent1, self.tangent2 = self._compute_local_coordinates()

    def _build_from_dem(self, E):
        nr, nc = E.shape
        jj, ii = np.meshgrid(np.arange(nc), np.arange(nr))
        vertices = np.column_stack([jj.ravel() * self.dx,
                                    ii.ravel() * self.dy,
                                    E.ravel()])

        def vid(i, j):
            return i * nc + j

        faces = []
        for i in range(nr - 1):
            for j in range(nc - 1):
                a, b = vid(i, j), vid(i, j + 1)
                c, d = vid(i + 1, j), vid(i + 1, j + 1)
                # winding chosen so a flat heightfield yields +z normals
                faces.append([a, b, c])
                faces.append([b, d, c])
        return vertices, np.array(faces, dtype=np.int64)


def load_dem(path, dx=None, dy=None):
    """Load a DEM -> (elevation[nrows,ncols], dx, dy). Supports .tif/.tiff (rasterio),
    .asc (ESRI ASCII grid), and .npy (requires explicit dx)."""
    ext = str(path).lower().rsplit('.', 1)[-1]
    if ext in ('tif', 'tiff'):
        try:
            import rasterio
        except ImportError as e:
            raise ImportError("Reading GeoTIFF DEMs needs rasterio: pip install rasterio") from e
        with rasterio.open(path) as ds:
            E = ds.read(1).astype(float)
            tr = ds.transform
            return E, abs(tr.a), abs(tr.e)
    if ext == 'asc':
        return _load_ascii_grid(path)
    if ext == 'npy':
        E = np.load(path).astype(float)
        if dx is None:
            raise ValueError(".npy DEMs carry no scale; pass dx (and dy) explicitly")
        return E, float(dx), float(dx if dy is None else dy)
    raise ValueError(f"Unsupported DEM format: .{ext} (use .tif/.tiff, .asc, or .npy)")


def _load_ascii_grid(path):
    header = {}
    with open(path) as fh:
        lines = fh.readlines()
    data_start = 0
    for k, line in enumerate(lines):
        parts = line.split()
        if parts and parts[0].lower() in (
                'ncols', 'nrows', 'xllcorner', 'yllcorner', 'xllcenter', 'yllcenter',
                'cellsize', 'nodata_value'):
            header[parts[0].lower()] = float(parts[1])
            data_start = k + 1
        else:
            break
    nrows = int(header['nrows'])
    ncols = int(header['ncols'])
    cell = float(header.get('cellsize', 1.0))
    vals = np.array([float(x) for line in lines[data_start:] for x in line.split()],
                    dtype=float)
    E = vals[:nrows * ncols].reshape(nrows, ncols)
    nodata = header.get('nodata_value')
    if nodata is not None:
        E = np.where(E == nodata, np.nan, E)
    return E, cell, cell
