"""Facet-to-facet view-factor generator for terrain / crater meshes.

The repo can only *read* view factors (crater.SelfHeatingList); this computes them for any mesh
with the CraterMesh interface (normals, areas, centroids, vertices, faces), with line-of-sight
occlusion, and writes them in that reader's sparse format so the crater radiative engine is
reused unchanged. First of the real-topography sub-projects; see
docs/superpowers/specs/2026-08-14-view-factor-generator-design.md.

Point-to-point view factor (facets as small patches at their centroids):
    F_ij = (cos_i * cos_j * A_j) / (pi * r^2)
with cos_i = n_i . u_ij, cos_j = n_j . (-u_ij), u_ij the unit c_i->c_j, r = |c_j - c_i|; zero
unless both cosines are positive and the ray i->j is unoccluded. The occlusion mask is
symmetrized so F satisfies reciprocity A_i F_ij = A_j F_ji exactly.

Memory note: builds dense [N,N] arrays -- fine for the crater/synthetic meshes here (N ~ 1e2-1e3).
Large DEM meshes (N ~ 1e4) will want a chunked/sparse pass; deferred to the DEM-loader sub-project.
"""
import numpy as np


def compute_view_factors(mesh, occlusion=True, lift=None):
    """Dense [N,N] view-factor matrix F[i,j] (i->j). Diagonal is 0."""
    normals = np.asarray(mesh.normals, float)
    centroids = np.asarray(mesh.centroids, float)
    areas = np.asarray(mesh.areas, float)
    N = len(normals)

    Fgeom = np.zeros((N, N))
    vis = np.ones((N, N), dtype=bool)          # ray visibility (True = not occluded)

    tm = None
    if occlusion:
        import trimesh
        tm = trimesh.Trimesh(vertices=np.asarray(mesh.vertices, float),
                             faces=np.asarray(mesh.faces), process=False)
        if lift is None:
            lift = 1e-3 * float(np.median(tm.edges_unique_length))

    for i in range(N):
        D = centroids - centroids[i]           # [N,3]  c_j - c_i
        r = np.linalg.norm(D, axis=1)          # [N]
        r_safe = np.where(r > 0, r, 1.0)
        U = D / r_safe[:, None]                # unit i->j
        cos_i = U @ normals[i]                 # n_i . u_ij
        cos_j = -(U * normals).sum(axis=1)     # n_j . (-u_ij)
        facing = (cos_i > 0) & (cos_j > 0) & (r > 0)
        Fgeom[i] = np.where(facing, cos_i * cos_j * areas / (np.pi * r_safe**2), 0.0)

        if occlusion and facing.any():
            js = np.where(facing)[0]
            origins = np.repeat((centroids[i] + lift * normals[i])[None, :], len(js), axis=0)
            first_hit = tm.ray.intersects_first(origins, U[js])
            occluded = first_hit != js         # blocked (or missed -> conservative)
            vis[i, js[occluded]] = False

    # Symmetric visibility -> exact reciprocity (Fgeom is reciprocal by construction).
    F = Fgeom * (vis & vis.T)
    np.fill_diagonal(F, 0.0)
    return F


def write_selfheating_list(F, fname, threshold=1e-8):
    """Write F in the sparse SelfHeatingList format: per facet, `n idx1..idxn vf1..vfn`
    with 1-based indices (SelfHeatingList subtracts 1 on read)."""
    N = F.shape[0]
    with open(fname, "w") as fh:
        for i in range(N):
            js = np.where(F[i] > threshold)[0]
            parts = [str(len(js))]
            parts += [str(int(j) + 1) for j in js]
            parts += ["%.12g" % F[i, j] for j in js]
            fh.write(" ".join(parts) + "\n")


def view_factor_matrix_from_file(fname, N):
    """Load a sparse view-factor file back into a dense [N,N] matrix (via SelfHeatingList)."""
    from crater import SelfHeatingList
    return SelfHeatingList(fname).as_view_matrix(N)
