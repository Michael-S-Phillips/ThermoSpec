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


def _point_vf_geometry(centroids, normals, areas):
    """Dense point-to-point geometric view factors Fgeom[a,b] = cos_a cos_b A_b/(pi r^2),
    zeroed where facets do not face each other. Reciprocal by construction (A_a Fgeom_ab =
    A_b Fgeom_ba). No occlusion."""
    n = len(centroids)
    Fg = np.zeros((n, n))
    for a in range(n):
        D = centroids - centroids[a]
        r = np.linalg.norm(D, axis=1)
        r_safe = np.where(r > 0, r, 1.0)
        U = D / r_safe[:, None]
        cos_a = U @ normals[a]
        cos_b = -(U * normals).sum(axis=1)
        facing = (cos_a > 0) & (cos_b > 0) & (r > 0)
        Fg[a] = np.where(facing, cos_a * cos_b * areas / (np.pi * r_safe**2), 0.0)
    return Fg


def _occluded_pairs_numpy(vertices, faces, centroids, normals, src, dst, lift_frac=1e-3, eps=1e-6):
    """Pure-numpy Moller-Trumbore occlusion: blocked[k] for each pair (src[k]->dst[k]), True if a
    third facet's triangle intersects the segment strictly between the two centroids. Dependency-
    free alternative to the trimesh ray index (rtree/libspatialindex). Contributed by the analysis
    session (claude_session_sync/scripts/numpy_occlusion.py); vendored here so the generator has no
    ray-engine dependency."""
    tris = vertices[faces]                                    # (N,3,3)
    v0 = tris[:, 0, :]; e1 = tris[:, 1, :] - v0; e2 = tris[:, 2, :] - v0
    med_edge = np.median(np.linalg.norm(np.diff(tris[:, [0, 1, 2, 0], :], axis=1), axis=2))
    lift = lift_frac * med_edge
    origins = centroids[src] + lift * normals[src]
    seg = centroids[dst] - origins
    seglen = np.linalg.norm(seg, axis=1)
    dirs = seg / seglen[:, None]
    blocked = np.zeros(len(src), dtype=bool)
    for k in range(len(src)):
        o, d, L = origins[k], dirs[k], seglen[k]
        p = np.cross(d, e2)
        det = np.einsum('ij,ij->i', e1, p)
        ok = np.abs(det) > eps
        inv = np.zeros(len(faces)); inv[ok] = 1.0 / det[ok]
        t0 = o - v0
        u = np.einsum('ij,ij->i', t0, p) * inv
        q = np.cross(t0, e1)
        v = np.einsum('j,ij->i', d, q) * inv
        t = np.einsum('ij,ij->i', e2, q) * inv
        hit = ok & (u >= -eps) & (v >= -eps) & (u + v <= 1 + eps) & (t > eps) & (t < L - 1e-3)
        hit[src[k]] = False; hit[dst[k]] = False
        blocked[k] = hit.any()
    return blocked


def compute_view_factors(mesh, occlusion=True, refine=False, lift=None, occlusion_backend='numpy'):
    """Dense [N,N] view-factor matrix F[i,j] (i->j), diagonal 0.

    refine=True integrates over the mesh's subdivided sub-facets (mesh.sub_* + sub_face_index)
    for accuracy when facets are closely spaced (r ~ facet size); the point approximation on raw
    facets over-counts near-field view factors. Occlusion is evaluated at the facet level.
    """
    normals = np.asarray(mesh.normals, float)
    centroids = np.asarray(mesh.centroids, float)
    areas = np.asarray(mesh.areas, float)
    N = len(normals)

    # --- geometric view factors (point, or sub-facet integrated) ---
    if refine:
        sc = np.asarray(mesh.sub_centroids, float)
        sn = np.asarray(mesh.sub_normals, float)
        sa = np.asarray(mesh.sub_areas, float)
        Fsub = _point_vf_geometry(sc, sn, sa)          # [Nsub, Nsub]
        smap = mesh.sub_face_index                      # facet -> sub-facet indices
        Fgeom = np.zeros((N, N))
        for i in range(N):
            ai = np.asarray(smap[i])
            # F[i,j] = (1/A_i) sum_{a in i} A_a sum_{b in j} Fsub[a,b]
            contrib = (sa[ai][:, None] * Fsub[ai]).sum(axis=0) / areas[i]   # per sub-facet b
            for j in range(N):
                Fgeom[i, j] = contrib[np.asarray(smap[j])].sum()
        np.fill_diagonal(Fgeom, 0.0)
    else:
        Fgeom = _point_vf_geometry(centroids, normals, areas)

    # --- facet-level occlusion ---
    vis = np.ones((N, N), dtype=bool)
    if occlusion:
        verts = np.asarray(mesh.vertices, float)
        faces = np.asarray(mesh.faces)
        if occlusion_backend == 'numpy':
            src, dst = np.nonzero(Fgeom)                 # only facing pairs need testing
            blk = _occluded_pairs_numpy(verts, faces, centroids, normals, src, dst,
                                        lift_frac=(1e-3 if lift is None else lift))
            vis[src[blk], dst[blk]] = False
        elif occlusion_backend == 'trimesh':
            import trimesh
            tm = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
            lft = 1e-3 * float(np.median(tm.edges_unique_length)) if lift is None else lift
            D = centroids[None, :, :] - centroids[:, None, :]
            r = np.linalg.norm(D, axis=2)
            for i in range(N):
                js = np.where(Fgeom[i] > 0)[0]
                if len(js) == 0:
                    continue
                dirs = D[i, js] / r[i, js][:, None]
                origins = np.repeat((centroids[i] + lft * normals[i])[None, :], len(js), axis=0)
                first_hit = tm.ray.intersects_first(origins, dirs)
                vis[i, js[first_hit != js]] = False
        else:
            raise ValueError(f"Unknown occlusion_backend: {occlusion_backend!r} (use 'numpy' or 'trimesh')")

    F = Fgeom * (vis & vis.T)                            # symmetric mask -> exact reciprocity
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


class ViewFactorList:
    """SelfHeatingList-compatible view of a dense view-factor matrix, for injecting generated
    view factors straight into the crater engine (CraterRadiativeTransfer) without a temp file.
    Duck-types crater.SelfHeatingList: exposes `indices`, `view_factors`, and `as_view_matrix`."""

    def __init__(self, F, threshold=1e-8):
        F = np.asarray(F, float)
        self._F = F
        N = F.shape[0]
        idx, vfs = [], []
        for i in range(N):
            js = np.where(F[i] > threshold)[0]
            idx.append(js)
            vfs.append(F[i, js])
        self.indices = np.array(idx, dtype=object)
        self.view_factors = np.array(vfs, dtype=object)

    def as_view_matrix(self, N):
        return self._F
