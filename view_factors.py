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

Occlusion backends (`occlusion_backend=`): 'numpy' (default, dependency-free, O(N^3) full scan --
fine to ~1e3 facets), 'numba' (grid-DDA line-walk, bit-identical but ~O(neighbours) per ray -- 80-115x
faster, makes nx>=40 DEM meshes tractable; needs numba), 'auto' (numba if available else numpy), or
'trimesh'. Memory note: builds dense [N,N] arrays -- fine to N ~ few 1e3; N ~ 1e4 still wants a
chunked/sparse pass.
"""
import warnings
import numpy as np

try:
    import numba
    _HAS_NUMBA = True
except Exception:                                            # numba optional -> numpy fallback
    _HAS_NUMBA = False


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


def _build_face_xy_grid(tris, cell):
    """Bin faces into a uniform XY grid by their bounding box -> CSR (cell_start, cell_faces).
    A face is listed in every cell its XY AABB covers."""
    fxmin = tris[:, :, 0].min(1); fxmax = tris[:, :, 0].max(1)
    fymin = tris[:, :, 1].min(1); fymax = tris[:, :, 1].max(1)
    ox = float(fxmin.min()); oy = float(fymin.min())
    ncx = int((float(fxmax.max()) - ox) / cell) + 1
    ncy = int((float(fymax.max()) - oy) / cell) + 1
    ixmn = np.clip(((fxmin - ox) / cell).astype(np.int64), 0, ncx - 1)
    ixmx = np.clip(((fxmax - ox) / cell).astype(np.int64), 0, ncx - 1)
    iymn = np.clip(((fymin - oy) / cell).astype(np.int64), 0, ncy - 1)
    iymx = np.clip(((fymax - oy) / cell).astype(np.int64), 0, ncy - 1)
    counts = np.zeros(ncx * ncy, np.int64)
    for f in range(len(tris)):
        for ix in range(ixmn[f], ixmx[f] + 1):
            base = ix * ncy
            for iy in range(iymn[f], iymx[f] + 1):
                counts[base + iy] += 1
    cell_start = np.zeros(ncx * ncy + 1, np.int64); cell_start[1:] = np.cumsum(counts)
    cell_faces = np.empty(int(cell_start[-1]), np.int64)
    cursor = cell_start[:-1].copy()
    for f in range(len(tris)):
        for ix in range(ixmn[f], ixmx[f] + 1):
            base = ix * ncy
            for iy in range(iymn[f], iymx[f] + 1):
                c = base + iy
                cell_faces[cursor[c]] = f; cursor[c] += 1
    return cell_start, cell_faces, ox, oy, ncx, ncy


if _HAS_NUMBA:
    @numba.njit(cache=True)
    def _occ_kernel(v0, e1, e2, origins, dirs, seglen, src, dst,
                    cell_start, cell_faces, ox, oy, cell, ncx, ncy, eps):
        """Per-ray occlusion via a 2D grid line-walk (Amanatides-Woo DDA): only faces in the cells
        the segment's XY projection crosses are Moller-Trumbore tested. Result is identical to the
        full scan -- a triangle can only intersect the segment inside a cell the segment crosses that
        the triangle's AABB also covers (and it is binned into all such cells)."""
        M = src.shape[0]
        blocked = np.zeros(M, dtype=np.bool_)
        for k in range(M):
            ox0 = origins[k, 0]; oy0 = origins[k, 1]; oz0 = origins[k, 2]
            dx = dirs[k, 0]; dy = dirs[k, 1]; dz = dirs[k, 2]; L = seglen[k]
            ex = ox0 + dx * L; ey = oy0 + dy * L
            s_f = src[k]; d_f = dst[k]
            cx = int((ox0 - ox) / cell); cy = int((oy0 - oy) / cell)
            cxe = int((ex - ox) / cell); cye = int((ey - oy) / cell)
            ddx = ex - ox0; ddy = ey - oy0
            stepx = 1 if ddx >= 0 else -1
            stepy = 1 if ddy >= 0 else -1
            if ddx != 0.0:
                tMaxX = (ox + (cx + (1 if stepx > 0 else 0)) * cell - ox0) / ddx
                tDeltaX = cell / abs(ddx)
            else:
                tMaxX = 1e30; tDeltaX = 1e30
            if ddy != 0.0:
                tMaxY = (oy + (cy + (1 if stepy > 0 else 0)) * cell - oy0) / ddy
                tDeltaY = cell / abs(ddy)
            else:
                tMaxY = 1e30; tDeltaY = 1e30
            hit = False
            for _step in range(ncx + ncy + 2):
                if 0 <= cx < ncx and 0 <= cy < ncy:
                    c = cx * ncy + cy
                    for idx in range(cell_start[c], cell_start[c + 1]):
                        f = cell_faces[idx]
                        if f == s_f or f == d_f:
                            continue
                        e1x = e1[f, 0]; e1y = e1[f, 1]; e1z = e1[f, 2]
                        e2x = e2[f, 0]; e2y = e2[f, 1]; e2z = e2[f, 2]
                        px = dy * e2z - dz * e2y; py = dz * e2x - dx * e2z; pz = dx * e2y - dy * e2x
                        det = e1x * px + e1y * py + e1z * pz
                        if -eps < det < eps:
                            continue
                        inv = 1.0 / det
                        tx = ox0 - v0[f, 0]; ty = oy0 - v0[f, 1]; tz = oz0 - v0[f, 2]
                        u = (tx * px + ty * py + tz * pz) * inv
                        if u < -eps or u > 1 + eps:
                            continue
                        qx = ty * e1z - tz * e1y; qy = tz * e1x - tx * e1z; qz = tx * e1y - ty * e1x
                        vv = (dx * qx + dy * qy + dz * qz) * inv
                        if vv < -eps or u + vv > 1 + eps:
                            continue
                        tt = (e2x * qx + e2y * qy + e2z * qz) * inv
                        if tt > eps and tt < L - 1e-3:
                            hit = True
                            break
                if hit or (cx == cxe and cy == cye):
                    break
                if tMaxX < tMaxY:
                    tMaxX += tDeltaX; cx += stepx
                else:
                    tMaxY += tDeltaY; cy += stepy
            blocked[k] = hit
        return blocked


def _occluded_pairs_numba(vertices, faces, centroids, normals, src, dst, lift_frac=1e-3, eps=1e-6):
    """numba grid-DDA occluder: bit-identical to _occluded_pairs_numpy but O(neighbours) per ray
    instead of O(N) -- makes occlusion tractable for large DEM meshes (nx>=40)."""
    tris = vertices[faces]
    v0 = np.ascontiguousarray(tris[:, 0, :])
    e1 = np.ascontiguousarray(tris[:, 1, :] - tris[:, 0, :])
    e2 = np.ascontiguousarray(tris[:, 2, :] - tris[:, 0, :])
    med_edge = np.median(np.linalg.norm(np.diff(tris[:, [0, 1, 2, 0], :], axis=1), axis=2))
    lift = lift_frac * med_edge
    origins = np.ascontiguousarray(centroids[src] + lift * normals[src])
    seg = centroids[dst] - origins
    seglen = np.linalg.norm(seg, axis=1)
    dirs = np.ascontiguousarray(seg / seglen[:, None])
    cell = max(2.0 * float(med_edge), 1e-12)
    cell_start, cell_faces, ox, oy, ncx, ncy = _build_face_xy_grid(tris, cell)
    return _occ_kernel(v0, e1, e2, origins, dirs, seglen,
                       np.asarray(src, np.int64), np.asarray(dst, np.int64),
                       cell_start, cell_faces, float(ox), float(oy), float(cell),
                       int(ncx), int(ncy), float(eps))


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
        backend = occlusion_backend
        if backend == 'auto':                            # numba if available (big DEM meshes), else numpy
            backend = 'numba' if _HAS_NUMBA else 'numpy'
        if backend in ('numpy', 'numba'):
            src, dst = np.nonzero(Fgeom)                 # only facing pairs need testing
            fn = _occluded_pairs_numba if backend == 'numba' else _occluded_pairs_numpy
            if backend == 'numba' and not _HAS_NUMBA:
                raise RuntimeError("occlusion_backend='numba' requires numba (pip install numba)")
            blk = fn(verts, faces, centroids, normals, src, dst,
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

    # --- physical safety guard for the far-field point kernel (CS audit 2026-08-27) ---
    # F = cos_i cos_j A_j/(pi r^2) is far-field only; for adjacent facets (r/sqrt(A) <~ 0.5) it can give
    # F_ij > 1 and row sums > 1 -- unphysical, and it makes the radiosity Neumann series diverge
    # (spectral radius > 1) for bright/low-emissivity surfaces. Cap each pair at F_ij <= 1 in the
    # reciprocity-preserving symmetric quantity M_ij = A_i F_ij = A_j F_ji (cap at min(A_i,A_j), which
    # is symmetric so reciprocity survives). This is a NO-OP on well-resolved meshes (production PSR
    # row-sum <= 0.09; refine=True also avoids it). A remaining row-sum > 1 signals the kernel is
    # inadequate for this mesh -> warn (use refine=True or a near-field contour integral). We do NOT
    # row-normalize: that would break reciprocity / energy conservation.
    M = areas[:, None] * F
    cap = np.minimum(areas[:, None], areas[None, :])
    n_over = int(np.count_nonzero(M > cap + 1e-12))
    if n_over:
        M = np.minimum(M, cap)
        F = M / areas[:, None]
        np.fill_diagonal(F, 0.0)
        warnings.warn(f"view_factors: capped {n_over} facet pair(s) with F_ij>1 (far-field kernel "
                      f"breakdown for adjacent facets); consider refine=True or a near-field integral",
                      RuntimeWarning)
    rmax = float(F.sum(1).max())
    if rmax > 1.0 + 1e-9:
        warnings.warn(f"view_factors: max row-sum {rmax:.3f} > 1 -> radiosity spectral radius may exceed "
                      f"1 (divergence for bright surfaces); mesh too steep/coarse for the far-field "
                      f"kernel -- use refine=True or a near-field contour integral", RuntimeWarning)
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
