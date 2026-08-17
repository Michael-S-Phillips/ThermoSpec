"""Build a ThermoSpec Mie optical-property endmember table from oriented mineral n,k.

Python replacement for the external Mie step (`compile_mie_results.py` parsed `.print` files from a
Fortran spheres code). Given complex refractive-index components n,k, this orientation-averages them
(for random-grain anisotropic minerals), resamples onto a target wavenumber grid, runs Mie theory
(miepython) for a single grain radius, and writes the 5-column table
`[lambda_um, g, Cext_um2, Csca_um2, ssalb]` (exactly the columns `rte_disort._load_constants` reads:
wns=1e4/col0, g=col1, Cext=col2 [um^2], Csca=col3 [um^2], ssalb=col4) plus a matching `_wn_bounds.txt`.

IMPORTANT (comparability): the model's particle volume Vp=(4/3)*pi*cfg.radius^3 and the optical depth
Et = n_p * Cext both use the SAME grain size, and (with cfg.scale_Et=False) the absolute Cext matters.
So `--radius` here MUST equal the `cfg.radius` used at run time, and to be comparable to the enstatite
endmember both should use the grain size enstatite was built at (its Cext magnitudes ~0.7-2.0e3 um^2 are
consistent with r ~= 14-15 um). Grain size is therefore a documented modelling choice, not recovered
from the enstatite table itself.

Orientation averaging: `eps` (default) averages the dielectric function eps_j=(n_j+i k_j)^2 over the
principal axes then m=sqrt(mean(eps)) -- the standard isotropic average for randomly oriented grains;
`nk` averages n and k directly. Either is an approximation (noted as a methods caveat).
"""
import argparse

import numpy as np
import miepython


def orientation_average(n_cols, k_cols, method):
    """n_cols,k_cols: [nwave, 3] principal-axis components -> (n, k) isotropic average [nwave]."""
    if method == "nk":
        return n_cols.mean(axis=1), k_cols.mean(axis=1)
    eps = (n_cols + 1j * k_cols) ** 2            # dielectric function per axis
    m = np.sqrt(eps.mean(axis=1))                # isotropic (Reuss-like) dielectric average
    return m.real, np.abs(m.imag)


def bin_average(wn_in, val_in, edges):
    """Bin-average val_in(wn_in) into the bins defined by `edges` (len nbin+1). Empty bins fall back
    to nearest input point. Matches Preprocessing/resample_optical_constants.py."""
    out = np.zeros(len(edges) - 1)
    lo = np.minimum(edges[:-1], edges[1:])
    hi = np.maximum(edges[:-1], edges[1:])
    for i in range(len(out)):
        m = (wn_in >= lo[i]) & (wn_in < hi[i])
        if np.any(m):
            out[i] = val_in[m].mean()
        else:
            out[i] = val_in[np.argmin(np.abs(wn_in - 0.5 * (lo[i] + hi[i])))]
    return out


def mie_table(lam_um, n, k, radius_um):
    """Per-band Mie: -> (g, Cext_um2, Csca_um2, ssalb). miepython uses m = n - i k (absorbing ->
    negative imaginary part); x = 2*pi*r/lambda; cross-section = efficiency * geometric area pi r^2."""
    m = n - 1j * k
    x = 2.0 * np.pi * radius_um / lam_um
    geom = np.pi * radius_um ** 2                 # um^2
    qext, qsca, qback, g = miepython.efficiencies_mx(m, x)   # miepython 3.x: (m, size parameter)
    Cext = qext * geom
    Csca = qsca * geom
    with np.errstate(divide="ignore", invalid="ignore"):
        ssalb = np.where(qext > 0, qsca / qext, 0.0)
    return g, Cext, Csca, ssalb


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nk", required=True, help="oriented n,k file: cols wn, lambda_um, n1,n2,n3, k1,k2,k3")
    ap.add_argument("--grid-bounds", required=True, help="target wavenumber bin-bounds file (e.g. enst_300K_wn_bounds.txt)")
    ap.add_argument("--out", required=True, help="output mie_combined .txt")
    ap.add_argument("--out-bounds", required=True, help="output wn_bounds .txt (copy of the target grid)")
    ap.add_argument("--radius-um", type=float, default=14.0, help="grain radius in um (MUST match cfg.radius at run time)")
    ap.add_argument("--orient", choices=["eps", "nk"], default="eps", help="orientation-average method")
    args = ap.parse_args()

    d = np.loadtxt(args.nk)
    wn_in = d[:, 0]
    n_cols, k_cols = d[:, 2:5], d[:, 5:8]
    n_iso, k_iso = orientation_average(n_cols, k_cols, args.orient)

    edges = np.loadtxt(args.grid_bounds)
    centers_wn = 0.5 * (edges[:-1] + edges[1:])
    n_g = bin_average(wn_in, n_iso, edges)
    k_g = bin_average(wn_in, k_iso, edges)
    lam_um = 1.0e4 / centers_wn

    g, Cext, Csca, ssalb = mie_table(lam_um, n_g, k_g, args.radius_um)

    table = np.column_stack([lam_um, g, Cext, Csca, ssalb])
    table = table[table[:, 0].argsort()]         # sort by lambda, as compile_mie_results.py does
    np.savetxt(args.out, table, fmt="%.6f", delimiter="\t")
    np.savetxt(args.out_bounds, edges, fmt="%.6e", delimiter="\t")

    print(f"wrote {args.out}  ({len(table)} bands, {table[:,0].min():.2f}-{table[:,0].max():.2f} um, "
          f"radius={args.radius_um} um, orient={args.orient})")
    print(f"  g {g.min():.3f}-{g.max():.3f}  Cext {Cext.min():.1f}-{Cext.max():.1f} um^2  "
          f"ssalb {ssalb.min():.3f}-{ssalb.max():.3f}")
    assert np.all((ssalb >= 0) & (ssalb <= 1.0000001)), "ssalb out of [0,1]"
    assert np.all(Cext >= Csca - 1e-9), "Cext < Csca"
    assert np.all(np.isfinite(table)), "non-finite entries"
    print("  sanity OK: ssalb in [0,1], Cext>=Csca, all finite")


if __name__ == "__main__":
    main()
