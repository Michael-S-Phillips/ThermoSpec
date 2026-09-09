"""1D layered cold-trap thermal model with annual + lunation forcing.

Purpose: predict the buried-ice surface signature dT_B(depth) for a PERMANENTLY
SHADOWED crater floor under seasonal forcing, so the 3D seasonal probe has a
falsifiable target rather than an unconstrained number.

Physics (matches the 3D model's cold-trap configuration):
  - heat equation  rho*cp dT/dt = d/dz( k(z) dT/dz )
  - k(z): k_dust above the ice table, k_ice below (step contrast)
  - surface: eps*sigma*T^4 = Q_wall(t) + conduction from below
             Q_wall(t) = Q0 + dQ_ann*sin(2 pi t/P_ann) + dQ_lun*sin(2 pi t/P_lun)
             (a shadowed floor's only radiative input is wall IR; no direct beam)
  - base: fixed geothermal flux F_geo (Neumann), matching the 3D bottom BC
  - IC: equilibrium conductive profile dT/dz = F_geo/k(z)  (the eqic CC implemented)

Numerics: backward Euler, harmonic-mean interface conductivity, Newton on the
nonlinear radiative surface term. Geometric grid resolves the 2.2 cm lunation
skin depth at the top and reaches 5 m at the base.
"""
import numpy as np

SIGMA = 5.670374419e-8
P_LUN = 29.530588 * 86400.0
P_ANN = 365.25 * 86400.0
F_GEO = 0.018
EPS = 1.0


def make_grid(z_max=5.0, n=220, dz0=2.0e-3):
    """Geometric grid: cell faces at 0 = f_0 < ... < f_n = z_max."""
    # bisect for the growth ratio r solving dz0*(r^n - 1)/(r - 1) = z_max
    def total(rr):
        return dz0*n if abs(rr-1.0) < 1e-12 else dz0*(rr**n - 1)/(rr - 1)
    lo, hi = 1.0, 1.2
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if total(mid) < z_max: lo = mid
        else: hi = mid
    r = 0.5*(lo+hi)
    dz = dz0 * r**np.arange(n)
    faces = np.concatenate([[0.0], np.cumsum(dz)])
    centres = 0.5*(faces[:-1] + faces[1:])
    return faces, centres, dz


def k_profile(centres, ice_depth_m, k_dust=5.5e-4, k_ice=2.0):
    k = np.full(len(centres), k_dust)
    if ice_depth_m is not None:
        k[centres >= ice_depth_m] = k_ice
    return k


def rhocp_profile(centres, ice_depth_m, rho_d=1100.0, cp_d=825.0,
                  rho_i=920.0, cp_i=800.0):
    rc = np.full(len(centres), rho_d*cp_d)
    if ice_depth_m is not None:
        rc[centres >= ice_depth_m] = rho_i*cp_i
    return rc


def equilibrium_ic(centres, k, T_surf0):
    """dT/dz = F_geo/k integrated from the surface: the eqic profile."""
    dz = np.gradient(centres)
    return T_surf0 + np.cumsum(F_GEO/k * dz)


def run(ice_depth_m, dQ_ann, Q0=0.148, dQ_lun=0.0, years=6.0,
        steps_per_lun=240, z_max=5.0, n=220, T_surf0=41.0, k_ice=2.0,
        return_series=False):
    """Integrate to periodic steady state. Returns surface-T diagnostics."""
    faces, centres, dz = make_grid(z_max, n)
    k = k_profile(centres, ice_depth_m, k_ice=k_ice)
    rc = rhocp_profile(centres, ice_depth_m)
    T = equilibrium_ic(centres, k, T_surf0)

    # interface conductances k_iface/d between adjacent cell centres (harmonic mean)
    d_c = np.diff(centres)
    k_if = 2.0*k[:-1]*k[1:]/(k[:-1] + k[1:])
    G = k_if/d_c                                  # W/m2/K between cells i and i+1

    dt = P_LUN/steps_per_lun
    nsteps = int(np.ceil(years*P_ANN/dt))
    N = len(centres)

    Ts_hist = np.empty(nsteps); t_hist = np.empty(nsteps)
    for s in range(nsteps):
        t = (s+1)*dt
        Q = Q0 + dQ_ann*np.sin(2*np.pi*t/P_ANN) + dQ_lun*np.sin(2*np.pi*t/P_LUN)
        Told = T.copy()
        # Newton on the surface node (radiative nonlinearity); interior linear
        for _ in range(40):
            A = np.zeros((3, N))                   # tridiagonal bands
            b = np.zeros(N)
            cap = rc*dz/dt
            A[1] = cap
            b = cap*Told
            # interior conduction
            A[1, :-1] += G;  A[2, :-1] = -G        # upper diag (i, i+1)
            A[1, 1:]  += G;  A[0, 1:]  = -G        # lower diag (i, i-1)
            # surface radiative flux, linearized: eps sig T^4 ~ 4 eps sig Ts^3 T - 3 eps sig Ts^4
            Ts = T[0]
            A[1, 0] += 4*EPS*SIGMA*Ts**3
            b[0] += Q + 3*EPS*SIGMA*Ts**4
            # base geothermal influx
            b[-1] += F_GEO
            ab = np.zeros((3, N))
            ab[0, 1:] = A[0, 1:]; ab[1] = A[1]; ab[2, :-1] = A[2, :-1]
            from scipy.linalg import solve_banded
            Tnew = solve_banded((1, 1), ab, b)
            if abs(Tnew[0] - T[0]) < 1e-10:
                T = Tnew; break
            T = Tnew
        Ts_hist[s] = T[0]; t_hist[s] = t

    # last full year = periodic steady state
    last = t_hist >= (years-1.0)*P_ANN
    out = dict(T_mean=float(Ts_hist[last].mean()), T_min=float(Ts_hist[last].min()),
               T_max=float(Ts_hist[last].max()),
               T_amp=float(0.5*(Ts_hist[last].max()-Ts_hist[last].min())),
               T_final_profile=T, centres=centres, k=k)
    if return_series:
        out["t"] = t_hist; out["Ts"] = Ts_hist; out["last"] = last
    return out
