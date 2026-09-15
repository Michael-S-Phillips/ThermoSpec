#!/usr/bin/env python3
"""Science acceptance gates for new ThermoSpec-3D PSR runs — standalone.

Any new run dropped into the sync data tree is tested against the physics
established this week. Each gate is a hard pass/fail with a printed number, so a
returning PI can see whether a weekend run is trustworthy WITHOUT re-deriving
the diagnostics. Run:  python3 tools/check_science_gates.py [file.npz ...]
With no arguments it scans the sync tree for prod_*/seasonal_* floor files.

GATES
  G1 shadow integrity   Floor facets blocked at every azimuth must not exceed the
                        temperature their wall-IR view supports. Catches the summer
                        direct-beam leak (floor 198 K where wall IR supports ~73 K).
  G2 IC drainage        Conductive flux through the dust cap must approach F_geo.
                        Catches the initial-condition artifact (cap was passing 20x).
  G3 monotone column    A converged cold-trap column rises monotonically with depth.
                        Catches the undrained interior maximum (109 K at node 159).
  G4 forcing regime     A seasonal run must actually excite the annual wave: report
                        the sun-elevation span and whether it crosses the horizon.
Exit 0 only if every applicable gate passes on every file.
"""
import glob, json, os, sys
import numpy as np

SIGMA = 5.670374419e-8
F_GEO = 0.018
SYNC = "/Users/phillipsm/Documents/Research/Publications/artemis-thermal-modeling/claude_session_sync"

def az_min_horizon(cen, idx, naz=36):
    """Minimum-over-azimuth horizon angle (deg) within the mesh, for facets idx."""
    edges = np.linspace(-180, 180, naz + 1)
    out = np.empty(len(idx))
    for k, i in enumerate(idx):
        p = cen[i]
        dx, dy = cen[:, 0] - p[0], cen[:, 1] - p[1]
        dxy = np.hypot(dx, dy)
        ok = dxy > 1.0
        az = np.degrees(np.arctan2(dx[ok], dy[ok]))
        el = np.degrees(np.arctan((cen[ok, 2] - p[2]) / dxy[ok]))
        hz = np.full(naz, -90.0)
        b = np.clip(np.digitize(az, edges) - 1, 0, naz - 1)
        np.maximum.at(hz, b, el)
        out[k] = hz.min()
    return out

def gate_shadow(d):
    """G1: blocked facets must be radiatively consistent with wall IR only."""
    if not {"Tsurf", "centroids", "sunelev_out", "areas"} <= set(d.files):
        return None, "missing keys"
    cen = np.asarray(d["centroids"]); T = d["Tsurf"]; se = d["sunelev_out"]
    if se.max() <= 0:
        return True, f"sun never above horizontal (max {se.max():+.2f} deg) — gate not applicable"
    ti = int(np.argmax(T.mean(0)))
    lowz = np.where(cen[:, 2] < np.percentile(cen[:, 2], 35))[0]     # candidate floor
    hz = az_min_horizon(cen, lowz)
    blocked = lowz[hz > se.max()]
    if blocked.size == 0:
        return True, "no fully-blocked low facets in mesh"
    A = np.asarray(d["areas"]); other = np.setdiff1d(np.arange(len(cen)), blocked)
    Tw = ((A[other] * T[other, ti] ** 4).sum() / A[other].sum()) ** 0.25
    F = np.asarray(d["F_rowsum"])[blocked].mean() if "F_rowsum" in d.files else 0.06
    # A shadowed floor is heated by wall IR AND by the geothermal flux from below.
    # Omitting F_geo makes an equilibrium floor look like a violation (false positive
    # caught on a synthetic equilibrium-IC run: 41 K flagged against a 20 K bound).
    T_allowed = ((F * SIGMA * Tw ** 4 + F_GEO) / SIGMA) ** 0.25
    Tb = T[blocked, ti].mean()
    ok = Tb <= T_allowed * 1.5 + 5.0
    return ok, (f"blocked floor facets n={blocked.size}: mean {Tb:.1f} K vs allowed "
                f"{T_allowed:.1f} K from wall IR + F_geo (eff wall {Tw:.1f} K, F={F:.4f})")

def gate_ic(d, bf=None):
    """G2: dust-cap conductive flux must approach F_geo."""
    if "T_crater_out" not in d.files:
        return None, "no T_crater_out"
    if "depth_m" in d.files and "k_profile" in d.files:
        z, k = np.asarray(d["depth_m"]), np.asarray(d["k_profile"])
    elif bf is not None:
        z, k = bf
    else:
        return None, "no depth_m/k_profile (and no backfill supplied)"
    prof = d["T_crater_out"].mean(axis=(1, 2))
    n = min(len(prof), len(z), len(k))
    prof, z, k = prof[:n], z[:n], k[:n]
    good = z > 1e-6
    q = (k * np.gradient(prof, z))[good] * 1e3
    cap = q[: max(3, int(0.3 * good.sum()))]
    ratio = float(np.mean(np.abs(cap)) / (F_GEO * 1e3))
    return ratio < 3.0, f"mean |cap flux| = {np.mean(np.abs(cap)):.1f} mW/m2 = {ratio:.1f}x F_geo (need <3x)"

def gate_monotone(d):
    """G3: converged cold-trap column rises monotonically with depth."""
    if "T_crater_out" not in d.files:
        return None, "no T_crater_out"
    prof = d["T_crater_out"].mean(axis=(1, 2))
    i = int(np.argmax(prof))
    frac = i / (len(prof) - 1)
    return frac > 0.95, (f"profile peak at node {i}/{len(prof)-1} (frac {frac:.2f}); "
                         f"T {prof[0]:.1f} -> peak {prof.max():.1f} -> base {prof[-1]:.1f} K")

def gate_forcing(d, fname=""):
    """G4: does the run excite the annual wave?

    Only meaningful for runs that CLAIM to be seasonal. Legacy single-epoch
    production runs are 30-day tiled-lunation by design, so flagging them here
    is noise, not a finding — skip unless the filename marks it seasonal."""
    if "sunelev_out" not in d.files:
        return None, "no sunelev_out"
    if not any(tag in os.path.basename(fname).lower()
               for tag in ("seasonal", "annual", "multiyear", "eqic")):
        se = d["sunelev_out"]
        return None, (f"single-epoch run by design (sun {se.min():+.2f}..{se.max():+.2f} deg); "
                      f"gate applies only to seasonal runs")
    se = d["sunelev_out"]
    span = float(se.max() - se.min())
    crosses = bool(se.min() < 0 < se.max())
    t = np.asarray(d["t_out"]).ravel() if "t_out" in d.files else None
    dur_d = (t.max() - t.min()) / 86400.0 if t is not None and t.size > 1 else float("nan")
    ok = crosses and dur_d > 300
    return ok, (f"sun elev {se.min():+.2f}..{se.max():+.2f} deg (span {span:.2f}), "
                f"crosses horizon={crosses}, duration {dur_d:.0f} d (need >300 d for annual)")

def _floor_mask(d):
    if "polygon_floor" in d.files and np.asarray(d["polygon_floor"]).sum() > 0:
        return np.asarray(d["polygon_floor"], bool), "polygon"
    if "floor_elev_p20" in d.files:
        return np.asarray(d["floor_elev_p20"], bool), "elev-p20"
    return None, None

def _per_lunation_means(d):
    """Time-weighted floor-mean BT per lunation index -> {L: mean}. Needs a convergence npz."""
    if "T_surf_crater_history" not in d.files or "t_history" not in d.files:
        return None, None
    mask, tag = _floor_mask(d)
    if mask is None:
        return None, None
    H = np.asarray(d["T_surf_crater_history"]); t = np.asarray(d["t_history"], float)
    fm = H[mask].mean(axis=0)
    P = 2551443.0
    lun = ((t - t.min()) / P).astype(int)
    means = {}
    for L in np.unique(lun):
        sel = lun == L
        if sel.sum() < 2:
            continue
        ts, fs = t[sel], fm[sel]
        if ts[-1] > ts[0]:
            means[L] = float(np.sum(0.5 * (fs[:-1] + fs[1:]) * np.diff(ts)) / (ts[-1] - ts[0]))  # trapezoid
        else:
            means[L] = float(fs.mean())
    return means, tag

def gate_pair_drift(d_dry, d_ice):
    """PAIRED matched-season convergence (CS 2026-09-03): the A1 decider. G5 gates each run separately, but
    two runs can both pass while their DIFFERENCE still drifts (opposite-sign drifts). Compute the ice-dry
    per-lunation floor-mean difference and test its matched-season (L vs L+12) drift; fail if >0.1 K/yr."""
    md_dry, tag = _per_lunation_means(d_dry)
    md_ice, _ = _per_lunation_means(d_ice)
    if md_dry is None or md_ice is None:
        return None, "need two convergence npz with floor history"
    diff = {L: md_ice[L] - md_dry[L] for L in md_dry if L in md_ice}
    drifts = [diff[L + 12] - diff[L] for L in sorted(diff) if (L + 12) in diff]
    if not drifts:
        return None, f"run too short for matched-season differential drift (<2 yr; {len(diff)} shared lunations)"
    md = float(np.mean(drifts))
    last = diff[max(diff)]
    return abs(md) < 0.1, (f"differential (ice-dry) matched-season drift {md:+.3f} K/yr over {len(drifts)} pairs "
                           f"(need |drift|<0.1); last-lunation ice-dry {last:+.3f} K; {tag} floor")

def gate_drift(d):
    """G5 matched-season convergence (CS 2026-09-02): a seasonal run whose floor still drifts
    year-over-year has not equilibrated, so no ice-signal magnitude/sign can be claimed from it. On the
    per-facet floor history, take the time-weighted floor-mean per lunation and compare lunation L to L+12
    (~1 yr). Fail if |mean matched-season drift| > 0.1 K/yr. Needs a convergence npz (T_surf_crater_history);
    skips otherwise. Prefers the polygon_floor mask over elev-p20."""
    means, tag = _per_lunation_means(d)
    if means is None:
        return None, "no crater history/floor mask (not a convergence file)"
    drifts = [means[L + 12] - means[L] for L in sorted(means) if (L + 12) in means]
    if not drifts:
        return None, f"run too short for matched-season drift (<2 yr; {len(means)} lunations)"
    md = float(np.mean(drifts))
    return abs(md) < 0.1, (f"matched-season (L vs L+12) drift {md:+.3f} K/yr over {len(drifts)} pairs "
                           f"(need |drift|<0.1); {tag} floor")

def main():
    # PAIRED differential-drift check (CS 2026-09-03, the A1 decider):
    #   python3 tools/check_science_gates.py --pair <dry_convergence.npz> <ice_convergence.npz>
    if len(sys.argv) >= 4 and sys.argv[1] == "--pair":
        d_dry = np.load(sys.argv[2], allow_pickle=True); d_ice = np.load(sys.argv[3], allow_pickle=True)
        ok, msg = gate_pair_drift(d_dry, d_ice)
        tag = "SKIP" if ok is None else ("PASS" if ok else "FAIL")
        print(f"[{tag}] PAIR drift  {msg}")
        print(f"\nOVERALL: {'PASS — differential converged' if ok else ('SKIP' if ok is None else 'FAIL — differential not converged')}")
        return 0 if ok else 1
    files = sys.argv[1:] or sorted(
        glob.glob(os.path.join(SYNC, "data", "**", "prod_*psr_floor*.npz"), recursive=True) +
        glob.glob(os.path.join(SYNC, "data", "**", "seasonal_*.npz"), recursive=True) +
        glob.glob(os.path.join(SYNC, "data", "**", "*convergence*.npz"), recursive=True) +
        glob.glob(os.path.join(SYNC, "data", "**", "*thermal*.npz"), recursive=True))
    if not files:
        print("no run files found — nothing to check"); return 0
    # optional backfilled depth/k grids for older outputs that lack them
    bf = None
    cand = glob.glob(os.path.join(SYNC, "data", "backfill_depth_k_*.npz"))
    if cand:
        z = np.load(cand[0])
        keys = [k for k in z.files if k.endswith("_depth_m")]
        if keys:
            base = keys[0][:-8]
            bf = (np.asarray(z[base + "_depth_m"]), np.asarray(z[base + "_k_profile"]))

    any_fail = False
    for f in files:
        d = np.load(f, allow_pickle=True)
        print(f"\n=== {os.path.basename(f)} ===")
        for name, res in [("G1 shadow  ", gate_shadow(d)),
                          ("G2 IC drain", gate_ic(d, bf)),
                          ("G3 monotone", gate_monotone(d)),
                          ("G4 forcing ", gate_forcing(d, f)),
                          ("G5 drift   ", gate_drift(d))]:
            ok, msg = res if res else (None, "n/a")
            tag = "SKIP" if ok is None else ("PASS" if ok else "FAIL")
            if ok is False: any_fail = True
            print(f"  [{tag}] {name}  {msg}")
    print(f"\nOVERALL: {'FAIL — see gates above' if any_fail else 'all applicable gates PASS'}")
    return 1 if any_fail else 0

if __name__ == "__main__":
    sys.exit(main())
