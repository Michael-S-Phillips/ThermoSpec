"""Wall-time benchmark for the LOD-ADI 3D conduction step.

Answers the Phase-1 gate: is an implicit 3D conduction solve fast enough at lunar cadence
(~800k steps for 8 lunations)? Times step() at a few grid sizes and projects the full run.
The vertical resolution (nz~100) matches the real 1D column; lateral (nx=ny) is what we trade.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from adi3d import ADI3DSolver

K = 7.7e-10          # regolith diffusivity m^2/s (k/(rho cp) ~ 7.4e-4/(1100*825))
DT = 29.5306 * 24 * 3600 / 100000.0   # ~25.5 s, the lunar diurnal step
STEPS_LUNAR = 100000 * 8              # 8 lunations at tsteps_day=100000


def bench(nx, ny, nz, n_time=200):
    solver = ADI3DSolver((nx, ny, nz), (0.05, 0.05, 0.005), K, DT)
    T = np.random.default_rng(0).random((nx, ny, nz))
    solver.step(T)  # warm up
    t0 = time.perf_counter()
    for _ in range(n_time):
        T = solver.step(T)
    dt = time.perf_counter() - t0
    rate = n_time / dt
    unknowns = nx * ny * nz
    full_h = STEPS_LUNAR / rate / 3600.0
    print(f"  {nx:3d}x{ny:3d}x{nz:3d} = {unknowns:8d} unknowns | {rate:8.1f} steps/s "
          f"| 800k-step lunar run ~ {full_h:6.2f} h")
    return rate


if __name__ == "__main__":
    print("LOD-ADI 3D conduction wall-time (CPU, single thread), dt=%.1fs" % DT)
    print("1D baseline for reference: Hapke core ~1e4 steps/s on 177 nodes.\n")
    for (nx, ny, nz) in [(16, 16, 100), (32, 32, 100), (48, 48, 100), (64, 64, 100),
                         (100, 100, 100)]:
        bench(nx, ny, nz)
