"""Audit of the crater iterative radiosity (multiple-scattering) solver for physical correctness.
Checks: (1) energy conservation on a CLOSED enclosure, (2) energy balance on an OPEN mesh (absorbed +
escaped-to-space = input), (3) convergence within max_iter across albedo, (4) all facets solved / no
NaN on a real crater mesh, (5) thermal self-heating = single-bounce VM@(emσT^4), (6) reciprocity of VM."""
import os, sys, numpy as np
ROOT="/Users/phillipsm/Documents/Software/ThermoSpec-3D"; sys.path.insert(0,ROOT)
from crater import (compute_multiple_scattered_sunlight, CraterMesh, SelfHeatingList,
                    CraterRadiativeTransfer)

def scattered_to_total(Fscat, direct):
    return Fscat + direct   # T = scattered + direct

print("="*70)
print("(1) ENERGY CONSERVATION — closed enclosure (row sums = 1, reciprocity)")
# synthetic closed enclosure: equal areas, symmetric F, rows sum to 1, zero diagonal
N=40; rng=np.arange(N)
M=np.abs(np.add.outer(np.sin(rng), np.cos(rng)))+0.1; np.fill_diagonal(M,0)
F=0.5*(M+M.T); F/=F.sum(1,keepdims=True)     # symmetric + row-normalized -> reciprocity w/ equal areas
for A in (0.1,0.3,0.6):
    Alb=np.full((N,1),A); illum=np.zeros((N,1)); illum[:N//2]=1.0
    cos=np.full((N,1),0.7); Fsun=1361.0
    Fscat=compute_multiple_scattered_sunlight(Alb,Fsun,illum,cos,F,max_iter=500,tol=1e-10)
    direct=Fsun*illum*cos; T=scattered_to_total(Fscat,direct)
    absorbed=((1-A)*T).sum(); inp=direct.sum()   # equal areas -> sum is energy up to A_facet const
    print(f"  A={A}: input={inp:.1f}  absorbed={absorbed:.1f}  ratio={absorbed/inp:.6f}  (closed: expect 1.000000)")

print("="*70)
print("(2) ENERGY BALANCE — real crater mesh (open: absorbed + escaped = input)")
MESH=os.path.join(ROOT,"Roughness_files","new_crater2.txt")
VF=os.path.join(ROOT,"Roughness_files","new_crater2_selfheating_list.txt")
mesh=CraterMesh(MESH); sh=SelfHeatingList(VF); N=len(mesh.normals)
rt=CraterRadiativeTransfer(mesh,sh); VM=rt.view_matrix
A_f=mesh.areas; rowsum=VM.sum(1)
print(f"  mesh: {N} facets; VM row-sum mean {rowsum.mean():.3f} (open, <1 = sky view)")
for A in (0.1,0.3):
    Alb=np.full((N,1),A)
    sun=np.array([0.3,0.0,0.95]); sun/=np.linalg.norm(sun)
    cosv=np.clip(mesh.normals@sun,0,None)[:,None]; illum=(cosv>0).astype(float)
    Fsun=1361.0
    Fscat=compute_multiple_scattered_sunlight(Alb,Fsun,illum,cosv,VM,max_iter=500,tol=1e-12)
    direct=Fsun*illum*cosv; T=Fscat+direct; B=A*T           # radiosity leaving each facet
    absorbed=(A_f[:,None]*(1-A)*T).sum()
    escaped=(A_f[:,None]*(1-rowsum[:,None])*B).sum()          # radiosity fraction to sky
    inp=(A_f[:,None]*direct).sum()
    print(f"  A={A}: input={inp:.3e}  absorbed+escaped={absorbed+escaped:.3e}  ratio={(absorbed+escaped)/inp:.6f} (expect 1.0)")

print("="*70)
print("(3) CONVERGENCE within max_iter across albedo")
for A in (0.1,0.5,0.9):
    Alb=np.full((N,1),A)
    sun=np.array([0.3,0,0.95]); sun/=np.linalg.norm(sun)
    cosv=np.clip(mesh.normals@sun,0,None)[:,None]; illum=(cosv>0).astype(float)
    G=Alb*(1361.0*illum*cosv); direct=1361.0*illum*cosv
    for it in range(1,201):
        Gn=Alb*(direct+VM@G)
        if np.allclose(Gn,G,rtol=1e-8,atol=1e-8): break
        G=Gn
    print(f"  A={A}: converged in {it} iters (production max_iter=100){'  <-- WARN: exceeds 100' if it>100 else ''}")

print("="*70)
print("(4) ALL FACETS SOLVED / no NaN (real mesh)")
Q_dir,Q_scat,Q_self,cos=rt.compute_fluxes(np.array([0.3,0,0.95]),
    (np.clip(mesh.normals@(np.array([0.3,0,0.95])/np.linalg.norm([0.3,0,0.95])),0,None)>0).astype(float),
    0.95*5.67e-8*np.full(N,250.0)**4, 0.1, 0.95, 1361.0, 1)
print(f"  Q_dir finite {np.isfinite(Q_dir).all()} ({(Q_dir>0).sum()}/{N} lit); Q_scat finite {np.isfinite(Q_scat).all()}; Q_self finite {np.isfinite(Q_self).all()}")
print(f"  every facet has a Q value: {Q_dir.shape[0]==N and Q_scat.shape[0]==N and Q_self.shape[0]==N}")

print("="*70)
print("(5) THERMAL SELF-HEATING = single-bounce VM@(emσT^4)")
T4=np.full(N,250.0); tf=0.95*5.67e-8*T4**4
Qself_direct=VM@tf
print(f"  max|Q_self - VM@therm_flux| = {np.max(np.abs(Q_self-Qself_direct)):.2e}  (0 = matches single-bounce)")
print("  NOTE: thermal is single-bounce (no thermal multiple-scattering); ~(1-em) error per omitted bounce")

print("="*70)
print("(6) VM RECIPROCITY  A_i F_ij == A_j F_ji")
resid=np.max(np.abs(A_f[:,None]*VM - (A_f[:,None]*VM).T))
print(f"  max|A_i F_ij - A_j F_ji| = {resid:.2e} (0 = exact reciprocity)")
