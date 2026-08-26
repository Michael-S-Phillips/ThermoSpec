"""Audit that the radiosity view factors are GEOMETRY-AWARE: if scene geometry blocks the line of
sight between two facets, their view factor must be 0 (no radiative exchange). Tests:
 (A) engineered valley with vs without a tall central blocking wall -> across-valley F drops to 0;
 (B) real bowl: occlusion=True zeros a set of pairs vs occlusion=False, and the zeroed pairs are the
     geometrically-blocked (far/across-bowl) ones, not near neighbours."""
import os, sys, numpy as np
ROOT="/Users/phillipsm/Documents/Software/ThermoSpec-3D"; sys.path.insert(0,ROOT)
from topography import DEMMesh
from view_factors import compute_view_factors

def bowl(n=21, R=9.0, depth=6.0):
    ax=np.arange(n)-(n-1)/2.0; X,Y=np.meshgrid(ax,ax); d=np.sqrt(X**2+Y**2)
    return np.where(d<R, -depth*(1-(d/R)**2), 0.0)

print("="*70)
print("(A) ENGINEERED VALLEY: two inward slopes, with vs without a central blocking wall")
n=21; ax=np.arange(n)-(n-1)/2.0; X,Y=np.meshgrid(ax,ax)
# V-shaped valley along x: z rises with |x| (two slopes facing each other across x=0)
Zvalley = np.abs(X)*1.2
# add a tall thin wall at the valley floor center (x=0 column) that blocks the two slopes
Zwall = Zvalley.copy(); Zwall[:, n//2] += 14.0
for label, Z in [("no wall (open valley)", Zvalley), ("with central wall (blocker)", Zwall)]:
    m=DEMMesh(Z, dx=1.0, dy=1.0, origin="centroid")
    Focc=compute_view_factors(m, occlusion=True,  refine=False)
    Fgeo=compute_view_factors(m, occlusion=False, refine=False)
    cen=m.centroids
    # pick a facet on the far-left slope and one on the far-right slope (should face each other)
    left =int(np.argmax(-cen[:,0])); right=int(np.argmax(cen[:,0]))
    print(f"  {label}:")
    print(f"    geom-only  F(left<->right) = {Fgeo[left,right]:.3e}  (no-occlusion baseline)")
    print(f"    occluded   F(left<->right) = {Focc[left,right]:.3e}")
    print(f"    total pairs zeroed by occlusion: {int((Fgeo>0).sum()-(Focc>0).sum())}")

print("="*70)
print("(B) REAL BOWL: occlusion removes across-bowl (blocked) pairs, keeps near neighbours")
m=DEMMesh(bowl(), dx=1.0, dy=1.0, origin="centroid")
Focc=compute_view_factors(m, occlusion=True,  refine=False)
Fgeo=compute_view_factors(m, occlusion=False, refine=False)
cen=m.centroids
blocked = (Fgeo>0) & (Focc==0)      # pairs the geometry occludes
n_block=int(blocked.sum())
print(f"  bowl {len(m.normals)} facets; pairs zeroed by occlusion: {n_block}")
if n_block:
    ii,jj=np.where(blocked)
    d_block=np.linalg.norm(cen[ii]-cen[jj],axis=1)
    kept=(Focc>0); ki,kj=np.where(kept); d_kept=np.linalg.norm(cen[ki]-cen[kj],axis=1)
    print(f"  mean separation of BLOCKED pairs: {d_block.mean():.2f}  vs KEPT pairs: {d_kept.mean():.2f}")
    print(f"  -> blocked pairs are farther apart (across the bowl, rim between them): "
          f"{'YES' if d_block.mean()>d_kept.mean() else 'NO'}")
print(f"  reciprocity of occluded VF: max|A_iF_ij-A_jF_ji| = "
      f"{np.max(np.abs(m.areas[:,None]*Focc-(m.areas[:,None]*Focc).T)):.2e}")
print(f"  self-view (diagonal) all zero: {np.all(np.diag(Focc)==0)}")
