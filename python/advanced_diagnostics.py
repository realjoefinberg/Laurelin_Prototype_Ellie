#!/usr/bin/env python3
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.fft import ifftn, fftn, fftfreq

# ─── CONFIG ──────────────────────────────────────────────────────────────────
# List of snapshot files to analyze
snapshots = ['../state_0001.h5', '../state_0050.h5', '../state_0100.h5']
Nkx, Nky, Nkz = 64, 64, 64
Lx, Ly, Lz = 2*np.pi, 2*np.pi, 2*np.pi   # match your box lengths

# Orders for structure functions
orders = [1,2,3,4,5,6]

# Separation lags (in grid units)
lags = np.array([1,2,4,8,16,32])

# ─── HELPERS ─────────────────────────────────────────────────────────────────
def load_fields(fn):
    """Load vx,vy,vz,bx,by,bz in k-space and return real-space arrays."""
    with h5py.File(fn,'r') as f:
        fields_k = {key: f[key][:] for key in ['vx','vy','vz','bx','by','bz']}
    # inverse FFT to get real-space complex fields
    fields_r = {k: np.real(ifftn(vk)) for k,vk in fields_k.items()}
    return fields_r

def compute_cross_helicity(fields_r):
    """H_c = ∫ v·B dV"""
    v = np.stack([fields_r['vx'], fields_r['vy'], fields_r['vz']], axis=0)
    B = np.stack([fields_r['bx'], fields_r['by'], fields_r['bz']], axis=0)
    dot = np.sum(v*B, axis=0)
    return np.mean(dot)   # volume‐averaged

def compute_increments(arr, lag, axis):
    """Compute increments arr(x+lag) - arr(x) along the given axis."""
    return np.take(arr, range(lag,arr.shape[axis]), axis=axis) \
         - np.take(arr, range(0,arr.shape[axis]-lag), axis=axis)

# ─── MAIN ────────────────────────────────────────────────────────────────────
# 1) Cross‑helicity time series
# … after computing BdotGradV and VdotGradB for a chosen snapshot fr …

# Compute the local coupling fields
C1 = bx*BdotGradV[0] + by*BdotGradV[1] + bz*BdotGradV[2]
C2 = vx*VdotGradB[0] + vy*VdotGradB[1] + vz*VdotGradB[2]

# Flatten and take real part
c1_flat = np.real(C1).ravel()
c2_flat = np.real(C2).ravel()

# Plot PDFs
plt.figure()
plt.hist(c1_flat, bins=200, density=True, alpha=0.6, label=r'$B\cdot(B\!\cdot\!\nabla v)$')
plt.hist(c2_flat, bins=200, density=True, alpha=0.6, label=r'$v\cdot(v\!\cdot\!\nabla B)$')
plt.xlabel('Local coupling value')
plt.ylabel('PDF')
plt.legend()
plt.title('PDF of local coupling terms')
plt.tight_layout()
plt.savefig('coupling_pdf.png')
plt.show()


# 2) PDF of increments (at a single lag) for vx
lag = lags[2]  # e.g. lag=4
fr = load_fields(snapshots[-1])
dv = compute_increments(fr['vx'], lag, axis=0).ravel()
plt.figure()
plt.hist(dv, bins=100, density=True)
plt.xlabel(f'Δv_x at lag={lag}')
plt.ylabel('PDF')
plt.title('PDF of velocity increments')
plt.savefig('pdf_vx_increments.png')

# 3) Structure functions S_p(lag) for vx
S = {p: [] for p in orders}
for lag in lags:
    inc = compute_increments(fr['vx'], lag, axis=0).ravel()
    for p in orders:
        S[p].append(np.mean(np.abs(inc)**p))
plt.figure()
for p in orders:
    plt.loglog(lags, S[p], '-o', label=f'p={p}')
plt.xlabel('Lag (grid units)')
plt.ylabel('S_p(lag)')
plt.legend()
plt.title('Structure functions of v_x')
plt.savefig('structure_functions.png')

# 4) Coupling rates (local v·(B·∇B) and B·(v·∇v))
# Compute spectral proxies: ⟨B·(B·∇)v⟩ and ⟨v·(v·∇)B⟩
# Approximate with finite differences
# 4) Time‑series of coupling rates across snapshots
dx, dy, dz = Lx/Nkx, Ly/Nky, Lz/Nkz

# helpers for derivatives
def d_dx(f):
    return (np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0)) / (2*dx)
def d_dy(f):
    return (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2*dy)
def d_dz(f):
    return (np.roll(f, -1, axis=2) - np.roll(f, 1, axis=2)) / (2*dz)

coupling1_ts = []
coupling2_ts = []

for fn in snapshots:
    fr = load_fields(fn)
    vx, vy, vz = fr['vx'], fr['vy'], fr['vz']
    bx, by, bz = fr['bx'], fr['by'], fr['bz']

    # B·(B·∇)v
    BdotGradV = (
      bx*d_dx(vx) + by*d_dy(vx) + bz*d_dz(vx),
      bx*d_dx(vy) + by*d_dy(vy) + bz*d_dz(vy),
      bx*d_dx(vz) + by*d_dy(vz) + bz*d_dz(vz),
    )
    # v·(v·∇)B
    VdotGradB = (
      vx*d_dx(bx) + vy*d_dy(bx) + vz*d_dz(bx),
      vx*d_dx(by) + vy*d_dy(by) + vz*d_dz(by),
      vx*d_dx(bz) + vy*d_dy(bz) + vz*d_dz(bz),
    )

    coupling1_ts.append(np.mean(bx*BdotGradV[0] +
                                by*BdotGradV[1] +
                                bz*BdotGradV[2]))
    coupling2_ts.append(np.mean(vx*VdotGradB[0] +
                                vy*VdotGradB[1] +
                                vz*VdotGradB[2]))

# Plot the time‐series
plt.figure()
plt.plot(coupling1_ts, '-o', label=r'⟨B·(B·∇)v⟩')
plt.plot(coupling2_ts, '-o', label=r'⟨v·(v·∇)B⟩')
plt.xlabel('Snapshot index')
plt.ylabel('Coupling rate')
plt.legend()
plt.title('Energy coupling vs time')
plt.tight_layout()
plt.savefig('coupling_time_series.png')
plt.show()

print(f"⟨B·(B·∇)v⟩ = {coupling1:.6e}")
print(f"⟨v·(v·∇)B⟩ = {coupling2:.6e}")
