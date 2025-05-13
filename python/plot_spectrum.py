#!/usr/bin/env python3
import numpy as np
import h5py
import matplotlib.pyplot as plt

# Parameters: change as needed
FILENAME = '../state_0100.h5'  # pick a snapshot
Nkx, Nky, Nkz = 64, 64, 64     # your grid sizes

# 1) Load k-space fields from HDF5
with h5py.File(FILENAME, 'r') as f:
    vx_k = f['vx'][:]  # complex64 array shape (Nkx,Nky,Nkz)
    vy_k = f['vy'][:]
    vz_k = f['vz'][:]
    bx_k = f['bx'][:]
    by_k = f['by'][:]
    bz_k = f['bz'][:]

# 2) Compute energy density per mode: 0.5*(|v|^2 + |B|^2)
E_k = 0.5 * ( np.abs(vx_k)**2 + np.abs(vy_k)**2 + np.abs(vz_k)**2
            + np.abs(bx_k)**2 + np.abs(by_k)**2 + np.abs(bz_k)**2 )

# 3) Build k-grid and compute |k| for each point
kx = np.fft.fftfreq(Nkx) * Nkx
ky = np.fft.fftfreq(Nky) * Nky
kz = np.fft.fftfreq(Nkz) * Nkz
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
Kmag = np.sqrt(KX**2 + KY**2 + KZ**2).ravel()

# 4) Bin into spherical shells
E_flat = E_k.ravel()
k_bins = np.arange(0.5, np.max(Kmag)+1, 1.0)
bin_inds = np.digitize(Kmag, k_bins)
spectrum = np.zeros(len(k_bins))
counts   = np.zeros(len(k_bins))
for i, eb in zip(bin_inds, E_flat):
    if 0 < i < len(k_bins):
        spectrum[i] += eb
        counts[i]   += 1
# Avoid division by zero
nonzero = counts > 0
spectrum[nonzero] /= counts[nonzero]

# 5) Plot E(k)
ks = k_bins
plt.loglog(ks[nonzero], spectrum[nonzero], '-o')
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel(r'$E(|\mathbf{k}|)$')
plt.title('1D Energy Spectrum')
plt.grid(True, which='both')
plt.tight_layout()
plt.savefig('energy_spectrum.png')
plt.show()
