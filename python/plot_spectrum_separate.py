#!/usr/bin/env python3
import numpy as np
import h5py
import matplotlib.pyplot as plt

# ─── Configuration ─────────────────────────────────────────────────────────
# Path to the HDF5 snapshot you want to analyze
FILENAME = '../state_0100.h5'

# Your grid sizes (must match what you ran)
Nkx, Nky, Nkz = 64, 64, 64

# ─── Load k‑space fields ────────────────────────────────────────────────────
with h5py.File(FILENAME, 'r') as f:
    vx_k = f['vx'][:]
    vy_k = f['vy'][:]
    vz_k = f['vz'][:]
    bx_k = f['bx'][:]
    by_k = f['by'][:]
    bz_k = f['bz'][:]

# ─── Compute per‑mode energy densities ───────────────────────────────────────
Ek_density = 0.5 * ( np.abs(vx_k)**2 + np.abs(vy_k)**2 + np.abs(vz_k)**2 )
Eb_density = 0.5 * ( np.abs(bx_k)**2 + np.abs(by_k)**2 + np.abs(bz_k)**2 )

# ─── Build |k| grid ─────────────────────────────────────────────────────────
kx = np.fft.fftfreq(Nkx) * Nkx
ky = np.fft.fftfreq(Nky) * Nky
kz = np.fft.fftfreq(Nkz) * Nkz
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
Kmag = np.sqrt(KX**2 + KY**2 + KZ**2).ravel()

# ─── Bin into spherical shells ──────────────────────────────────────────────
Ekin_spectrum = []
Emag_spectrum = []
k_bins = np.arange(0.5, np.max(Kmag)+1, 1.0)
counts = np.zeros(len(k_bins))
Ekin_spectrum = np.zeros(len(k_bins))
Emag_spectrum = np.zeros(len(k_bins))

inds = np.digitize(Kmag, k_bins)
for idx, ek, em in zip(inds, Ek_density.ravel(), Eb_density.ravel()):
    if 0 < idx < len(k_bins):
        Ekin_spectrum[idx] += ek
        Emag_spectrum[idx] += em
        counts[idx]       += 1

# Avoid division by zero
mask = counts > 0
Ekin_spectrum[mask] /= counts[mask]
Emag_spectrum[mask] /= counts[mask]

# ─── Plotting ───────────────────────────────────────────────────────────────
ks = k_bins

plt.loglog(ks[mask], Ekin_spectrum[mask], '-o', label='Kinetic')
plt.loglog(ks[mask], Emag_spectrum[mask], '-o', label='Magnetic')

# Optional: overlay a k^-5/3 slope for reference
C = Ekin_spectrum[mask][10] * ks[mask][10]**(5/3)
plt.loglog(ks[mask], C * ks[mask]**(-5/3), '--', label=r'$k^{-5/3}$')

plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel(r'$E(k)$')
plt.title('1D Kinetic & Magnetic Energy Spectra')
plt.legend()
plt.grid(True, which='both')
plt.tight_layout()
plt.savefig('energy_spectrum_separate.png')
plt.show()
