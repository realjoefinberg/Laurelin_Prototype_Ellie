#!/usr/bin/env python3
import h5py
import numpy as np

# Load the very first snapshot
with h5py.File('../state_0001.h5', 'r') as f:
    vx_k = f['vx'][:]  # complex k-space array

# Take the mid‐plane in z
sl = vx_k.shape[2] // 2
slice_mag = np.abs(vx_k[:, :, sl])

print("vx k-space slice at k_z mid-plane:")
print("  min:", slice_mag.min())
print("  max:", slice_mag.max())
print(" mean:", slice_mag.mean())
