#!/usr/bin/env python3
import glob, os
import h5py
import numpy as np
import matplotlib.pyplot as plt

def load_snapshot(path, field):
    with h5py.File(path, 'r') as f:
        return f[field][:]

def compute_total(field, dx, dy, dz):
    return np.sum(field) * dx * dy * dz

def main():
    # grid & box (must match parameters.f90)
    Nkx, Nky, Nkz = 64, 64, 64
    Lx, Ly, Lz = 1.0, 1.0, 1.0
    dx, dy, dz = Lx/Nkx, Ly/Nky, Lz/Nkz

    files = sorted(glob.glob('../state_*.h5'))
    if not files:
        print("No snapshots found.")
        return

    # 1) Mass conservation plot
    masses = []
    steps  = []
    for fn in files:
        step = int(os.path.basename(fn)[6:10])
        rho = load_snapshot(fn, 'rho')
        masses.append(compute_total(rho, dx, dy, dz))
        steps.append(step)

    plt.figure()
    plt.plot(steps, masses, 'o-')
    plt.xlabel('Step')
    plt.ylabel('Total mass')
    plt.title('Mass conservation check')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    plt.savefig('mass_conservation.png')
    print("Saved mass_conservation.png")

    # 2) Density slices for first & last
    for step in (steps[0], steps[-1]):
        fn = f'../state_{step:04d}.h5'
        rho = load_snapshot(fn, 'rho')
        sl  = Nkz//2

        plt.figure()
        plt.imshow(rho[:,:,sl].real, origin='lower', aspect='equal')
        plt.colorbar(label='ρ')
        plt.title(f'Density slice at step {step}')
        plt.tight_layout()
        plt.show()
        outname = f'density_slice_{step:04d}.png'
        plt.savefig(outname)
        print(f"Saved {outname}")

if __name__=='__main__':
    main()
