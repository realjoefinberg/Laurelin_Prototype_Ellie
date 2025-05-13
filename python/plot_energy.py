#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Load the CSV (skip header)
data = np.loadtxt('../energy.csv', delimiter=',', skiprows=1)
steps, Ekin, Emag = data.T

plt.plot(steps, Ekin, label='Kinetic Energy')
plt.plot(steps, Emag, label='Magnetic Energy')
plt.xlabel('Time Step')
plt.ylabel('Energy')
plt.title('Energy Evolution')
plt.legend()
plt.tight_layout()
plt.savefig('energy_evolution.png')
plt.show()
