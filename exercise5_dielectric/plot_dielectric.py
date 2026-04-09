"""
Exercise-5: Dielectric function of crystalline Si
Plot: (1) Matter current density J_z(t)
      (2) Re(epsilon_z)
      (3) Im(epsilon_z)

Data files (from SALMON output in rt/ directory):
  Si_rt.data       -- real-time current density
  Si_rt_energy.data -- total energy vs time
  Si_response.data  -- Fourier-transformed dielectric function
"""

import numpy as np
import matplotlib.pyplot as plt

# ---- Load data ----
rt       = np.loadtxt('Si_rt.data',       comments='#')
response = np.loadtxt('Si_response.data', comments='#')

# Si_rt.data: col1=Time[fs], col16=Jm_z[1/fs/Ang^2]
t    = rt[:, 0]
Jm_z = rt[:, 15]

# Si_response.data: col1=Energy[eV], col10=Re(eps_z), col13=Im(eps_z)
e_eV   = response[:, 0]
re_eps = response[:, 9]
im_eps = response[:, 12]

# Mask: 1-8 eV
mask = (e_eV >= 1.0) & (e_eV <= 8.0)
e  = e_eV[mask]
re = re_eps[mask]
im = im_eps[mask]

# ---- Plot ----
fig, axes = plt.subplots(3, 1, figsize=(9, 11))
fig.suptitle('Si RT-TDDFT Exercise-5: Dielectric Function\n'
             '(LDA/PZ, 4×4×4 k-grid, 12×12×12 r-grid)', fontsize=13)

# (1) Matter current density
ax = axes[0]
ax.plot(t, Jm_z, color='steelblue', lw=1.0)
ax.set_xlabel('Time (fs)', fontsize=12)
ax.set_ylabel('$J_z$ (1/fs/Å²)', fontsize=12)
ax.set_title('Matter Current Density $J_z(t)$ after impulse', fontsize=12)
ax.axhline(0, color='gray', lw=0.5, ls='--')
ax.grid(alpha=0.3)

# (2) Re(epsilon)
ax = axes[1]
ax.plot(e, re, color='royalblue', lw=1.5, label='Re(ε_z)')
ax.axhline(0, color='gray', lw=0.5, ls='--')
ax.axvline(2.4, color='green', lw=1.0, ls=':', label='LDA gap ≈ 2.4 eV')
ax.set_ylabel('Re(ε)', fontsize=12)
ax.set_title('Dielectric Function — Real Part', fontsize=12)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (3) Im(epsilon)
ax = axes[2]
ax.plot(e, im, color='tomato', lw=1.5, label='Im(ε_z)')
ax.axhline(0, color='gray', lw=0.5, ls='--')
ax.axvline(2.4, color='green', lw=1.0, ls=':', label='LDA gap ≈ 2.4 eV')
ax.set_ylabel('Im(ε)', fontsize=12)
ax.set_xlabel('Energy (eV)', fontsize=12)
ax.set_title('Dielectric Function — Imaginary Part', fontsize=12)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('results/Si_ex05_result.png', dpi=150, bbox_inches='tight')
print("Saved: results/Si_ex05_result.png")
plt.show()
