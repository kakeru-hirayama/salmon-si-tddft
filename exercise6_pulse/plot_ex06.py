import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ---- データ読み込み ----
rt   = np.loadtxt("Si_rt.data",        comments='#')
ps   = np.loadtxt("Si_pulse.data",     comments='#')
eng  = np.loadtxt("Si_rt_energy.data", comments='#')

# rt.data: col1=time[fs], col16=Jm_z[1/fs/Å²]
time  = rt[:, 0]
Jm_z  = rt[:, 15]

# pulse.data: col1=energy[eV], col10=|Jm_z|²[1/Å⁴]
energy   = ps[:, 0]
power_z  = ps[:, 9]

# rt_energy.data: col1=time[fs], col3=Eall-Eall0[eV]  (= 1 conventional unit cell)
e_time = eng[:, 0]
excit  = eng[:, 2]

# ---- パルスパラメータ ----
omega_eV = 1.55   # 光子エネルギー [eV]

# ---- プロット ----
fig = plt.figure(figsize=(14, 4.5))
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35)

# --- (1) Matter current density ---
ax1 = fig.add_subplot(gs[0])
ax1.plot(time, Jm_z, color='steelblue', lw=0.8)
ax1.set_xlabel("Time [fs]", fontsize=12)
ax1.set_ylabel(r"$J_{m,z}$ [1/fs·Å²]", fontsize=12)
ax1.set_title("Matter Current Density", fontsize=13)
ax1.axhline(0, color='k', lw=0.5, ls='--')
ax1.set_xlim(time[0], time[-1])

# --- (2) Power Spectrum (log scale) ---
ax2 = fig.add_subplot(gs[1])
mask = (energy >= 0.5) & (energy <= 30.0)
ax2.semilogy(energy[mask], power_z[mask], color='crimson', lw=1.0)
# 高次高調波の位置にマーク (n×ω)
for n in range(1, 9):
    x = n * omega_eV
    if x <= 15:
        ax2.axvline(x, color='gray', lw=0.6, ls=':', alpha=0.7)
        ax2.text(x, power_z[mask].max()*0.5, f"{n}ω",
                 ha='center', va='bottom', fontsize=7, color='gray')
ax2.set_xlabel("Energy [eV]", fontsize=12)
ax2.set_ylabel(r"$|J_{m,z}|^2$ [1/Å⁴]", fontsize=12)
ax2.set_title("Power Spectrum", fontsize=13)
ax2.set_xlim(0.5, 15)

# --- (3) Excitation energy per unit cell ---
ax3 = fig.add_subplot(gs[2])
ax3.plot(e_time, excit, color='darkorange', lw=1.0)
ax3.set_xlabel("Time [fs]", fontsize=12)
ax3.set_ylabel("Excitation energy [eV/cell]", fontsize=12)
ax3.set_title("Excitation Energy per Unit Cell", fontsize=13)
ax3.axhline(0, color='k', lw=0.5, ls='--')
ax3.set_xlim(e_time[0], e_time[-1])

fig.suptitle("SALMON Exercise-6: Si under laser pulse  "
             r"($I=10^{12}$ W/cm², $\hbar\omega=1.55$ eV, $T_p=10.672$ fs)",
             fontsize=11, y=1.02)

plt.savefig("ex06_results.png", dpi=150, bbox_inches='tight')
print("Saved: ex06_results.png")
plt.show()
