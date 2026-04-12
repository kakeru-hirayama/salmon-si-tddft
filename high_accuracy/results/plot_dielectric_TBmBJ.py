"""
High Accuracy: Dielectric function of crystalline Si
Functional: TBmBJ, k-grid: 24x24x24, r-grid: 24x24x24, nstate: 64

Compares with Tutorial version (LDA, 4x4x4) and experimental data.

Data file: Si_response_TBmBJ.data
  col1  = Energy [eV]
  col10 = Re(eps_z)
  col13 = Im(eps_z)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ---- Load high accuracy data (TBmBJ) ----
resp_ha = np.loadtxt('Si_response_TBmBJ.data', comments='#')
e_ha    = resp_ha[:, 0]
re_ha   = resp_ha[:, 9]
im_ha   = resp_ha[:, 12]

# ---- Load tutorial data (LDA) for comparison ----
try:
    import os
    lda_path = os.path.join(os.path.dirname(__file__),
                            '../../exercise5_dielectric/results')
    resp_lda = np.loadtxt(os.path.join(lda_path, '../rt/Si_response.data'),
                          comments='#')
    e_lda    = resp_lda[:, 0]
    re_lda   = resp_lda[:, 9]
    im_lda   = resp_lda[:, 12]
    has_lda  = True
except Exception:
    has_lda = False

# ---- Energy mask: 1 - 8 eV ----
def mask(e_arr, emin=1.0, emax=8.0):
    m = (e_arr >= emin) & (e_arr <= emax)
    return m

m_ha = mask(e_ha)

# ---- Experimental reference points (Aspnes & Studna 1983) ----
# Si direct gap E0 = 3.4 eV, E1 ~ 3.4 eV, E2 ~ 4.3 eV
exp_direct_gap = 3.4   # eV
exp_E2         = 4.27  # eV

# ---- Plot ----
fig = plt.figure(figsize=(10, 8))
fig.suptitle('Si Dielectric Function: High Accuracy (TBmBJ, 24×24×24)\nvs Tutorial (LDA, 4×4×4)',
             fontsize=13, y=0.98)

gs_layout = gridspec.GridSpec(2, 1, hspace=0.35)

# (1) Im(epsilon)
ax1 = fig.add_subplot(gs_layout[0])
ax1.plot(e_ha[m_ha], im_ha[m_ha], color='tomato', lw=2.0,
         label='Im(ε) — TBmBJ (High Acc.)')
if has_lda:
    m_l = mask(e_lda)
    ax1.plot(e_lda[m_l], im_lda[m_l], color='tomato', lw=1.2,
             ls='--', alpha=0.5, label='Im(ε) — LDA (Tutorial)')
ax1.axhline(0, color='gray', lw=0.5, ls='--')
ax1.axvline(2.80, color='orangered', lw=1.0, ls=':',
            label='TBmBJ onset ≈ 2.80 eV')
ax1.axvline(exp_direct_gap, color='green', lw=1.0, ls='-.',
            label=f'Exp. direct gap ≈ {exp_direct_gap} eV')
ax1.axvline(exp_E2, color='purple', lw=1.0, ls='-.',
            label=f'Exp. E₂ ≈ {exp_E2} eV')
ax1.set_ylabel('Im(ε)', fontsize=12)
ax1.set_title('Imaginary Part — absorption spectrum', fontsize=11)
ax1.legend(fontsize=9, loc='upper right')
ax1.grid(alpha=0.3)
ax1.set_xlim(1, 8)
ax1.set_ylim(-5, 55)

# (2) Re(epsilon)
ax2 = fig.add_subplot(gs_layout[1])
ax2.plot(e_ha[m_ha], re_ha[m_ha], color='royalblue', lw=2.0,
         label='Re(ε) — TBmBJ (High Acc.)')
if has_lda:
    ax2.plot(e_lda[m_l], re_lda[m_l], color='royalblue', lw=1.2,
             ls='--', alpha=0.5, label='Re(ε) — LDA (Tutorial)')
ax2.axhline(0, color='gray', lw=0.5, ls='--')
ax2.axvline(exp_direct_gap, color='green', lw=1.0, ls='-.',
            label=f'Exp. direct gap ≈ {exp_direct_gap} eV')
ax2.axhline(11.7, color='gray', lw=0.8, ls=':',
            label='Exp. Re(ε)→ ~11.7 (static)')
ax2.set_ylabel('Re(ε)', fontsize=12)
ax2.set_xlabel('Energy (eV)', fontsize=12)
ax2.set_title('Real Part', fontsize=11)
ax2.legend(fontsize=9, loc='upper right')
ax2.grid(alpha=0.3)
ax2.set_xlim(1, 8)
ax2.set_ylim(-25, 45)

plt.savefig('Si_dielectric_TBmBJ.png', dpi=150, bbox_inches='tight')
print("Saved: Si_dielectric_TBmBJ.png")
plt.show()
