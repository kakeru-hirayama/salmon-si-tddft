"""
Slack-ready: Si dielectric function (High Accuracy, TBmBJ)
With parameter table and interpretation annotations.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ---- Data ----
resp = np.loadtxt('Si_response_TBmBJ.data', comments='#')
e    = resp[:, 0]
re   = resp[:, 9]
im   = resp[:, 12]

m = (e >= 1.0) & (e <= 8.0)
e_m, re_m, im_m = e[m], re[m], im[m]

# ---- Figure ----
fig = plt.figure(figsize=(12, 12))
fig.patch.set_facecolor('white')

# GridSpec: table + 2 plots, with room at top for title
gs = gridspec.GridSpec(3, 1, height_ratios=[1.0, 2.2, 2.2], hspace=0.42,
                       top=0.80, bottom=0.06, left=0.10, right=0.95)

fig.text(0.5, 0.97,
         'Si Dielectric Function (High Accuracy)',
         ha='center', va='center', fontsize=14, fontweight='bold')
fig.text(0.5, 0.93,
         'TBmBJ + 16x16x16 k-grid  |  Ishikawa Lab, K. Hirayama  |  SALMON RT-TDDFT',
         ha='center', va='center', fontsize=11, color='#444444')

# ---- Top: parameter table ----
ax_table = fig.add_subplot(gs[0])
ax_table.axis('off')
table_data = [
    ['Functional',         'TBmBJ (Tran-Blaha modified Becke-Johnson)'],
    ['k-grid',             '16x16x16 = 4,096 pts  (Tutorial: 4x4x4 = 64 pts)'],
    ['r-grid',             '16x16x16  (Tutorial: 12x12x12)'],
    ['nstate',             '64  (16 occupied + 48 empty)'],
    ['SCF convergence',    '68 iterations  (density change < 1e-9)'],
    ['Indirect band gap',  '1.033 eV  (Exp: 1.17 eV)'],
    ['Completed',          'GS: 2026-04-09 22:27  /  Linear response: 2026-04-11 21:07'],
]
tbl = ax_table.table(
    cellText=table_data,
    colLabels=['Parameter', 'Value'],
    cellLoc='left', loc='center',
    colWidths=[0.25, 0.72],
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1, 1.5)
for (r, c), cell in tbl.get_celld().items():
    if r == 0:
        cell.set_facecolor('#2a6fa8')
        cell.set_text_props(color='white', fontweight='bold')
    else:
        cell.set_facecolor('#f7f9fc' if r % 2 == 0 else 'white')
    cell.set_edgecolor('#cccccc')

# ---- Im(epsilon) ----
ax1 = fig.add_subplot(gs[1])
ax1.plot(e_m, im_m, color='#e74c3c', lw=2.2, label='Im(eps_z)  TBmBJ', zorder=3)
ax1.axhline(0, color='gray', lw=0.6, ls='--', zorder=1)
ax1.axvline(2.80, color='#e67e22', lw=1.3, ls=':', zorder=2)
ax1.axvline(3.40, color='#27ae60', lw=1.3, ls='-.', zorder=2)
ax1.axvline(4.27, color='#8e44ad', lw=1.3, ls='-.', zorder=2)

ax1.annotate('TBmBJ onset\n~2.80 eV\n(direct gap start)',
             xy=(2.80, 4), xytext=(1.9, 22),
             fontsize=8.5, color='#e67e22',
             arrowprops=dict(arrowstyle='->', color='#e67e22', lw=1.2), ha='center')
ax1.annotate('Exp. E1 ~3.4 eV\n(direct gap)',
             xy=(3.40, 26), xytext=(3.95, 8),
             fontsize=8.5, color='#27ae60',
             arrowprops=dict(arrowstyle='->', color='#27ae60', lw=1.2), ha='center')
ax1.annotate('Calc. peak\n~3.97 eV',
             xy=(3.97, 47.9), xytext=(5.3, 42),
             fontsize=8.5, color='#e74c3c',
             arrowprops=dict(arrowstyle='->', color='#e74c3c', lw=1.2), ha='center')
ax1.annotate('Exp. E2 ~4.27 eV',
             xy=(4.27, 44), xytext=(4.27, 30),
             fontsize=8.5, color='#8e44ad',
             arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=1.2), ha='center')

ax1.text(6.3, 48,
         'LDA onset: 2.4 eV\n-> TBmBJ: 2.80 eV\n-> Exp: ~3.4 eV\n(improved!)',
         fontsize=8.5, va='top', ha='center',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#e8f5e9',
                   edgecolor='#a5d6a7', alpha=0.9))

ax1.set_ylabel('Im(eps)', fontsize=12)
ax1.set_title('Im(eps): Absorption spectrum  [light absorbed at this energy]', fontsize=11, pad=6)
ax1.set_xlim(1, 8); ax1.set_ylim(-3, 55)
ax1.grid(alpha=0.3); ax1.tick_params(labelbottom=False)

# ---- Re(epsilon) ----
ax2 = fig.add_subplot(gs[2])
ax2.plot(e_m, re_m, color='#2980b9', lw=2.2, zorder=3)
ax2.axhline(0, color='gray', lw=0.6, ls='--', zorder=1)
ax2.axhline(11.7, color='#7f8c8d', lw=1.0, ls=':', zorder=2)
ax2.axvline(3.40, color='#27ae60', lw=1.3, ls='-.', zorder=2)

ax2.annotate('Static dielectric\nRe(eps->0) ~12.2\n(Exp ~11.7)',
             xy=(1.5, 12.2), xytext=(2.6, 28),
             fontsize=8.5, color='#7f8c8d',
             arrowprops=dict(arrowstyle='->', color='#7f8c8d', lw=1.2), ha='center')
ax2.annotate('Sign change ~4.0 eV\n(near Im peak)',
             xy=(4.02, 0), xytext=(5.0, 12),
             fontsize=8.5, color='#2980b9',
             arrowprops=dict(arrowstyle='->', color='#2980b9', lw=1.2), ha='center')

ax2.set_ylabel('Re(eps)', fontsize=12)
ax2.set_xlabel('Energy (eV)', fontsize=12)
ax2.set_title('Re(eps): Related to refractive index and reflectance', fontsize=11, pad=6)
ax2.set_xlim(1, 8); ax2.set_ylim(-22, 40)
ax2.grid(alpha=0.3)


plt.savefig('Si_dielectric_TBmBJ_slack.png', dpi=150, bbox_inches='tight',
            facecolor='white')
print("Saved: Si_dielectric_TBmBJ_slack.png")
plt.show()
