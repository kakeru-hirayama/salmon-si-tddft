"""
フェルミ-ディラック分布フィッティング → 電子温度 Te 算出
堀井修士論文 4.4節の手順に基づく

【何をやってるか】
  レーザー照射後の電子エネルギー分布 f(ε) を、
  フェルミ-ディラック分布 f_FD(ε; T, μ) にフィットして
  電子温度 Te を逆算する。

【必要なファイル】
  - Si_eigen_after.data : パルス照射後の (エネルギー, 占有数) per k点
  - dos.data            : GS計算から得た状態密度 DOS(ε)
  - Si_rt_energy.data   : 励起エネルギーの時間発展（Te推定のクロスチェック用）

【SALMONの出力ファイルについて】
  tddft_pulse計算後に以下を確認:
    Si_eigen.data (照射後の占有数が含まれる場合)
    または analysis セクションに yn_out_occ = 'y' を追加して再実行
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import Boltzmann, elementary_charge

# ---- 定数 ----
kB_eV = Boltzmann / elementary_charge   # ボルツマン定数 [eV/K]

# ==============================================================================
# 1. フェルミ-ディラック分布関数
# ==============================================================================
def fermi_dirac(eps, mu, T):
    """
    f(ε) = 1 / (1 + exp((ε - μ) / kB·T))
    eps : エネルギー [eV]
    mu  : 化学ポテンシャル [eV]
    T   : 温度 [K]
    """
    x = (eps - mu) / (kB_eV * T)
    # オーバーフロー防止
    x = np.clip(x, -500, 500)
    return 1.0 / (1.0 + np.exp(x))


# ==============================================================================
# 2. SALMONのeigen dataを読み込む（k点平均）
# ==============================================================================
def load_eigen_averaged(filepath, nstate=64):
    """
    Si_eigen.data を読み込んで k点平均した
    (エネルギー[eV], 平均占有数) を返す。

    ファイル形式:
      k= N, spin= 1
        io  esp[eV]  occ
        ...
    """
    energies_all = []
    occs_all = []

    with open(filepath, 'r') as f:
        lines = f.readlines()

    current_eps  = []
    current_occ  = []

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('k='):
            # 前のk点データを保存
            if current_eps:
                energies_all.append(current_eps)
                occs_all.append(current_occ)
            current_eps = []
            current_occ = []
        else:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    eps = float(parts[1])
                    occ = float(parts[2])
                    current_eps.append(eps)
                    current_occ.append(occ)
                except ValueError:
                    pass

    # 最後のk点
    if current_eps:
        energies_all.append(current_eps)
        occs_all.append(current_occ)

    # k点平均（全k点で同じ軌道インデックスのエネルギーと占有数を平均）
    energies_all = np.array(energies_all)  # shape: (nk, nstate)
    occs_all     = np.array(occs_all)      # shape: (nk, nstate)

    eps_avg = energies_all.mean(axis=0)    # k平均エネルギー
    occ_avg = occs_all.mean(axis=0)        # k平均占有数（0〜2）

    return eps_avg, occ_avg


# ==============================================================================
# 3. DOSを読み込む
# ==============================================================================
def load_dos(filepath):
    """
    dos.data を読み込む。
    形式: Energy[eV]  DOS[1/eV]
    """
    data = np.loadtxt(filepath, comments='#')
    eps_dos = data[:, 0]
    dos     = data[:, 1]
    return eps_dos, dos


# ==============================================================================
# 4. フェルミフィット実行
# ==============================================================================
def fit_fermi(eps, occ_per_state, label='', plot=True):
    """
    各軌道の占有数（0〜2） → 規格化して f(ε)（0〜1）にしてフィット。
    占有数の最大値が2なので 2 で割る。

    eps          : エネルギー配列 [eV]
    occ_per_state: 占有数配列（0〜2）
    """
    f_data = occ_per_state / 2.0   # 0〜1 に規格化

    # 価電子帯（occ > 0.1）と伝導帯（occ < 0.9）の境界付近を使ってフィット
    # バンドギャップ付近は飛ばす（2.5 eV 以下 or 3.5 eV 以上）
    mask_vb = (eps < 4.5) & (f_data > 0.05)    # 価電子帯側
    mask_cb = (eps > 2.5) & (f_data < 0.95)    # 伝導帯側
    mask = mask_vb | mask_cb

    eps_fit  = eps[mask]
    f_fit    = f_data[mask]

    if len(eps_fit) < 5:
        print(f"[{label}] フィット用データ点が少なすぎます（{len(eps_fit)}点）")
        return None, None

    # 初期値：μ ≈ フェルミエネルギー付近、T ≈ 5000 K
    eps_fermi_guess = eps[np.argmin(np.abs(f_data - 0.5))]
    p0 = [eps_fermi_guess, 5000.0]

    try:
        popt, pcov = curve_fit(fermi_dirac, eps_fit, f_fit,
                               p0=p0, bounds=([-np.inf, 1], [np.inf, 1e7]),
                               maxfev=10000)
        mu_fit, Te_fit = popt
        perr = np.sqrt(np.diag(pcov))
        print(f"[{label}]")
        print(f"  化学ポテンシャル μ  = {mu_fit:.4f} eV")
        print(f"  電子温度 Te         = {Te_fit:.0f} K  ({Te_fit * kB_eV * 1000:.1f} meV)")
        print(f"  フィット誤差（T）   = ±{perr[1]:.0f} K")
        print(f"  参考（堀井論文）    10¹² W/cm² → ~1,330 K")
    except RuntimeError as e:
        print(f"[{label}] フィット失敗: {e}")
        return None, None

    # ---- プロット ----
    if plot:
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.scatter(eps, f_data, s=15, alpha=0.6, color='tomato', label='TDDFT occupation')
        eps_plot = np.linspace(eps.min(), eps.max(), 500)
        ax.plot(eps_plot, fermi_dirac(eps_plot, mu_fit, Te_fit),
                color='royalblue', lw=2,
                label=f'Fermi-Dirac fit\n$T_e$ = {Te_fit:.0f} K\n$\\mu$ = {mu_fit:.3f} eV')
        ax.axvline(mu_fit, color='royalblue', lw=0.8, ls='--', alpha=0.5)
        ax.axhline(0.5, color='gray', lw=0.5, ls=':')
        ax.set_xlabel('Energy (eV)', fontsize=12)
        ax.set_ylabel('Occupation f(ε)', fontsize=12)
        ax.set_title(f'Fermi-Dirac Fit  [{label}]\n'
                     f'TBmBJ + 16×16×16,  I = 10¹² W/cm²', fontsize=11)
        ax.legend(fontsize=10)
        ax.grid(alpha=0.3)
        ax.set_ylim(-0.05, 1.15)
        plt.tight_layout()
        outname = f'fermi_fit_{label.replace(" ","_")}.png'
        plt.savefig(outname, dpi=150, bbox_inches='tight')
        print(f"  グラフ保存: {outname}")
        plt.show()

    return mu_fit, Te_fit


# ==============================================================================
# 5. 励起エネルギーからの Te 概算（クロスチェック用）
# ==============================================================================
def estimate_Te_from_excitation(rt_energy_file, dos_file):
    """
    照射後の励起エネルギー ΔE から、Te を粗く見積もる。
    ΔE ≈ ∫ ε・DoS・[f_FD(ε;Te) - f_FD(ε;0)] dε
    （厳密なフィットの前のクロスチェック用）
    """
    data = np.loadtxt(rt_energy_file, comments='#')
    delta_E = data[-1, 2]   # 最終時刻の励起エネルギー [eV/cell]
    print(f"\n[励起エネルギーチェック]")
    print(f"  ΔE (照射後) = {delta_E:.4f} eV/unit cell")
    print(f"  堀井論文 10¹² W/cm²: ~1.17 eV/cell")
    return delta_E


# ==============================================================================
# 6. メイン実行
# ==============================================================================
if __name__ == '__main__':
    import os

    # ---- パス設定 ----
    # ※ 計算完了後に Si_eigen_after.data をサーバーからダウンロードして配置
    eigen_after = 'Si_eigen_after.data'     # パルス照射後の占有数
    eigen_gs    = '../../high_accuracy_gs/Si_eigen.data'   # GS（比較用）
    dos_file    = '../../high_accuracy_gs/dos.data'         # DoS
    rt_energy   = 'Si_rt_energy.data'       # 励起エネルギー（ダウンロード後）

    # ---- GS占有数の確認（基底状態：全部0か2のはず）----
    if os.path.exists('../../../high_accuracy_gs/Si_eigen.data'.replace('../../../', os.path.expanduser('~') + '/SALMON/calc/')):
        pass

    # ---- 励起エネルギーのクロスチェック ----
    if os.path.exists(rt_energy):
        estimate_Te_from_excitation(rt_energy, dos_file)
    else:
        print(f"[INFO] {rt_energy} が見つかりません。計算完了後に配置してください。")

    # ---- フェルミフィット ----
    if os.path.exists(eigen_after):
        eps, occ = load_eigen_averaged(eigen_after)
        fit_fermi(eps, occ, label='after pulse (10^12 W/cm2)')
    else:
        print(f"\n[INFO] {eigen_after} が見つかりません。")
        print("  計算完了後に以下を実行:")
        print("    scp hirayama@192.168.200.20:~/SALMON/calc/high_accuracy_rt_pulse/Si_eigen.data ./Si_eigen_after.data")
        print("    scp hirayama@192.168.200.20:~/SALMON/calc/high_accuracy_rt_pulse/Si_rt_energy.data ./")
        print("  その後このスクリプトを再実行してください。")

        # ---- デモ：GS占有数でフィットしてみる（T ≈ 0 に収束するはず）----
        print("\n[デモ] GS占有数（T≈0）でフィットを試します...")
        gs_path = os.path.join(os.path.dirname(__file__),
                               '../../exercise5_dielectric/gs/Si_eigen.data')
        # GS は occupation が全部0か2なのでフィットできないことを確認
        print("  ※ GS は T≈0 なのでフィット不能（正常）")
