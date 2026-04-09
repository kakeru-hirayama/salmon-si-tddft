# High Accuracy Calculations

具志堅修士論文の計算条件を参考に、収束した精度での計算を行うためのinputファイル群。

## 変更点（Tutorial比較）

| パラメータ | Tutorial | High Accuracy | 根拠 |
|---|---|---|---|
| k-grid | 4×4×4 | **24×24×24** | 具志堅論文Fig.3.10: 16×16×16以上で収束 |
| r-grid | 12×12×12 | **24×24×24** | 具志堅論文Fig.3.10: 16×16×16で収束確認 |
| 汎関数 | PZ (LDA) | **TBmBJ** | バンドギャップ精度: LDA~2.6eV → TBmBJ~3.2eV |

## ディレクトリ構成

```
high_accuracy/
├── gs/
│   ├── Si_gs.inp         # DFT基底状態 (TBmBJ, 24×24×24)
│   └── calc_gs.sh        # PBS投入スクリプト (walltime=24h)
├── rt_response/
│   ├── Si_rt_response.inp  # TDDFT線形応答 → 誘電関数 ε(ω)
│   └── calc_rt_response.sh # PBS投入スクリプト (walltime=48h)
└── rt_pulse/
    ├── Si_rt_pulse.inp     # TDDFT パルス照射 → HHG, 励起エネルギー
    └── calc_rt_pulse.sh    # PBS投入スクリプト (walltime=48h)
```

## 実行順序

1. **GS計算**（restart/ が生成される）
   ```bash
   cd gs/
   cp /path/to/Si_rps.dat ./
   qsub calc_gs.sh
   ```

2. **RT応答計算**（GS の restart/ を参照）
   ```bash
   cd ../rt_response/
   cp /path/to/Si_rps.dat ./
   ln -s ../gs/restart ./restart   # または cp -r
   qsub calc_rt_response.sh
   ```

3. **RTパルス計算**（GS の restart/ を参照）
   ```bash
   cd ../rt_pulse/
   cp /path/to/Si_rps.dat ./
   ln -s ../gs/restart ./restart
   qsub calc_rt_pulse.sh
   ```

## 計算コスト見積もり

k-grid が 4→24 倍になると、k点数は (24/4)³ = 216 倍。
ただし並列化（nproc_k=4）で分散するため、実壁時間は数十倍程度になる見込み。
サーバーの空き状況・ノード数に応じて walltime を調整すること。

## 参考文献

- 具志堅英雄 修士論文 (2025): c-Si への TBmBJ + 24×24×24 grid で実験値再現
- Tran & Blaha (2009) PRL: TBmBJ汎関数の提案
- Aspnes & Studna (1983) PRB: Siの誘電関数実験値
