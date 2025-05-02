
import math
import numpy as np
import matplotlib.pyplot as plt

def plot_river(sz, xs, zs, Hs):

    # 結果のプロット
    plt.plot(xs, zs, label="$z_b$")        # 河床高のプロット
    plt.plot(xs, Hs, label="$H$")          # 水位のプロット
    plt.xlabel("$x$ (m)")                   # x軸ラベル
    plt.ylabel("$z_b, H$ (m)")             # y軸ラベル
    plt.legend()                            # 凡例の表示
    plt.grid()                              # グリッドの表示
    plt.show()                              # グラフの表示

    # 計算結果の表形式出力
    print("|$x$<br>[m]|$z$<br>[m]|$H$<br>[m]|")
    print("|---|---|---|")
    for i in range(sz):
        print(f"|{xs[i]:.0f}|{zs[i]:.1f}|{Hs[i]:.3f}|")


def main(n, So, B, Q, h, z, dx, L):

    # 配列の初期化
    sz = L // dx + 1                    # 計算点の数
    xs = np.linspace(0, L, sz)         # x座標の配列
    zs = xs * So                        # 河床高の配列
    Hs = np.empty(sz)                   # 水位の配列

    # 計算に使用する定数の準備
    p = 10 / 3                          # Manning式のべき乗数
    qq = (Q / B)**2                     # 単位幅流量の2乗
    qq_2g = qq / 19.6                   # 2gで割った単位幅流量の2乗
    K = n * n * qq * dx / 2             # Manning式の係数

    # 下流端の境界条件
    H = z + h; Hs[0] = H

    # 逐次計算による水面形の計算
    cnst = H + qq_2g / (h * h) + K / math.pow(h, p)
    for i, z in zip(range(1, sz), zs[1:]):
        cnst -= z
        # 水深の初期推定値として直下の計算点の水深を当てる
        while True:
            vv_2g = qq_2g / (h * h)        # 速度水頭
            Sfdx_2 = K / math.pow(h, p)    # 損失水頭
            er = h + vv_2g - Sfdx_2 - cnst # 残差
            if abs(er) < 1e-4: break       # 収束判定
            # ニュートン法による水深の補正
            h -= er / (1 + (p * Sfdx_2 - 2 * vv_2g) / h)
        H = h + z; Hs[i] = H               # 水位の計算と保存
        cnst = H + vv_2g + Sfdx_2          # 次の計算点のための定数更新

    # 計算結果のプロット    
    plot_river(sz, xs, zs, Hs)


if __name__ == "__main__":

    # 計算条件の設定
    n = 0.02      # Manning の粗度係数
    So = 1 / 2000 # 河床勾配
    B = 700       # 川幅 [m]
    Q = 7000      # 流量 [cu.m/s]
    h = 5         # 下流端の水深 [m]
    z = 0         # 下流端の河床高 [m]
    dx = 500      # 計算距離刻み [m]
    L = 5000      # 流路長 [m]

    main(n, So, B, Q, h, z, dx, L)
