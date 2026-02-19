import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib as mpl

mpl.rcParams.update({
    "font.size": 18,          # fonte base
    "axes.titlesize": 22,     # título
    "axes.labelsize": 20,     # labels dos eixos
    "legend.fontsize": 18,    # legenda
    "xtick.labelsize": 16,    
    "ytick.labelsize": 16,
})

# --- arquivo com dados experimentais ---
# ainda n tenho dados experimentais para rapidez

# --- todas as séries experimentais pertencem a Q² = 0.1 ---
Q2 = 0
s = 5.26e+03 # GeV
plt.figure(figsize=(8,6))
meson = "phi"
# --- curvas teóricas para Q²=0.1 ---
try:
    filename = f"csv/{meson}_rapidez_5.26e+03GeV.csv"
    data = pd.read_csv(filename)

    Y = data["y"]
    rap_GLC = data["d_sigma_dy_GLC"]
    rap_BG = data["d_sigma_dy_BG"]
    plt.plot(Y, rap_GLC, linestyle="-", color="red",
         linewidth=1.8,
         label=rf"GLC")

    plt.plot(Y, rap_BG, linestyle="--", color="red",
         linewidth=1.2,
         label=rf"BG")

    plt.text(Y.iloc[-1]*0.9,
         rap_GLC.iloc[-1]*1.1,
         "GLC",
         color="red",
         fontsize=9,
         va="center")

except Exception as e:
    print(f"Aviso: não consegui carregar {meson}_rapidez_5.26e+03GeV.csv ({e})")
# --- ajustes de gráfico ---
plt.xlim(-8, 8)
#plt.ylim(0, 140)
plt.xlabel("Y")
plt.ylabel("dσ/dY [nb]")
if meson == "Jpsi":
    meson = "psi"
else:
    meson = meson
plt.title(f" Distribuição de rapidez do $\{meson}$ em pp a {s/1e+03:.2f} TeV")
plt.grid(True)
plt.legend()
plt.savefig(f"plots/{meson}_rapidez_5,26TeV_plot_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png", dpi=300)
plt.show()