import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib as mpl

mpl.rcParams.update({
    "font.size": 18,          # fonte base
    "axes.titlesize": 22,     # título
    "axes.labelsize": 20,     # labels dos eixos
    "legend.fontsize": 18,    # legenda
    "xtick.labelsize": 16,    
    "ytick.labelsize": 16,
})
filename = "Jpsi_overlap_r.csv"
plot = pd.read_csv(filename)

r = plot["r"]
overlap_GLC = plot["overlap_GLC"]
overlap_BG = plot["overlap_BG"]

# Crie a figura ANTES de plotar
plt.figure(figsize=(10, 6))

plt.plot(r, overlap_GLC, linestyle="--", color="red", linewidth=1.8, label="GLC")
plt.plot(r, overlap_BG, linestyle="-", color="black", linewidth=1.2, label="BG")

plt.xlabel(r"$r$ [fm]")
plt.xlim(0.001, 1)
plt.ylim(0, 0.025)
plt.xscale("log")
plt.ylabel(r"Overlap $r \Psi_V \Psi_{\gamma}$")
plt.title(r"Sobreposição de funções de onda ($\gamma p \to J/\psi p$)")
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend(loc="upper right", fontsize=8)
plt.tight_layout()
plt.savefig(f"plots/Jpsi_overlap_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.png", dpi=300)
plt.show()