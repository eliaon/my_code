import pandas as pd
import matplotlib.pyplot as plt

filename = "Jpsi_overlap_r.csv"
plot = pd.read_csv(filename)

r = plot["r"]
overlap_GLC = plot["overlap_GLC"]
overlap_BG = plot["overlap_BG"]

# Crie a figura ANTES de plotar
plt.figure(figsize=(8, 6))

plt.plot(r, overlap_GLC, linestyle="--", color="black", linewidth=1.8, label="GLC")
plt.plot(r, overlap_BG, linestyle="-", color="black", linewidth=1.2, label="BG")

plt.xlabel(r"$r$ [fm]")
plt.xlim(0.01, 1)
plt.ylim(0, 0.025)
plt.xscale("log")
plt.ylabel(r"Overlap $r \Psi^*_V \Psi_{\gamma}$")
plt.title(r"$J/\psi$ sobreposição de funções de onda ($\gamma^*p \to J/\psi p$)")
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend(loc="upper right", fontsize=8)
plt.tight_layout()
plt.savefig("Jpsi_overlap_plot.png", dpi=300)
plt.show()