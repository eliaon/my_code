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
#exp_data = pd.read_csv("phi_gamap_paper_all.csv")

# --- todas as séries experimentais pertencem a Q² = 0.1 ---
Q2 = 0

# --- símbolos e cores diferentes para distinguir datasets experimentais ---
markers = ["o", "s", "^", "v", "D", "P", "X"]
colors = ["black", "gray", "orange", "green", "blue", "purple", "red"]

plt.figure(figsize=(8,6))

# --- plotar todos os datasets experimentais ---
#unique_datasets = sorted(exp_data["dataset"].unique())

#for i, dataset in enumerate(unique_datasets):
#    subset = exp_data[exp_data["dataset"] == dataset]

#    plt.errorbar(subset["W"], subset["Sigma"], yerr=subset["Error"],
#                 fmt=markers[i % len(markers)],
#                 color=colors[i % len(colors)],
#                 ecolor="lightgray", capsize=3,
#                 label=f"Exp. conjunto {dataset}")

# --- curvas teóricas para Q²=0.1 ---
try:
    filename = f"csv/phi_sigma_Q2={Q2}.csv"
    data = pd.read_csv(filename)

    W = data["W"]
    sigma_GLC = data["sigma_GLC"]
    sigma_BG = data["sigma_BG"]

    plt.plot(W, sigma_GLC, linestyle="-", color="red", linewidth=1.8, label=f"GLC")
    plt.plot(W, sigma_BG, linestyle="--", color="red", linewidth=1.2, label=f"BG")

    plt.annotate(
    f"$Q^2={Q2}$",
    xy=(W.iloc[-1], sigma_GLC.iloc[-1]),
    xytext=(10,0),
    textcoords="offset points",
    color="red"
)

except Exception as e:
    print(f"Aviso: não consegui carregar sigma_Q2={Q2}.csv ({e})")

#plotar dado experimental de H1 para Q2=2
plt.errorbar(70, 960.0, yerr=280.0, fmt="o", color="blue", ecolor="blue", capsize=3, label="H1")

# --- ajustes de gráfico ---
plt.yscale("log")
plt.xscale("log")
plt.xlim(30, 1000)
plt.ylim(10, 5000)
plt.xlabel(r"$W$ [GeV]")
plt.ylabel(r"$\sigma$ [nb]")
plt.title(rf"$\phi$ produção exclusiva ($\gamma p \to \phi p$) — $Q^2={Q2}$")
plt.grid(True, which="both", ls="--", alpha=0.6)

# --- legenda sem duplicação ---
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc="lower right", fontsize=8)

plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
plt.savefig(f"plots/sigma_phi_Q2={Q2}_{timestamp}.png", dpi=600, bbox_inches='tight')
plt.show()
