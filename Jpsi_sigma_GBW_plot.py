import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime


# --- arquivo com dados experimentais ---
exp_data = pd.read_csv("sigma_gammap_jpsi.csv")

# --- todas as séries experimentais pertencem a Q² = 0.1 ---
Q2 = 0.0

# --- símbolos e cores diferentes para distinguir datasets experimentais ---
markers = ["o", "s", "^", "v", "D", "P", "X"]
colors = ["black", "gray", "orange", "green", "blue", "purple", "red"]

plt.figure(figsize=(8,6))

# --- plotar todos os datasets experimentais ---
unique_datasets = sorted(exp_data["dataset"].unique())
    
for i, dataset in enumerate(unique_datasets):
    subset = exp_data[exp_data["dataset"] == dataset]

    plt.errorbar(subset["W"], subset["Sigma"], yerr=subset["Error"],
                 fmt=markers[i % len(markers)],
                 color=colors[i % len(colors)],
                 ecolor="lightgray", capsize=3,
                 label=f"Exp. conjunto {dataset}")

# --- curvas teóricas para Q²=0.1 ---

    filename = f"Jpsi_sigma_Q2=0.0.csv"
    data = pd.read_csv(filename)

    W = data["W"]
    sigma_GLC = data["sigma_GLC"]
    sigma_BG = data["sigma_BG"]

    plt.plot(W, sigma_GLC, linestyle="-", color="red", linewidth=1.8, label=f"GLC Q²={Q2}")
    plt.plot(W, sigma_BG, linestyle="--", color="red", linewidth=1.2, label=f"BG Q²={Q2}")

    plt.text(W.iloc[-1]*1.05, sigma_GLC.iloc[-1], f"Q²={Q2}", 
             color="red", fontsize=9, va="center")


# --- ajustes de gráfico ---
plt.yscale("log")
plt.xlim(10, 350)
plt.xlabel(r"$W$ [GeV]")
plt.ylabel(r"$\sigma$ [nb]")
plt.title(r"$J/\psi$ produção exclusiva ($\gamma p \to J/\psi p$) — $Q^2 =0.0$")
plt.grid(True, which="both", ls="--", alpha=0.6)

# --- legenda sem duplicação ---
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc="lower right", fontsize=8)

plt.tight_layout()
# Salva com timestamp no nome
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
plt.savefig(f"plots/sigma_Jpsi_Q2=0.0_{timestamp}.png", dpi=300)

plt.show()
