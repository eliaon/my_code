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
exp_data = pd.read_csv("csv/sigma_gammap_jpsi.csv")

Q2 = 0

plt.figure(figsize=(8,6))

# --- mapa dataset -> experimento ---
exp_map = {
    0: "H1",
    1: "H1",
    2: "ALICE",
    3: "LHCb"
}

# estilos por experimento
exp_style = {
    "H1":   dict(color="blue",   marker="o"),
    "ALICE": dict(color="black",  marker="s"),
    "LHCb": dict(color="purple", marker="^")
}

# --- plot experimental ---
unique_datasets = sorted(exp_data["dataset"].unique())

for dataset in unique_datasets:

    subset = exp_data[exp_data["dataset"] == dataset]

    exp_name = exp_map.get(dataset, f"set {dataset}")
    style = exp_style.get(exp_name, dict(color="gray", marker="x"))

    plt.errorbar(
        subset["W"],
        subset["Sigma"],
        yerr=subset["Error"],
        fmt=style["marker"],
        color=style["color"],
        ecolor=style["color"],
        capsize=3,
        label=exp_name
    )

# --- curvas teóricas ---
try:
    filename = "csv/Jpsi_sigma_Q2=0.csv"
    data = pd.read_csv(filename)

    data = data.sort_values("W")

    W = data["W"]
    sigma_GLC = data["sigma_GLC"]
    sigma_BG = data["sigma_BG"]

    plt.plot(W, sigma_GLC, "-", color="red", lw=1.8, label=f"GLC")
    plt.plot(W, sigma_BG, "--", color="red", lw=1.2, label=f"BG")

    # anotação segura (não quebra layout)
    plt.text(
        0.72, 0.88,
        f"$Q^2={Q2}$",
        transform=plt.gca().transAxes,
        color="red",
        fontsize=9
    )

except Exception as e:
    print("Falha ao carregar curva teórica:", e)


# --- ajustes de gráfico ---
plt.yscale("log")
plt.xscale("log")
plt.xlim(300, 10000)
plt.ylim(1e1, 2e3)

plt.xlabel(r"$W$ [GeV]")
plt.ylabel(r"$\sigma$ [nb]")
plt.title(r"$J/\psi$ produção exclusiva ($\gamma p \to J/\psi p$)")

plt.grid(True, which="both", ls="--", alpha=0.6)

# legenda sem duplicação
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc="lower right", fontsize=8)

plt.tight_layout()

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
plt.savefig(f"plots/sigma_Jpsi_Q2=0_{timestamp}.png", dpi=300)

plt.show()

