import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# --- arquivo com dados experimentais ---
# ainda n tenho dados experimentais para rapidez

# --- todas as séries experimentais pertencem a Q² = 0.1 ---
Q2 = 0

plt.figure(figsize=(8,6))

# --- curvas teóricas para Q²=0.1 ---
try:
    filename = f"csv/phi_d_sigma_dy_Q2=0.csv"
    data = pd.read_csv(filename)

    Y = data["y"]
    rap_GLC = data["d_sigma_dy_GLC"]
    rap_BG = data["d_sigma_dy_BG"]

    plt.plot(Y, rap_GLC, linestyle="-", color="red", linewidth=1.8, label=r"GLC \sqrt{{s}}=13 TeV")
    plt.plot(Y, rap_BG, linestyle="--", color="red", linewidth=1.2, label=r"BG \sqrt{{s}}=13 TeV")

    plt.text(Y.iloc[-1]*0.9, rap_GLC.iloc[-1]*1.1, f"Q²={Q2}", 
             color="red", fontsize=9, va="center")
except Exception as e:
    print(f"Aviso: não consegui carregar phi_d_sigma_dy_Q2=0.csv ({e})")
# --- ajustes de gráfico ---
#plt.xlim(2, 4.5)
#plt.ylim(0, 14)
plt.xlabel("Y")
plt.ylabel("dσ/dY [nb]")
plt.title("phi - Rapidez")
plt.grid(True)
plt.legend()
plt.savefig(f"plots/phi_rapidez_plot_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png", dpi=300)
plt.show()