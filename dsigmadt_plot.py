import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

def plot_dsigma_dt(W):
    
    plot_w = pd.read_csv(f"csv/Jpsi_dsigma_dt_W={W}GeV.csv")
    if W == "100":
        plt.plot(plot_w['t'], plot_w['dsigma_dt_GLC']*5, label='GLC')
        plt.plot(plot_w['t'], plot_w['dsigma_dt_BG']*5, label='BG', linestyle='--')
    else:
        plt.plot(plot_w['t'], plot_w['dsigma_dt_GLC'], label=f'GLC (W={W} GeV)')
        plt.plot(plot_w['t'], plot_w['dsigma_dt_BG'], label=f'BG (W={W} GeV)', linestyle='--')
    


w_values = ["75", "100"]
plt.figure(figsize=(7,5))
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

for W_str in w_values:
    plot_dsigma_dt(W_str)

plt.xlim(0, 2.5)
plt.ylim(1e-3, 1e4)
plt.xscale('linear')
plt.yscale('log')
plt.xlabel(r'$|t|$ (GeV$^{-2}$)')
plt.ylabel(r'$d\sigma/dt$ (GeV$^{-2}$)')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(f"plots/dsigma_dt_comparison_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.png")
plt.show()
