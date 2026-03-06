import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

def load_set(xstr):
    ipsat = pd.read_csv(f"csv/N_ipsat_x={xstr}.csv")
    gbw  = pd.read_csv(f"csv/N_GBW_x={xstr}.csv")
    iim  = pd.read_csv(f"csv/N_IIM_x={xstr}.csv")
    return ipsat, gbw, iim


fig, axes = plt.subplots(1,2, figsize=(12,5), sharey=True)

for ax, xstr in zip(axes, ["1e-4","1e-2"]):

    ipsat, gbw, iim = load_set(xstr)

    ax.plot(ipsat['r2'], ipsat['N_ipsat'], label='ipsat')
    ax.plot(gbw['r2'], gbw['N_GBW'], '--', label='GBW')
    ax.plot(iim['r2'], iim['N_IIM'], ':', label='IIM')

    ax.set_title(f'x={xstr}')
    ax.set_xlim(0,100)
    ax.set_ylim(0,1.2)
    ax.grid()
    ax.legend()

axes[0].set_xlabel(r'$r^2$')
axes[1].set_xlabel(r'$r^2$')
axes[0].set_ylabel(r'$N_p$')

plt.tight_layout()
plt.savefig(f"plots/N_gbw_iim_ipsat_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.png")
plt.show()

