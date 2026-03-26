import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

def load_set(xstr):
    ipsat = pd.read_csv(f"csv/sigma_dipolo_ipsat_x={xstr}.csv")
    return ipsat


fig, axes = plt.subplots(1,2, figsize=(12,5), sharey=True)

for ax, xstr in zip(axes, ["0.000100","0.010000"]):

    ipsat = load_set(xstr)

    ax.plot(ipsat['r'], ipsat['sigma_dipolo'], label='ipsat')

    ax.set_title(f'x={xstr}')
    ax.set_xlim(0, 2)
    ax.set_ylim(0,60)
    ax.grid()
    ax.legend()

axes[0].set_xlabel(r'$r$(fm)')
axes[1].set_xlabel(r'$r$(fm)')
axes[0].set_ylabel(r'$\sigma_{qq}$(mb)')
#axes[0].set_yscale('log')
#axes[1].set_yscale('log')
#axes[0].set_xscale('log')
#axes[1].set_xscale('log')
plt.tight_layout()
plt.savefig(f"plots/sigmaqq_ipsat_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.png")
plt.show()

