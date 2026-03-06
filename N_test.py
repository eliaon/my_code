import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

def load_set(xstr):
    ipsat = pd.read_csv(f"csv/N_ipsat_x={xstr}.csv")
    ipsat['x'] = xstr
    return ipsat

x_values = ["0.000100", "0.001000", "0.010000"]

plt.figure(figsize=(7,5))

for xstr in x_values:
    ipsat = load_set(xstr)
    plt.plot(ipsat['r'], ipsat['N_ipsat'], label=f'x={xstr}')

plt.xlim(8e-3,2e0)
plt.ylim(1e-3,2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$r$')
plt.ylabel(r'$N_p$')
plt.grid()
plt.legend()

plt.tight_layout()
plt.savefig(f"plots/N_x-ipsat_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.png")
plt.show()