import pandas as pd
import matplotlib.pyplot as plt

filename = "csv/N_gbw_r.csv"
data = pd.read_csv(filename)
plt.plot(data['r'], data['N_gbw'], label='N_gbw')
plt.xlabel('$r²$')
plt.xlim(0, 100)
plt.ylabel('N_cgc')
plt.title(filename)
plt.legend()
plt.grid()
plt.savefig('plots/N_gbw_plot.png')
plt.show()