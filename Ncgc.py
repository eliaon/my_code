import pandas as pd
import matplotlib.pyplot as plt

filename = "csv/N_GBW.csv"
data = pd.read_csv(filename)
plt.plot(data['r²'], data['N_GBW'], label='N_GBW')
plt.xlabel('$r²$')
plt.xlim(0, 100)
plt.ylabel('N_GBW')
plt.title(filename)
plt.legend()
plt.grid()
plt.savefig('plots/N_gbw_plot.png')
plt.show()