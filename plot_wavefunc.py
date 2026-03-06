import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

filename = "csv/Jpsi_wavefunctions_GLC.csv"
data = pd.read_csv(filename)
r = data["r"]
phi_T_BG = data["phi_T"]
phi_L_BG = data["phi_L"]
z = data["z"]

r_vals = np.sort(r.unique())
z_vals = np.sort(z.unique())

R, Z = np.meshgrid(r_vals, z_vals)

phi_T = phi_T_BG.values.reshape(len(z_vals), len(r_vals))
phi_L = phi_L_BG.values.reshape(len(z_vals), len(r_vals))

fig = plt.figure(figsize=(14, 6))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(R, Z, phi_T, cmap='viridis', linewidth=0, antialiased=True)

ax.set_title('J/ψ Transverse Wavefunction (GLC Model)')
ax.set_xlabel('r [GeV$^{-1}$]')
#plt.xlim(0, 1)
#plt.ylim(-1, 1)
ax.set_ylabel('z')
ax.set_zlabel(r'$\phi_T$ [GeV$^{-1}$]')
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
plt.tight_layout()
plt.savefig(f"plots/Jpsi_wavefunction_T_GLC_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png", dpi=300)
plt.show()


