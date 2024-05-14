import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline, splrep, splev
# import scienceplots
plt.style.use('science')

# QFT and quasiparticle QFT results for Z
rs_qft = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
z_qft = [1.0, 0.9224, 0.8608210502584686, 0.749144526156203, 0.6622891781023675, 0.5861045365581065, 0.5208831611860456, 0.46342605505136064, 0.37562045192664967]
z_qft_err = [1e-8, 0.0002, 0.0028199354623565345, 0.00832375521694881, 0.010795966185016277, 0.009732720423881042, 0.0100925784804457435, 0.01569005557909372, 0.027285334913524847]
z_qft_quasi= [1.0, 0.9224, 0.8611633434933357, 0.7496024658544531, 0.6709968875211331, 0.6040765444587273, 0.5492968494379177, 0.5025532635810275, 0.457480762325804]
z_qft_quasi_err = [1e-8, 0.0002, 0.001040187306472799, 0.01779348428847196, 0.010931222054177377, 0.015168873975929115, 0.023233316150618702, 0.057910434448675205, 0.10363210442228994]


# Kristjan & Chen 2022, T=T_F/25
rs_old = [1.0, 2.0, 3.0, 4.0]
z_old = [0.8725, 0.7984, 0.7219, 0.6571]
z_error_old = [0.0002, 0.0002, 0.0002, 0.0002]

# RPA
rs_rs = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
z_rs = [1.0, 0.922378, 0.8601, 0.7642, 0.6927, 0.6367, 0.5913, 0.5535, 0.52244, 0.49496]

# BF-RMC
rs_qmc = [1.0, 2.0, 3.99, 5.0, 10.0]
z_rmc_bf = [0.84, 0.77, 0.64, 0.58, 0.40]
z_rmc_bf_err = [0.02, 0.01, 0.01, 0.01, 0.01]

# SJ-VMC
rs_qmc = [1.0, 2.0, 3.99, 5.0, 10.0]
z_vmc_sj = [0.894, 0.82, 0.69, 0.61, 0.45]
z_vmc_sj_err = [0.009, 0.01, 0.01, 0.02, 0.01]

# BF-VMC
rs_qmc = [1.0, 2.0, 3.99, 5.0, 10.0]
z_vmc_bf = [0.86, 0.78, 0.65, 0.59, 0.41]
z_vmc_bf_err = [0.01, 0.01, 0.01, 0.02, 0.01]

# Define the data for plotting

# Plot the data
fig, ax = plt.subplots(figsize=(8,6))

# # Fit a polynomial curve to the QFT data
# coeffs_qft = np.polyfit(rs_qft, z_qft, 2)
# curve_qft = np.polyval(coeffs_qft, rs_qft)

# # Fit a polynomial curve to the Quasi-QFT data
# coeffs_quasi = np.polyfit(rs_qft, z_qft_quasi, 2)
# curve_quasi = np.polyval(coeffs_quasi, rs_qft)

xnew = np.linspace(0, 8, 300)

# Fit a smooth curve through all the RPA datapoints
tck = splrep(rs_rs, z_rs, s=0)
curve_rpa = splev(xnew, tck)

# Fit a spline curve to the QFT data
spl_qft = UnivariateSpline(rs_qft, z_qft, w=1/np.array(z_qft_err), k=2, s=None)
curve_qft = spl_qft(xnew)

# Fit a spline curve to the Quasi-QFT data
spl_quasi = UnivariateSpline(rs_qft, z_qft_quasi, w=1/np.array(z_qft_quasi_err), k=2, s=None)
curve_quasi = spl_quasi(xnew)

markersize = 4
capsize = 3
colors = ['blue', 'red', 'green', 'black', 'purple', 'brown', 'orange', 'gray']

def plot(x, y, y_err, fmt, markersize, label, color, markerfacecolor):
    plt.errorbar(x, y, yerr=y_err, fmt=fmt, markersize=markersize, label=label, color=color,
                 capsize=3, elinewidth=1.5, markerfacecolor=markerfacecolor)

# RPA
# plt.plot(rs_rs, z_rs, '-', linewidth=2, markersize=markersize, label='RPA', color=colors[3])
plt.plot(xnew, curve_rpa, '-', linewidth=1, markersize=markersize, label='RPA', color=colors[3])

# BF-RMC
plot(rs_qmc, z_rmc_bf, z_rmc_bf_err, '*', 4, 'BF-RMC', colors[4], 'none')

# SJ-VMC
plot(rs_qmc, z_vmc_sj, z_vmc_sj_err, 'D', 4, 'SJ-VMC', colors[5], 'none')

# BF-VMC
plot(rs_qmc, z_vmc_bf, z_vmc_bf_err, 'P', 4, 'BF-VMC', colors[6], 'none')

# Kristjan & Chen 2022, T=T_F/25
plot(rs_old, z_old, z_error_old, 's', 4, 'VDMC, $T=T_F/25$', colors[2], 'none')

# Experimental value
plot(3.99, 0.58, 0.07, 'o', 6, 'Experiment', 'gray', 'none')

# QFT and quasiparticle QFT results for Z
# plot(rs_qft, z_qft_quasi, z_qft_quasi_err, '^', markersize, '$G_RW_0$', colors[1], colors[1])
# plt.plot(xnew, curve_quasi, '--', color=colors[1], linewidth=1, markersize=markersize)

# z_qft_err = [e*3.0 for e in z_qft_err]
plot(rs_qft, z_qft, z_qft_err, 'o', markersize, '$G_0W_0$', colors[0], colors[0])
plt.plot(xnew, curve_qft, '--', color=colors[0], linewidth=1, markersize=markersize)

# Set the plot title and axes labels
# plt.title('Z-Factor versus R_s', fontsize=16)
plt.xlabel('$r_s$', fontsize=18)
plt.ylabel('$Z$', fontsize=18)

# Set the axis limits
plt.xlim(0, 8.5)
# plt.ylim(0.3, 1)

# Add a grid
# plt.grid()

# set the tick label size
ax.tick_params(axis='both', which='major', labelsize=16)

# Add a legend
plt.legend(loc='best', fontsize=16, ncol=2)
plt.tight_layout()

# Save the plot as a PDF file
plt.savefig('z_factor.pdf')

# Show the plot
plt.show()
