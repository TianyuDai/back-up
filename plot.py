import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline
from scipy.special import kn
import scipy.integrate as integrate
"""
data1 = np.loadtxt("llbt_l_pureglue_20GeV_muscale1.000000alpha0.100000.txt")
data2 = np.loadtxt("llbt_l_pureglue_20GeV_muscale1.000000alpha0.100000_check.txt")

plt.figure()

x1 = np.array(data1.T[0]); 
y1 = np.array(data1.T[1]); 
err1 = np.array(data1.T[2]); 
plt.errorbar(x1, y1, yerr=err1, fmt='o-', markersize=1, label='original', capsize=2, elinewidth=0.5, linewidth=0.5)

x2 = np.array(data2.T[0]); 
y2 = np.array(data2.T[1]); 
err2 = np.array(data2.T[2]); 
plt.errorbar(x2, y2, yerr=err2, fmt='o-', markersize=1, label='changed', capsize=2, elinewidth=0.5, linewidth=0.5)

plt.legend()
#plt.yscale("log")
plt.xlabel("E(GeV)")
plt.ylabel("dN/dE")
plt.title("Large angle + langevin, T=300MeV, $\\alpha_s = 0.01$, $E_0 = 20GeV$, $\\alpha_s^2t=0.3^2$")
#plt.savefig("alpha001_muperp_comparison_log_more.pdf")
plt.show()

"""

data1 = np.loadtxt("20GeV_muscale1.000000alpha0.300000omega0.100000_gluemedium_inel.txt")
data2 = np.loadtxt("20GeV_muscale1.000000alpha0.300000omega0.250000_gluemedium_inel.txt")
data3 = np.loadtxt("20GeV_muscale1.000000alpha0.300000omega0.500000_gluemedium_inel.txt")
data4 = np.loadtxt("20GeV_muscale1.000000alpha0.300000omega1.000000_gluemedium_inel.txt")
data5 = np.loadtxt("20GeV_muscale1.000000alpha0.300000omega2.000000_gluemedium_inel.txt")
data6 = np.loadtxt("20GeV_muscale1.000000alpha0.300000omega4.000000_gluemedium_inel.txt")


plt.figure()

x1 = np.array(data1.T[0]); 
y1 = np.array(data1.T[1]); 
err1 = np.array(data1.T[2]); 
plt.errorbar(x1, y1, yerr=err1, fmt='o-', markersize=1, label='$\omega/T = 0.1$', capsize=2, elinewidth=0.5, linewidth=0.5)

x2 = np.array(data2.T[0]); 
y2 = np.array(data2.T[1]); 
err2 = np.array(data2.T[2]); 
plt.errorbar(x2, y2, yerr=err2, fmt='o-', markersize=1, label='$\omega/T = 0.25$', capsize=2, elinewidth=0.5, linewidth=0.5)

x3 = np.array(data3.T[0]); 
y3 = np.array(data3.T[1]); 
err3 = np.array(data3.T[2]); 
plt.errorbar(x3, y3, yerr=err3, fmt='o-', markersize=1, label='$\omega/T = 0.5$', capsize=2, elinewidth=0.5, linewidth=0.5)

x4 = np.array(data4.T[0]); 
y4 = np.array(data4.T[1]); 
err4 = np.array(data4.T[2]); 
plt.errorbar(x4, y4, yerr=err4, fmt='o-', markersize=1, label='$\omega/T = 1$', capsize=2, elinewidth=0.5, linewidth=0.5)

x5 = np.array(data5.T[0]); 
y5 = np.array(data5.T[1]); 
err5 = np.array(data5.T[2]); 
plt.errorbar(x5, y5, yerr=err5, fmt='o-', markersize=1, label='$\omega/T = 2$', capsize=2, elinewidth=0.5, linewidth=0.5)

x6 = np.array(data6.T[0]); 
y6 = np.array(data6.T[1]); 
err6 = np.array(data6.T[2]); 
plt.errorbar(x6, y6, yerr=err6, fmt='o-', markersize=1, label='$\omega/T = 4$', capsize=2, elinewidth=0.5, linewidth=0.5)

plt.axvline(x=20, color='black', linewidth=0.5)
plt.axvline(x=10.8014, linewidth=0.5, linestyle='dashed')
plt.axvline(x=10.6758, color='green', linewidth=0.5, linestyle='dashed')
plt.axvline(x=10.1105, color='orange', linewidth=0.5, linestyle='dashed')
plt.axvline(x=10.0711, color='red', linewidth=0.5, linestyle='dashed')
plt.axvline(x=10.0037, color='purple', linewidth=0.5, linestyle='dashed')
plt.axvline(x=9.86571, color='brown', linewidth=0.5, linestyle='dashed')

plt.legend()
plt.yscale("log")
plt.xlabel("E(GeV)")
plt.ylabel("dN/dE")
plt.title("Large angle + langevin, T=300MeV, $\\alpha_s = 0.3$, $E_0 = 20GeV$, $\\alpha_s^2t=0.3^2$")
plt.savefig("alpha03_gluemedium_inel.pdf")
plt.show()

