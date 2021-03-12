import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from  scipy.optimize import curve_fit

def pol_1(x, a, b):
    y = a + b*x
    return y

def pol_3(x,a,b,c,d):
    return a + b*x + c*x*x + d*x*x*x

nsteps = [10, 20, 50, 100, 200, 1000]


df = pd.read_csv("./data/gauss_wf_wp_129_1_10_0_3_30_6_1.txt", sep="\t")
print(df.head())
df.dtypes
df.describe()
df.loc[-1, ['tavg_state_energy']]

list_nsteps = [10, 20, 50, 100, 200]
averx_EU = [ 6.5953679520216994e+01, 6.5999999999999858e+01, 6.5999999999999986e+01, 6.6000000000000000e+01, 6.6000000000000000e+01]
averx_UCM = [ 3.1597744822801875e+01 , 3.3372174193127833e+01, 3.8695347803555599e+01, 4.7562401808328651e+01, 6.0940384627076881e+01]
list_nsteps_SSM = [10, 20, 50, 100, 200, 1000]
averx_SSM = [ 3.1604928333763809e+01, 3.3388175710852977e+01, 3.8737908489339816e+01, 4.7649005852721366e+01, 6.0955560903005633e+01, 5.4514916844856444e+01]
averx_SSM_short = [ 3.1604928333763809e+01, 3.3388175710852977e+01, 3.8737908489339816e+01, 4.7649005852721366e+01, 6.0955560903005633e+01]


tau = df["tau"].to_numpy()
averx = df["averx"].to_numpy()
np.transpose(tau)
# tau.tolist()
averx.dtype
tau.dtype
popt3,_ = curve_fit(pol_3, tau, averx)
popt1,_ = curve_fit(pol_1,  tau, averx)

# Plot for E(tau) in EU
sns.set_theme()
plt.figure(figsize=(8, 5))
plt.xlabel(r"$time$")
plt.ylabel(r'$averx (N)$')
plt.title(r'$averx$')
plt.plot(tau, averx)
# plt.plot(tau, pfit, c='firebrick', marker='o', label=r'A = 0.8')
# plt.plot(tau, pol_1(tau,*popt1), c='firebrick', )
# plt.show()
plt.savefig('./plots/Euler_averx(time).pdf')
#
sns.set_theme()
plt.figure(figsize=(8, 5))
plt.xlabel(r"$time step$")
plt.ylabel(r'$averx (N)$')
plt.title(r'$averx(time step) for Euler$')
plt.plot(list_nsteps, averx_EU)
plt.savefig('./plots/Euler_averx(nstep).pdf')

sns.set_theme()
plt.figure(figsize=(8, 5))
plt.xlabel(r"$time step$")
plt.ylabel(r'$averx (N)$')
plt.title(r'$averx(time step) for UCM and SSM$')
plt.plot(list_nsteps, averx_UCM, marker='+', c="r", label="UCM")
plt.plot(list_nsteps, averx_SSM_short, marker='+',c="g", label="SSM")
plt.legend()
# plt.show()
plt.savefig('./plots/UCM_SSM_averx(nstep).pdf')

sns.set_theme()
plt.figure(figsize=(8, 5))
plt.xlabel(r"$time step$")
plt.ylabel(r'$averx (N)$')
plt.title(r'$averx(time step) for Euler$')
plt.plot(list_nsteps_SSM, averx_SSM, marker='+', c='g')
plt.savefig('./plots/SSM_averx(nstep).pdf')
