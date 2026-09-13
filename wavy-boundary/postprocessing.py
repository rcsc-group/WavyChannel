import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sc

plt.rcParams['text.usetex'] = True

eps = 0.01
sep = 0.2

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,6))

n_points = 3500
t, l = [], []
for i in range(n_points):
    time_point = i * sep
    
    try:
        interface = pd.read_csv(f"interface/interface-{time_point:09.3f}.dat", sep = " ", header=None)#.iloc[::2]
    except FileNotFoundError: # Reached the end of stored times
        print(f"Interface at t={time_point:09.3f} not found")
        break

    # Extracting interface x-position from bottom of computational domain
    l_val = interface.iloc[interface.idxmin()[1]][0]
    if l_val > 0.9997/eps: # Avoiding noise near end of channel
        continue

    t.append(time_point)
    l.append(l_val)

t = np.array(t)
l = np.array(l)
late = l>50 

ax1.axhline(50, color="k", linestyle="--")
p = ax1.plot(t, l, label="DNS Data")
cur_color = np.array(mpl.colors.to_rgb(p[-1].get_color())) * 1

def square_root(t, a, b):
    return np.sqrt(a*t + b)

# Fitting to trace using nonlinear least squares
res, pcov = sc.optimize.curve_fit(square_root, t[late], l[late], bounds=((0.,-np.inf), (np.inf,np.inf)), maxfev=10000)
per = np.sqrt(np.diag(pcov))
print(f"Fitted growth rate: {res[0]:10.2e}")
print(f"Fit uncertainty:    {per[0]:10.2e}")

ax1.plot(t, square_root(t, *res), color=cur_color, linestyle=":", linewidth=6, label="Fitted trace")  # Plotting the fitted curve
ax2.plot(t[late], (l-square_root(t, *res))[late]/eps)                                                 # Plotting the residuals

lim = np.max(np.abs(l-square_root(t, *res))[late])/eps*1.5
ax2.set_ylim(-lim, lim)

ax1.set_xlabel("Time $t$")
ax1.set_ylabel("Interface position $\\ell$")
ax1.legend()

ax2.set_xlabel("Time $t$")
ax2.set_ylabel("Absolute error $/\\:\\varepsilon$")

fig.tight_layout()
fig.savefig("img-growth-rate.pdf")
print("Figure saved in 'img-growth-rate.pdf'")