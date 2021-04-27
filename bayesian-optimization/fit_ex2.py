import matplotlib.pyplot as plt 
import numpy as np 

data = np.loadtxt('bounds_ex2.dat')

fit = np.polyfit(data[:, 0], data[:, 1], 2)

polyfit = np.poly1d(fit)

x = np.linspace(2.5, 3.5, 1000)

fit_chi2 = polyfit(x)

# Find minimum
min_chi2 = 1000
min_chi2_energy = 0
for x_ in x:
    if polyfit(x_) < min_chi2:
        min_chi2 = polyfit(x_)
        min_chi2_energy = x_

# Low bound
low_low = -0.01
low_high = min_chi2_energy
while (low_high - low_low) > 1e-5:
    high_mid = (low_high + low_low)/2.
    current_chi2 = polyfit(high_mid)
    if current_chi2 < min_chi2 + 1.5:
        low_high = high_mid
    else:
        low_low = high_mid

# High bound
high_low = min_chi2_energy
high_high = 0.05
while (high_high - high_low) > 1e-5:
    high_mid = (high_high + high_low)/2.
    current_chi2 = polyfit(high_mid)
    if current_chi2 > min_chi2 + 1.5:
        high_high = high_mid
    else:
        high_low = high_mid

print(min_chi2_energy)
print(low_low, low_high)
print(high_high)

plt.scatter(data[:, 0], data[:, 1])
plt.plot(x, polyfit(x), 'r--')

plt.show()