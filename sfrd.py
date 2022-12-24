import numpy as np

cent = np.loadtxt("L250_N1290/mass-sfr.txt", unpack=True)
sfrsum = np.loadtxt("L250_N1290/mass-sfr-sum.txt", unpack=True)

mass,sfr = [cent[0], cent[-1]]
mass1, sfr1 = [sfrsum[0], sfrsum[-1]]

mask = (np.isinf(sfr) == 0) & (np.isnan(sfr) == 0)   #& (sfr > 1)# & (sfr < 1)
mask1 = (np.isinf(sfr1) == 0) & (np.isnan(sfr1) == 0)  & (sfr1>1)# & (sfr1<1)  


sfr_sum = np.sum(10**sfr[mask])
sfr_sum1 = np.sum(10**sfr1[mask1])
sfr_sum = sfr_sum1
boxsize_mpc_h = 250
boxsize_mpc = boxsize_mpc_h /0.67
#boxsize_mpc = 1000
volume = boxsize_mpc**3
sfrd = sfr_sum/volume
print(sfr_sum)
sfrd = np.log10(sfr_sum/volume)
print(sfrd)
