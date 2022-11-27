import numpy as np

boxsize = 250

logMmin = 8
dlogM = 0.2
logMmax = 15
nbins = (logMmax - logMmin)/dlogM
mass_bins = np.arange(logMmin, 15, dlogM)

sfr_sums = [0]*len(mass_bins)
sfr_counts = [0]*len(mass_bins)
tot_counts = [0]*len(mass_bins)

#mass, sfr = np.loadtxt("mass-sfr.txt", unpack=True)
mass1, sfr1 = np.loadtxt("mass-sfr-sum.txt", unpack=True)

#mask = (np.isinf(sfr) == 0) & (np.isnan(sfr) == 0) & (mass < 11)
mask1 = (np.isinf(sfr1) == 0) & (np.isnan(sfr1) == 0)

#masses = np.concatenate((mass[mask], mass1[mask1]) )
#sfrs = np.concatenate((sfr[mask], sfr1[mask1])) 

masses = mass1[mask1]
sfrs = sfr1[mask1]


for ihalo, mass in enumerate(masses):
    ibin = int((mass - logMmin)/dlogM)
    if ibin < nbins:
        tot_counts[ibin] += 1
        if sfrs[ihalo] != 0:
            sfr_sums[ibin] += 10**sfrs[ihalo]
            sfr_counts[ibin] += 1

f = open("emerge_sfrd_bins.txt", "w")
for ibin, m in enumerate(mass_bins):
    print(m, sfr_sums[ibin]/boxsize**3, file=f)
    
f.close()
