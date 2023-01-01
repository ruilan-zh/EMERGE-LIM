import numpy as np

dir1 = "emerge_data"
halo_types = ["sum", "cent", "sat"]

h = 0.6736
boxsize = 250 /h

logMmin = 9
dlogM = 0.3
logMmax = 15
nbins = (logMmax - logMmin)/dlogM
mass_bins = np.arange(logMmin, logMmax, dlogM)

for halo in halo_types:
    sfr_sums = [0]*len(mass_bins)
    sfr_counts = [0]*len(mass_bins)
    tot_counts = [0]*len(mass_bins)

    data = np.loadtxt(f"{dir1}/mass-sfr-{halo}.txt", unpack=True)
   # sums = np.loadtxt(f"{dir1}/mass-sfr-sum.txt", unpack=True)
   # sats = np.loadtxt(f"{dir1}/mass-sfr-sat.txt", unpack=True)

    mass1, sfr1 = data[0], data[-1]
   # mass_cent, sfr_cent = cents[0], cents[-1]
   # mass_sat, sfr_sat = sats[0], sats[-1]

    mask = (np.isinf(sfr1) == 0) & (np.isnan(sfr1) == 0) 
   # mask_sum = (np.isinf(sfr_sum) == 0) & (np.isnan(sfr_sum) == 0)  
   # mask_sat = (np.isinf(sfr_sat) == 0) & (np.isnan(sfr_sat) == 0) 


    masses  = mass1[mask]
    sfrs = sfr1[mask]
   # masses_cent  = mass_cent[mask_cent]
   # sfrs_cent = sfr_cent[mask_cent]
   # masses_sat  = mass_sat[mask_sat]
   # sfrs_sat = sfr_sat[mask_sat]

    for ihalo, mass in enumerate(masses):
        ibin = int((mass - logMmin)/dlogM)
        if ibin < nbins:
            tot_counts[ibin] += 1
            if sfrs[ihalo] != 0:
                sfr_sums[ibin] += 10**sfrs[ihalo]
                sfr_counts[ibin] += 1

    f = open(f"{dir1}/emerge_sfrd_bins_{halo}.txt", "w")
    for ibin, m in enumerate(mass_bins):
        print(m, sfr_sums[ibin]/boxsize**3, file=f)
        
    f.close()
