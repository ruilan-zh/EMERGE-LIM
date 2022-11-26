#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "constant.h"
#include "para.h"
#include "proto.h"

#include "spline.h"

double E_z(Cosmology cosm, double z);
double omega_m(Cosmology cosm, double z);
double rho_crit(Cosmology cosm, double z);
double rho_m(Cosmology cosm, double z);

double t_dyn(Cosmology cosm, double z);
double rvir_from_mvir(Cosmology cosm, double mvir_sol_h, double z);
double my_t2z(double t_yr, Cosmology cosm);
double my_z2t(double z, Cosmology cosm);
double rho_nfw(Cosmology cosm, double conc, double m_vir_sol_h, double r_vir, double r);
double baryon_conversion_efficiency(double M, double z);

double modelKlypin16(double M, double z, char *mdef);

double mu(double conc);
double modelKlypin16(double M, double z, char *mdef);
double lin_interp(double x_arr[], double y_arr[], double z, int ibin);

double *c_pseudo_interp_table(int n_table, double table_x[], double table_y[]);
double c_pseudo(Cosmology cosm, double ci, double zi, double zf, double *cubic_set, int n_table, double table_x[], double table_y[]);

double M_pseudo(Cosmology cosm, double Mi, double zi, double zf, double * cubic_set, int n_table, double table_x[], double table_y[], double mass_bins[], double z_bins[], double **conc_arr);

double c_change_mass_def(Cosmology cosm, double ci, double Delta_i, double Delta_f, double zi, double zf, double *cubic_set, int n_table, double table_x[], double table_y[]);

double M_dyn(double masses[], double redshifts[], int index, double z1, int equals);

double tau(Cosmology cosm, double M_star_msun_h, double z) ;
double mhalo2mstar(Cosmology cosm, double M_halo_msun_h, double z);
double conc_ishiyama(double logM, double z, double mass_bins[], double z_bins[], double **conc_arr);

double sfr2mstar(Cosmology cosm, double sfr, double z);

double M_pinocchio2vir(Cosmology cosm, double M_pinocchio, double z, double *cubic_set, int n_table, double table_x[], double table_y[], double mass_bins[], double z_bins[], double **conc_arr);

int main(int argc, char **argv)
{
	clock_t start, end;
	double cpu_time_used;

	time_t now;
	time(&now);
	printf("%s\n", ctime(&now));

	Cosmology cosm;
	cosm.omega = 0.3153;
	cosm.lambda = 0.6847;
	cosm.omegab = 0.0493;
	cosm.hubble = 0.6736;

	typedef struct
	{
	unsigned long long int id; //8 bytes
	double M; //8 bytes
	double q[3], x[3], v[3]; //8*9 = 72 bytes
	int N, pad; //4*2 = 8 bytes
	} cat;
	cat catdata;

	typedef struct
	{
		unsigned long long int name; //group ID
		int nick, ll, mw, mass, mam; // index within tree, linking list, merged with, mass at merger, mass of main halo is merges with
		double zme, zpe, zap; // merger redshift, peak collapse redshift, halo overtakes min mass redshift
	} hist;
	hist histdata;


	int irun = 22;
//	double pmass = 1.27346 * pow(10,9); // N1024 particle mass [Msun/h]
//	double pmass = 3.77321 * pow(10,8); // N1536 particle mass [Msun/h]
//	double pmass = 2.55136 * pow(10,8); // N1750 particle mass [Msun/h]
//	double pmass = 1.34544 * pow(10,7); // N1400 75Mpc/h

//	double pmass = 2.59293 * 1e7; // N1500 100 Mpc/h
//	double pmass = 3.3383 * 1e8; // N1600 250 Mpc/h (z=0.4)

	double pmass = 1.59183 * 1e8; // N4096 500 Mpc/h (z=1.5)

//	double pmass = 6.36966 * 1e8; // N1290 250 Mpc/h (z=1.47)


//	double min_part = 10;
//	double min_part = 8;
	double min_part = 32;

	double min_hmass = pmass*min_part; 
	printf("Min halo mass: %lf\n", log10(min_hmass));

	double zout = 1.5;
	
	double tdyn = t_dyn(cosm, zout);
	printf("tdyn: %g years\n", tdyn);

	double tout = my_z2t(zout, cosm);

	double t1 = tout - tdyn;
	printf("tout - tdyn: %g years\n", t1);

	double z1 = my_t2z(t1, cosm);
	printf("z(tout - tdyn): %lf\n", z1);

	char runname[16], filename[256];

//	sprintf(runname, "r000%d", irun);
	sprintf(runname, "test32");


//	char cat_dir[128] = "/mnt/data_cat5/rlzhang/Pinocchio/test/output/tests22/1.47/N1750";
	char cat_dir[128] = "/mnt/data_cat4/moriwaki";

	FILE *fp_cat;
	char fname_cat[64];

	sprintf(fname_cat, "%s/pinocchio.%.1f000.%s.catalog.out", cat_dir, zout, runname);
	printf("Opening file %s\n",fname_cat);

	fp_cat = fopen(fname_cat, "r");
	if (fp_cat==NULL)
	{
		fprintf(stderr, "cannot open file: %s\n", fname_cat);
		exit(1);
	}


	int dummy_cat, ninhead, NSlices_cat, nproc;

	int u_cat = fread(&ninhead, sizeof(int),1,fp_cat);
	ninhead /= sizeof(int);
	u_cat = fread(&nproc, sizeof(int), 1, fp_cat);
	if (ninhead==2)
		u_cat = fread(&NSlices_cat, sizeof(int), 1, fp_cat);
	u_cat = fread(&dummy_cat, sizeof(int), 1, fp_cat);
	printf("The file has been written by %d tasks\n", nproc);
	if (ninhead==2)
		printf("The box has been fragmented in %d slices\n",NSlices_cat);
	else
		printf("This old version catalog does not give the number of slices used\n");

	
	int arrlen = 100000;

	unsigned long long int *groupids; //8 bytes
	groupids = (unsigned long long int*) malloc(arrlen*sizeof(unsigned long long int));
	double *mass_cat;
	mass_cat = (double*) malloc(arrlen*sizeof(double));
	double **pos;
	pos = (double**) malloc(3 * sizeof(double*));
	for (int i = 0; i < 3; i++)
	{
		pos[i] = (double*) malloc(arrlen*sizeof(double));
	}


	int nblocks = 0;
	int ngood, igood;

	int n = 0;
	int nhalo = 0;

	while (1)
	{
		if (!fread(&dummy_cat, sizeof(int), 1, fp_cat))// determines number of bytes of cat; if 0 bytes then at end of file and !0 = 1 
			break;
			
				
		++nblocks;
		u_cat = fread(&ngood,sizeof(int),1,fp_cat);
		u_cat = fread(&dummy_cat,sizeof(int),1,fp_cat); 
		n+=ngood;
		//printf("  found %d halos in block N. %d\n",ngood,nblocks);

		if (n > arrlen)
		{
			arrlen = n;
			groupids = realloc(groupids, arrlen * sizeof(unsigned long long int));
			mass_cat = realloc(mass_cat, arrlen * sizeof(double));
	
			for (int i = 0; i < 3; i++)
			{
				pos[i] = realloc(pos[i], arrlen * sizeof(double));
			}

			if (!groupids)
			{
				fprintf(stderr, "Can't realloc\n");
				exit(1);
			}
		}
		
		for (igood=0; igood<ngood; igood++)
		{
			u_cat=fread(&dummy_cat,sizeof(int),1,fp_cat);
			u_cat=fread(&catdata,dummy_cat,1,fp_cat);
			u_cat=fread(&dummy_cat,sizeof(int),1,fp_cat);


			for (int i=0; i < 3; i++)
			{
				pos[i][nhalo] = catdata.x[i];
			}
			mass_cat[nhalo] = catdata.M;
			mass_cat[nhalo] = log10(mass_cat[nhalo]);

			groupids[nhalo] = catdata.id;

			nhalo++;
		}
	}
	fclose(fp_cat);
	printf("nhalo in cat: %d\n", nhalo);

	FILE *fp;

//	char dir[128] = "/mnt/data_cat5/rlzhang/Pinocchio/test/output/tests22/1.47/N1750";
	char dir[128] = "/mnt/data_cat4/moriwaki";
	sprintf(filename, "%s/pinocchio.%s.histories.out", dir, runname);
	printf("Opening file %s\n", filename);

	fp = fopen(filename, "r");

	if (fp==NULL)
	{
		printf("Error: history file %s not found\n", filename);
		return 1;
	}

	unsigned long long int ntrees_tot, nbranch_tot;
	int dummy, NSlices, u;

	u = fread(&dummy, sizeof(int), 1, fp);
	u = fread(&NSlices, sizeof(int), 1, fp);
	u = fread(&dummy, sizeof(int), 1, fp);

	printf("The run used %d slices\n", NSlices);

	ntrees_tot = nbranch_tot = 0;

	int ntrees, nbranch, thistree, nbranch_tree;

	int First = 1;


	FILE *fp_sfr;
	char fname_sfr[32];
	sprintf(fname_sfr, "mass-sfr.txt");
	fp_sfr = fopen(fname_sfr, "w");
	if (fp_sfr==NULL)
	{
		fprintf(stderr, "cannot open file: %s\n", fname_sfr);
		exit(1);
	}
	fprintf(fp_sfr, "# log10(mass[Msun/h]) log10(sfr[Msun/yr])\n");



	FILE *fp_sfr_sat;
	char fname_sfr_sat[32];
	sprintf(fname_sfr_sat, "mass-sfr-sat.txt");
	fp_sfr_sat = fopen(fname_sfr_sat, "w");

	if (fp_sfr_sat==NULL)
	{
		fprintf(stderr, "cannot open file: %s\n", fname_sfr_sat);
		exit(1);
	}

	fprintf(fp_sfr_sat, "# satellites\n");
	fprintf(fp_sfr_sat, "# log10(mass[Msun/h]) log10(sfr[Msun/yr])\n");

	
	double mhalo_now;

	int n_table = 999;
	double table_x[n_table];
	double table_y[n_table];
	double *cubic_set;
	cubic_set = c_pseudo_interp_table(n_table, table_x, table_y);

	double y2 = 0.35;
	double x2 = 11.20;
	double y1 = 1.60;
	double x1 = 10.0;
	double m = (y2-y1)/(x2-x1);
	double c = y1 - m*x1;

					
	FILE *fp_conc;
	char fname_conc[128];
	sprintf(fname_conc, "conc_diemer.txt");
	printf("Opening file %s\n", fname_conc);

	fp_conc = fopen(fname_conc, "r");
	if (fp_conc==NULL)
	{
		printf("Error: conc table file %s not found\n", fname_conc);
		return 1;
	}
	
	int nbins_m = 61;  
	int nbins_z = 57;

	double **conc_arr;
	conc_arr = (double**) malloc(nbins_z * sizeof(double));
	for (int ibin_z=0; ibin_z < nbins_z; ibin_z++)
	{
		conc_arr[ibin_z] = (double*) malloc(nbins_m * sizeof(double));
	}

	printf("nbins_m: %d\n", nbins_m);
	printf("nbins_z: %d\n", nbins_z);

	char str[1024];
	const char s[2] = " ";
	char *token;

	int ibin_z = 0;
	int ibin_m;
	while(fgets(str, 1024, fp_conc) != NULL)
	{	
		if (str[0] == '#') printf("%s\n", str);
		else
		{
//			printf("ibin_z: %d\n", ibin_z);
			
			char *rest = str;
			ibin_m = 0;
			while((token = strtok_r(rest, " ", &rest)))
			{
				//printf("ibin_m: %d\n", ibin_m);
				double value = strtod(token, &token);
				//printf("%lf\n", value);
			    conc_arr[ibin_z][ibin_m] = value;
				ibin_m ++;
			
			}
			if (ibin_m != nbins_m)
			{
				fprintf(stderr, "There are %d M bins in %s, but assigned %d bins\n",  ibin_m, fname_conc, nbins_m);
				exit(1);
			}
			ibin_z ++;
		}
	}
	printf("ibin_z: %d\n", ibin_z);

	if (ibin_z != nbins_z)
	{
		fprintf(stderr, "There are %d z bins in %s, but assigned %d bins\n",  ibin_z, fname_conc, nbins_z);
		exit(1);
	}
	fclose(fp_conc);


	double logMmin = 9.0;
	double dlogM = 0.1;
	printf("logMmin: %lf\n", logMmin);
	printf("dlogM: %lf\n", dlogM);

	double *mass_bins;
	mass_bins = (double*) malloc(nbins_m * sizeof(double));
	for (int ibin = 0; ibin < nbins_m; ibin++)
	{
		mass_bins[ibin] = logMmin + ibin*dlogM;
	}

	double zmin = 0.4;
	double dz = 0.1;
	printf("zmin for conc table: %lf\n", zmin);
	printf("dz: %lf\n", dz);

	double *z_bins;
	z_bins = (double*) malloc(nbins_z * sizeof(double));
	for (int ibin = 0; ibin < nbins_z; ibin++)
	{
		z_bins[ibin] = zmin + ibin*dz;
	}

	double log_eta_ha = 39.38;
	double e_ha = 0.951;  
	double chabrier_factor = 0.88;

	for (int ThisSlice = 0; ThisSlice < NSlices; ThisSlice++)
	{
		u = fread(&dummy, sizeof(int), 1, fp);
		u = fread(&ntrees, sizeof(int), 1, fp);
		u = fread(&nbranch, sizeof(int), 1, fp);
		u = fread(&dummy, sizeof(int), 1, fp);
		ntrees_tot += ntrees;
		nbranch_tot += nbranch;

		printf("Ntrees: %d\n", ntrees);
		printf("Nbranches: %d\n", nbranch);


		int imass;
		for (int itree = 0; itree < ntrees; itree++)
		{
			u = fread(&dummy, sizeof(int), 1, fp);
			u = fread(&thistree, sizeof(int), 1, fp);
			u = fread(&nbranch_tree, sizeof(int), 1, fp);
			u = fread(&dummy, sizeof(int), 1, fp);
			

			//printf("Tree: %d\n", thistree);
			//printf("Nbranches: %d\n", nbranch_tree);
		//
			double masses[nbranch_tree+1];
//			double mass_sums[nbranch_tree+1];
			double redshifts[nbranch_tree+1];
			int imw = 0;
			double M_out;
			int equals = 0;
			int index = -1;
			int found = 0;
			int mw_none = 1;
			int zmin_ok = 0;
			int mw_main[nbranch_tree];

			int *imw_sats;
			imw_sats = (int*) malloc(nbranch_tree* sizeof(int));
			for (int i = 0; i < nbranch_tree; i++)
			{
				imw_sats[i] = 0;
			}

			double **masses_sat;
			masses_sat = (double**) malloc((nbranch_tree) * sizeof(double));
			for (int ibranch=0; ibranch < nbranch_tree; ibranch++)
			{
				masses_sat[ibranch] = (double*) malloc((nbranch_tree) * sizeof(double));
			}

			double **redshifts_sat;
			redshifts_sat = (double**) malloc((nbranch_tree) * sizeof(double));
			for (int ibranch=0; ibranch < nbranch_tree; ibranch++)
			{
				redshifts_sat[ibranch] = (double*) malloc((nbranch_tree) * sizeof(double));
			}

			double M_host[nbranch_tree];

			start = clock();
			for (int ibranch = 0; ibranch < nbranch_tree; ibranch++)
			{
				u = fread(&dummy, sizeof(int), 1, fp);
				u = fread(&histdata, sizeof(hist), 1, fp);
				u = fread(&dummy, sizeof(int), 1, fp);

				double m1 = histdata.mass * pmass; //mass of halo at merger
				int mw = histdata.mw; // (larger) halo it merges with

				if (ibranch == 0)
				{
					M_out = m1;
					unsigned long long int id = histdata.name;
				}
				else if (mw == nbranch_tree) //merged with main halo
				{
					mw_none = 0;
					double m2 = histdata.mam *pmass; //mass of main halo it merges with at merger
					double zmerger = histdata.zme; // redshift at merger 
					double mass_sum = m1 + m2; //mass of resulting halo after merging
					//printf("mass1:%lf\n", log10(mass_sum));
					
					redshifts[imw] = zmerger;
					masses[imw] = m2;
					//mass_sums[imw] = mass_sum;
					/*
					imw++;
					redshifts[imw] = zmerger;
					masses[imw] = mass_sum;
					*/
					//masses[imw] = m2;

					
					// if found = 1 then already have redshift bin so don't need to execute if statements
					if (z1 == redshifts[imw] && found == 0)
					{
						// if z1 equals the merger redshift 
						index = imw;
						equals = 1;
						found = 1;
					}
					else if (z1 > redshifts[imw] && found == 0) 
					{
						// high redshift -> redshift
						// if z1 is higher than merger redshift then it falls in between this redshift and the one before 
						index = imw-1;
						found = 1;
						if (imw == 0) index = 0; //underestimates dynamical time
					}
					
					if (log10(M_out) > 10)
					{
						// before merging with main branch
						redshifts_sat[ibranch][imw_sats[ibranch]] = zmerger;
						masses_sat[ibranch][imw_sats[ibranch]] = m1;
						imw_sats[ibranch] ++;
						M_host[ibranch] = m2;
					}
									
					imw++;

				}
				else if (log10(M_out) > 10)
				{
					
					double zmerger = histdata.zme; // redshift at merger 
					//for halo it merges with
					double m2 = histdata.mam *pmass; //mass of main halo it merges with at merger
					double mass_sum = m1 + m2; //mass of resulting halo after merging
					
					redshifts_sat[mw][imw_sats[mw]] = zmerger;
					masses_sat[mw][imw_sats[mw]] = m2;
					imw_sats[mw]++;

					//for itself
					redshifts_sat[ibranch][imw_sats[ibranch]] = zmerger;
					masses_sat[ibranch][imw_sats[ibranch]] = m1;
					imw_sats[ibranch]++;

					M_host[ibranch] = m2;

				
					/*
					redshifts_sat[ibranch][0] = histdata.zap;
					masses_sat[ibranch][0] = min_hmass;
					*/
				}
				
				
			}
			end = clock();
     			cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
			//printf("after reading: %lf s\n", cpu_time_used);
			// high z -> low z: output redshift is at end of array
			redshifts[imw] = zout; 		
			masses[imw] = M_out;
			if (index == -1) index = imw-1; // default value, then in last bin (between second to last and last (output) redshift - z1 only higher than output redshift)
			double Md; // M(t - tdyn)
			double sfr = 0;
			if (mw_none == 0)
			{
				if (redshifts[0] > z1)
				{
				Md = M_dyn(masses, redshifts, index, z1, equals);

				/*
				Md = masses[imw-1];
				double t = my_z2t(redshifts[imw-1],cosm);
				tdyn = tout - t;
				*/

				//if (Md > 1E10)
				{
//				printf("Mout: %lf\n", log10(M_out));
				/*
				printf("Md: %lf\n", log10(Md));
				printf("Mpseudo: %lf\n", M_pseudo);
				*/

//				printf("Mpseudo_vir: %lf\n", log10(Mpseudo_vir));
				double dM_dt_dyn = (M_out - Md)/cosm.hubble/tdyn; 
				//printf("term1: %lf\n", dM_dt_dyn);
				//
				double logM_out = log10(M_out);
				double c179 = conc_ishiyama(logM_out, zout, mass_bins, z_bins, conc_arr);
				double Delta179 = 179;
				double x = omega_m(cosm, zout) - 1;
				double Delta_vir_crit = (18 * pow(M_PI,2)) + (82.0 * x) - (39.0 * pow(x,2)); // Bryan & Norman
				double Delta_vir_BN = Delta_vir_crit / omega_m(cosm,zout);
				double cBN = c_change_mass_def(cosm, c179, Delta179, Delta_vir_BN, zout, zout, cubic_set, n_table, table_x, table_y);

				double r_vir_out = rvir_from_mvir(cosm, M_out, zout);
				double rho_vir_out = rho_nfw(cosm, cBN, M_out, r_vir_out, r_vir_out); // [g/cm3] // depends on conc-mass relation
				double r_vir_d = rvir_from_mvir(cosm, Md, z1); //[cm]
//				r_vir_d = rvir_from_mvir(cosm, masses[imw-1], redshifts[imw-1]); //[cm]


				/*
				double dR=0;
				for (int i = index+1; i < imw; i++)
				{
					double z_a = redshifts[i];
					double r_vir_a = rvir_from_mvir(cosm, masses[i], redshifts[i]);
					double r_vir_b = rvir_from_mvir(cosm, masses[i+1], redshifts[i+1]);
			
					double logM_a = log10(masses[i]);
					c179 = conc_ishiyama(logM_a,redshifts[i] , mass_bins, z_bins, conc_arr);
					Delta179 = 179;
					x = omega_m(cosm, redshifts[i]) - 1;
					Delta_vir_crit = (18 * pow(M_PI,2)) + (82.0 * x) - (39.0 * pow(x,2)); // Bryan & Norman
					Delta_vir_BN = Delta_vir_crit / omega_m(cosm,z_a);
					cBN = c_change_mass_def(cosm, c179, Delta179, Delta_vir_BN, z_a, z_a, cubic_set, n_table, table_x, table_y);
					double rho_vir_a = rho_nfw(cosm, cBN, masses[i], r_vir_a, r_vir_a); // [g/cm3] // depends on conc-mass relation

					dR += (r_vir_b - r_vir_a)*r_vir_b*r_vir_b*rho_vir_a;
				}
				*/

				/*
				double logMd = log10(Md);
				c179 = conc_ishiyama(logMd, z1, mass_bins, z_bins, conc_arr);
				Delta179 = 179;
				x = omega_m(cosm, z1) - 1;
				Delta_vir_crit = (18 * pow(M_PI,2)) + (82.0 * x) - (39.0 * pow(x,2)); // Bryan & Norman
				Delta_vir_BN = Delta_vir_crit / omega_m(cosm,z1);
				cBN = c_change_mass_def(cosm, c179, Delta179, Delta_vir_BN, z1, z1, cubic_set, n_table, table_x, table_y);
				double rho_vir_d = rho_nfw(cosm, cBN, Md, r_vir_d, r_vir_d); // [g/cm3] // depends on conc-mass relation
				printf("rho_vir_out: %g\n", rho_vir_out);
				printf("rho_vir_d: %g\n", rho_vir_d);
				*/
				/*
				printf("mass_sum: %lf\n", log10(mass_sums[imw-1]));
				printf("m2: %lf\n", log10(masses[imw-1]));
				*/

				/*
				printf("zmerger: %lf\n", redshifts[imw-1]);
				printf("mass_sum: %lf\n", log10(masses[imw-1]));
				printf("mass_out: %lf\n", log10(masses[imw]));
				printf("r_merge: %lf\n", r_vir_d/kpc);
				printf("r_out: %lf\n", log10(r_vir_out));
				*/
				
				double dR_dt_dyn = (r_vir_out - r_vir_d) / tdyn; // "virial radius" calculated from fof mass
				double term2 = 4 * M_PI * r_vir_out*r_vir_out * rho_vir_out * dR_dt_dyn / Msun;  
//				double term2 = 4 * M_PI * r_vir_d*r_vir_d * rho_vir_d * dR_dt_dyn / Msun;  

				double dM_dt = dM_dt_dyn - term2;
				/*
				if ((M_out - mass_sums[imw-1])/tdyn < term2)
				{
					dM_dt = (mass_sums[imw-1] - masses[imw-1])/tdyn;
				}
				*/
//				if (dM_dt < 0) dM_dt = dM_dt_dyn - term2*0.8;
//				printf("term2: %lf\n", term2);
//				printf("dM_dt: %lf\n", dM_dt);
				double f_b = cosm.omegab / cosm.omega;
				double dmb_dt = f_b * dM_dt;
				double e = baryon_conversion_efficiency(M_out/cosm.hubble, zout);
				sfr = dmb_dt * e;
//				double log_lum_cent = (log10(sfr*chabrier_factor) + log_eta_ha)/e_ha;
//				log_lum_cent = log10(sfr/(4.4 *pow(10,-42)));
				//printf("%lf\n", log_lum_cent);
				if (log10(M_out) > 9 && sfr > 0)
				{
//				fprintf(fp_sfr, "%lf %lf\n", log10(M_out), log_lum_cent);
//				printf("pos
				fprintf(fp_sfr, "%lf %lf %lf %lf %lf %lf %lf\n", log10(M_out), pos[0][itree], pos[1][itree], pos[2][itree], cBN, r_vir_out, log10(sfr));
			//	fprintf(fp_sfr, "%lf %lf\n", log10(M_out), dM_dt);
				}
				}
				}
			}
			
			else
			{
				if (log10(M_out) > 9)
				{
			//	fprintf(fp_sfr, "%lf -inf\n", log10(M_out));
				}
			}
			
			
			start=clock();
			//printf("M_out: %lf\n", log10(M_out));
			
			if (log10(M_out) > 10)
			{
				fprintf(fp_sfr_sat, "# %lf\n", log10(M_out));
			for (int ibranch = 0; ibranch < nbranch_tree; ibranch++)
			{
				int nmw_sat = imw_sats[ibranch];
				double *redshifts_sat1 = redshifts_sat[ibranch];
				double *masses_sat1 = masses_sat[ibranch];

				
				//if ((nmw_sat > 1) && (log10(masses_sat1[nmw_sat-1]) > 9))
				if (nmw_sat > 1) 
				{
					double z_infall = redshifts_sat1[nmw_sat-1];
					double M_infall = masses_sat1[nmw_sat-1];
					//printf("%lf\n", redshifts_sat1[nmw_sat-1]);
					double tdyn = t_dyn(cosm, z_infall);
					double Om = omega_m(cosm, z_infall);
					double Ol = 1 - Om;
					double gz = (5.0/2.0) * Om / (pow(Om, 4.0/7.0) - Ol + (1 + Om/2)*(1 + Ol/70)); 
					double Dz = gz / (1+z_infall);
					double Om_0 = omega_m(cosm, 0);
					double gz_0 = (5.0/2.0) * cosm.omega / (pow(cosm.omega, 4.0/7.0) - cosm.omega + (1 + Om_0/2)*(1 + cosm.lambda/70)); 
					double Dz_0 = gz_0 / (1);
					Dz = Dz/Dz_0;
//					printf("Dz: %lf\n", Dz);
					double A = 0.195 * pow(Dz/0.6,2);
					double b = 0.92 * Dz;
					double c = 1.9;
					double eta = 1;
					double Mhost = M_host[ibranch];
//					printf("Mhost: %lf\n", log10(Mhost));
					double tau_merger = A * tdyn * pow(Mhost/M_infall,b) * exp(c * eta) / log(1 + Mhost/M_infall);



//					printf("%lf\n", M_infall/Mhost);
//					printf("tau_merger:%g\n", tau_merger/1e9);
					double t_merge = my_z2t(z_infall, cosm);
					double t1_sat = t_merge - tdyn;
					double z1_sat = my_t2z(t1_sat, cosm);
//					printf("z1=%lf\n", z1);
					//printf("tdyn=%lf\n", tdyn/1E9);

					double t = tout - t_merge;
					double t_Gyr = t/1E9;
					//printf("t=%lf Gyr\n", t/1E9);
					//printf("M=%lf\n", log10(masses_sat1[nmw_sat-1]));
					//printf("mhalo: %lf\n", log10(masses_sat1[nmw_sat-1]));
					double t3_Gyr = 2 * pow(1+z_infall, -3.0/2.0);

//					printf("t3: %lf\n", t3_Gyr);
					double t3 = t3_Gyr*1E9;
					if (z1_sat < 6.0 && z1_sat < redshifts_sat1[0])// && t < tau_merger)
					{
//					printf("t=%lf Gyr\n", t/1E9);
//						printf("t3: %lf\n", t3_Gyr);
						int index1;
						for (int iz=0; iz < nmw_sat; iz++)
						{
							if (z1_sat < redshifts_sat1[iz])
							{
								index1 = iz - 1;
								break;
							}
							else if (z1_sat == redshifts_sat1[iz])
							{
								index1 = iz;
								equals = 1;
								break;
							}
						}

						Md = M_dyn(masses_sat1, redshifts_sat1, index1, z1_sat, equals);
						//double Mpseudo = M_pseudo(cosm, Md, z1_sat, z_infall, cubic_set, n_table, table_x, table_y, mass_bins, z_bins, conc_arr);
						double dM_dt_dyn = (M_infall - Md)/cosm.hubble/tdyn; 


						//printf("term1: %lf\n", dM_dt_dyn);
						//
						double logM_infall = log10(M_infall);
						double c179 = conc_ishiyama(logM_infall, z_infall, mass_bins, z_bins, conc_arr);
						double Delta179 = 179;
						double x = omega_m(cosm, z_infall) - 1;
						double Delta_vir_crit = (18 * pow(M_PI,2)) + (82.0 * x) - (39.0 * pow(x,2)); // Bryan & Norman
						double Delta_vir_BN = Delta_vir_crit / omega_m(cosm,z_infall);
						double cBN = c_change_mass_def(cosm, c179, Delta179, Delta_vir_BN, z_infall, z_infall, cubic_set, n_table, table_x, table_y);

						double r_vir_infall = rvir_from_mvir(cosm, M_infall, z_infall);
						double rho_vir_infall = rho_nfw(cosm, cBN, M_infall, r_vir_infall, r_vir_infall); // [g/cm3] // depends on conc-mass relation
						double r_vir_d = rvir_from_mvir(cosm, Md, z1); //[cm]
						
						double dR_dt_dyn = (r_vir_infall - r_vir_d) / tdyn; // "virial radius" calculated from fof mass
						double term2 = 4 * M_PI * r_vir_infall*r_vir_infall * rho_vir_infall * dR_dt_dyn / Msun;  
						double dM_dt = dM_dt_dyn - term2; // [Msun/yr]


						//printf("dM_dt: %lf\n", dM_dt);
						double f_b = cosm.omegab / cosm.omega;
						double dmb_dt = f_b * dM_dt;
						double e = baryon_conversion_efficiency(M_infall/cosm.hubble, z_infall);
						double sfr_sat = dmb_dt * e;
												
						//printf("sfr_sat: %g\n", sfr_sat);
						
						double Mstar = sfr2mstar(cosm, sfr_sat, z_infall); //[Msun]
						//printf("Mstar: %g\n", Mstar);

						Mstar = mhalo2mstar(cosm, M_infall, z_infall);
						double frac = sfr_sat;
//						printf("Mstar: %g\n", Mstar);
						//printf("frac: %g\n", frac);
						//t3 /= frac;
						//printf("t3: %g\n", t3/1e9);


						//printf("%lf\n", f_b);
						double ssfr = sfr_sat / Mstar;
						double f_s = 0.03;
						double f_strip = 0;
						double eta1 = 2.5;
						double R = 0.4;

						double T_delay = (f_b - f_s*((0.1 * pow(1+z_infall,2)) + (eta1 * (1+R)) + 1) - f_strip)/(f_s * (1 - R + eta1) * ssfr);
//						printf("ssfr: %g\n", ssfr);
						double T_delay_Gyr = T_delay/1e9;
						//printf("z_infall: %lf\n", z_infall);
						//printf("T_delay: %g\n", T_delay_Gyr);
						t3 = T_delay;
						t3_Gyr = T_delay_Gyr;
						//double t_fade = 4 * pow(1 + z_infall, -3/2);
						//printf("T_delay: %g\n", T_delay_Gyr);
						double t_fade = 0.5;

						Mstar = pow(10,10.5);
						double a = 0.7 - 0.13*z_infall;
						double b = 0.38 + 1.14*z_infall - 0.19*z_infall*z_infall;
						double logsfr = a * (log10(Mstar) - 10.5) + b;
						double z0 = 0.05;
						double a0 = 0.7 - 0.13*z0;
						double b0 = 0.38 + 1.14*z0 - 0.19*z0*z0;
	
						double logsfr0 = a0 * (log10(Mstar) - 10.5) + b0;
						double t_quench = 4.4 * pow(10, logsfr0) / pow(10, logsfr);
						t_quench = 1.7*pow(1+zout, -3.0/2.0);

						/*
						printf("sfr1: %lf\n", pow(10,logsfr));
						printf("sfr_sat: %lf\n", sfr_sat);
						printf("t: %g\n", 4.4*pow(10,logsfr0)/pow(10,logsfr));
						*/

						
						//if (t > t3)
						{
						//sfr_sat *= exp(-(t_Gyr-t3_Gyr)/0.25);
						sfr_sat *= exp(-(t_Gyr)/t_quench);
						}
						
						if (log10(sfr_sat) > -1)
						{
						//fprintf(fp_sfr_sat, "%lf %g\n", log10(M_out), log10(sfr_sat));
						//printf("%lf %g\n", log10(M_out), log10(sfr_sat));
						}
						if (sfr_sat > 0)
						{
						//double log_lum_sat = (log10(sfr_sat*chabrier_factor) + log_eta_ha)/e_ha;
						//log_lum_sat = log10(sfr_sat/(4.4 *pow(10,-42)));
//						fprintf(fp_sfr_sat, "%g\n", log_lum_sat);
						}

						//printf("t: %lf\n", t_Gyr);
						//printf("t3: %lf\n", t3_Gyr);
					//	printf("sfr: %lf\n", log10(sfr_sat));
						if (sfr_sat> 0) sfr += sfr_sat;
					}

				}
			}
			fprintf(fp_sfr_sat, "%lf %g\n", log10(M_out), log10(sfr));
			end = clock();
     			cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
			//printf("assigning sats: %lf s\n", cpu_time_used);
			}
		
			

		}
	}
	printf("I found %Ld trees and %Ld branches in the file\n", ntrees_tot, nbranch_tot);
		

	fclose(fp);

	fclose(fp_sfr);
	fclose(fp_sfr_sat);

	time(&now);
	printf("%s\n", ctime(&now));
}

double E_z(Cosmology cosm, double z)
{

	double zp1 = 1.0 + z;

	double Ez = sqrt(cosm.omega*pow(zp1,3) + pow((1.0 - cosm.omega - cosm.lambda),2) + cosm.lambda); 
	return Ez;
}

double omega_m(Cosmology cosm, double z)
{
	double zp1 = 1.0 + z;

	double Ez = sqrt(cosm.omega*pow(zp1,3) + pow((1.0 - cosm.omega - cosm.lambda),2) + cosm.lambda); 
	double omega_m = cosm.omega * pow(zp1,3)/pow(Ez,2);

	return omega_m;
}

double rho_crit(Cosmology cosm, double z)
{
	// H0 [/s]
	// G [cm3 g-1 s-2]
	double Ez = E_z(cosm, z);

	double rho_crit = 3 * pow(Ez*cosm.hubble*H0,2)/(8*M_PI*G); 
	return rho_crit; // [g/cm3]
}

double rho_m(Cosmology cosm, double z)
{
	double rho_m = omega_m(cosm,z)*rho_crit(cosm,z);
	return rho_m; // [g/cm3]
}

double t_dyn(Cosmology cosm, double z)
{
	// H0               3.226e-18       // Hubble constant(/s) //
	double zp1 = 1.0 + z;
	
	double Ez = sqrt(cosm.omega*pow(zp1,3) + pow((1.0 - cosm.omega - cosm.lambda),2) + cosm.lambda); 
	double omega_m = cosm.omega * pow(zp1,3)/pow(Ez,2);
	double x = omega_m - 1; 

	double Delta_vir_crit = (18 * pow(M_PI,2)) + (82.0 * x) - (39.0 * pow(x,2)); // Bryan & Norman
	double Delta_vir = Delta_vir_crit/omega_m;

//	double Delta_vir = 35.4;
//	double rho_crit = 3 * pow(Ez*cosm.hubble,2)/(8*M_PI*G); 

	// tdyn = pow(3/(4 * G * M_PI * Delta_vir * rho_crit),0.5);
	double H  = Ez * (H0*cosm.hubble);
	double tdyn_s = pow(2/(Delta_vir*omega_m),0.5)/ H;//dynamical time in seconds
	double tdyn_yr = tdyn_s / yr;//dynamical time in years
	return tdyn_yr;
}




double rvir_from_mvir(Cosmology cosm, double mvir_sol_h, double z)
{
	// Need to consider units
//	double omega_m = omega_m(cosm,z);
	double x = omega_m(cosm,z) - 1; 
	double Delta_vir_crit = (18 * pow(M_PI,2)) + (82.0 * x) - (39.0 * pow(x,2));

	double Delta_vir = Delta_vir_crit/omega_m(cosm,z);
//	double Delta_vir = 35.4;
	
	double mvir = mvir_sol_h * Msun /cosm.hubble; //[g]
	double rho_mean = rho_m(cosm, z);

	double rvir = pow((3 * mvir)/(4*M_PI*rho_mean*Delta_vir), 1.0/3.0); // [cm]
	return rvir; // [cm]
}



double my_t2z(double t_yr, Cosmology cosm)
{
	double tnow_ref = 13.7957 * pow(10,9); //from cosmology calculator for Planck 2018 cosmology
	double tnow_cosm = ztotime((double)0, cosm)/ (H0*cosm.hubble*yr); // note this is not equal to 1 for z = 0
	double tnow_frac = tnow_cosm/tnow_ref; // to account for inaccuracy of code (this is only approx correction)
	
	double t_s = t_yr * yr;
	double t_cosm = t_s * (H0 * cosm.hubble);
	double z = timetoz(t_cosm,cosm);
	return z;
}

double my_z2t(double z, Cosmology cosm)
{
	double tnow_ref = 13.7957 * pow(10,9); //from cosmology calculator for Planck 2018 cosmology
	double tnow_cosm = ztotime((double)0, cosm) / (H0*cosm.hubble*yr);
	double tnow_frac = tnow_cosm/tnow_ref;

	double t_cosm = ztotime(z, cosm);
	double t_s = t_cosm / (H0 * cosm.hubble);

	double t_yr = t_s / yr;
	return t_yr;
}

double rho_nfw(Cosmology cosm, double conc, double m_vir_sol_h, double r_vir, double r)
{
	double m_vir = m_vir_sol_h * Msun / cosm.hubble; // [g]
	double rho_halo = m_vir / (4 * M_PI * pow(r_vir,3)/3);
	double A_NFW = log(1 + conc) - conc/(1+conc);

	double x = r/r_vir;

	double rho = rho_halo / (3 * A_NFW * x * pow(1/conc + x, 2));

	return rho; // [g/cm3]
}

double baryon_conversion_efficiency(double M, double z)
{
	//2018
	double M0 = 11.339;
	double Mz = 0.692;
	double e0 = 0.005;
	double ez = 0.689;
	double beta0 = 3.344;
	double betaz = -2.079;
	double gamma0 = 0.966;

	//2020
	/*
	double M0 = 11.34829;
	double Mz = 0.654238;
	double e0 = 0.009010;
	double ez = 0.596666;
	double beta0 = 3.094621;
	double betaz = -2.019841;
	double gamma0 = 1.107304;
	*/

	double a = 1 / (1 + z);

	
	double log10Mc = M0 + Mz*(1-a); // characteristic mass
	double Mc = pow(10, log10Mc);
	double eN = e0 + ez*(1-a); // normalisation
	double beta = beta0 + betaz*(1-a); 
	double gamma = gamma0;

	
	double e; // instantaneous baryon conversion efficiency
	e = 2 * eN / (pow(M/Mc,-beta) + pow(M/Mc,gamma));


	return e;
}

double modelKlypin16(double M, double z, char *mdef)
{
	
	/*
	if (M < 1E10)
	{
		fprintf(stderr, "Invalid mass for Klypin et al 2016 m-based model, %g\n", M);
		exit(0);
	}
	*/
	


	double z_bins[10] = {0.0, 0.35, 0.5, 1.0, 1.44, 2.15, 2.5, 2.9, 4.1, 5.4};
	double C0_bins[10];
	double gamma_bins[10];
	double M0_bins[10];
	if (!strcmp(mdef,"200c"))
	{
		static const double C0_bins1[10] = {7.4, 6.25, 5.65, 4.3, 3.53, 2.7, 2.42, 2.2, 1.92, 1.65};
		static const double gamma_bins1[10] = {0.120, 0.117, 0.115, 0.110, 0.095, 0.085, 0.08, 0.08, 0.08, 0.08};
		static const double M0_bins1[10] = {5.5E5, 1E5, 2E4, 900.0, 300.0, 42.0, 17.0, 8.5, 2.0, 0.3};
 		memcpy(C0_bins, C0_bins1, sizeof(C0_bins1));
 		memcpy(gamma_bins, gamma_bins1, sizeof(gamma_bins1));
 		memcpy(M0_bins, M0_bins1, sizeof(M0_bins1));
	}
	else if (!strcmp(mdef, "vir"))
	{
		printf("2\n");
		static const double C0_bins1[10] = {9.75, 7.25, 6.5, 4.75, 3.8, 3.0, 2.65, 2.42, 2.1, 1.86};
		static const double gamma_bins1[10] = {0.110, 0.107, 0.105, 0.1, 0.095, 0.085, 0.08, 0.08, 0.08, 0.08};
		static const double M0_bins1[10] = {5E5, 2.2E4, 1E4, 1000.0, 210.0, 43.0, 18.0, 9.0, 1.9, 0.42};
	
		memcpy(C0_bins, C0_bins1, sizeof(C0_bins1));
 		memcpy(gamma_bins, gamma_bins1, sizeof(gamma_bins1));
 		memcpy(M0_bins, M0_bins1, sizeof(M0_bins1));

	}
	else
	{
		fprintf(stderr, "Invalid mass definition for Klypin et al 2016 m-based model, %s\n", mdef);
		exit(0);
	}

	int bin;
	for (int ibin = 0; ibin < 10; ibin++)
	{
		if (z < z_bins[ibin])
		{ 
			bin = ibin-1; 
			break;
		}
		else if (z > 5.4) 
		{
			printf("z = %lf > 5.4\n", z);
			exit(0);
		}
	}

	// In my other code I interpolated the concentrations. Here I'm interpolating the parameters
	double C0 = lin_interp(z_bins, C0_bins, z, bin);
	double gamma = lin_interp(z_bins, gamma_bins, z, bin);
	double M0 = lin_interp(z_bins, M0_bins, z, bin);
	M0 *= 1E12;

	double c = C0 * pow(M/1E12, -gamma) * (1.0 + pow(M/M0, 0.4));

	return c;

}

double lin_interp(double x_arr[], double y_arr[], double z, int ibin)
{
	double m = (y_arr[ibin+1] - y_arr[ibin]) / (x_arr[ibin+1] - x_arr[ibin]);
	double c = y_arr[ibin] - m*x_arr[ibin];

	double interp_val = m*z +c;
	return interp_val;
}

double mu(double conc)
{
	double mu = log(1+conc) - conc/(1.0+conc);
	return mu;
}


double *c_pseudo_interp_table(int n_table, double table_x[], double table_y[])
{
	
	for (int i = 0; i < n_table; i++)
	{
		double power = 4.0 - (8.0/n_table)*i; 	
		table_x[i] = pow(10,power);
		table_y[i] = mu(table_x[i])/pow(table_x[i],3);
	}
	double *cubic_set;
	cubic_set = spline_cubic_set(n_table, table_y, table_x, 3, 0, 3, 0);
	return cubic_set;
}

double c_pseudo(Cosmology cosm, double ci, double zi, double zf, double *cubic_set, int n_table, double table_x[], double table_y[])
{
	//xDelta in COLOSSUS

	double rho_i = rho_m(cosm, zi); // [g/cm3]
	double rho_f = rho_m(cosm, zf); // [g/cm3]
	double y;
	y = (mu(ci)/pow(ci,3))*(rho_f/rho_i);
	double ypval;
	double yppval;
	double c2 = spline_cubic_val(n_table, table_y, table_x, cubic_set, y, &ypval, &yppval);
	return c2;
}

double M_pseudo(Cosmology cosm, double Mi_pinocchio, double zi, double zf, double *cubic_set, int n_table, double table_x[], double table_y[], double mass_bins[], double z_bins[], double **conc_arr)
{

	double logMi_pinocchio = log10(Mi_pinocchio);
	double c35 = conc_ishiyama(logMi_pinocchio, zi, mass_bins, z_bins, conc_arr);
//	printf("%lf\n", ci);
	double Delta35 = 35;
	double f_a = 0.25;
	double Delta_pinocchio = 3/(4*M_PI*pow(f_a,3));

	double x = omega_m(cosm,zf) - 1; 
	double Delta_vir_crit = (18 * pow(M_PI,2)) + (82.0 * x) - (39.0 * pow(x,2)); // Bryan & Norman
	double Delta_vir_BN = Delta_vir_crit / omega_m(cosm,zf);
//	
	double ci_pinocchio = c_change_mass_def(cosm, c35, Delta35, Delta_pinocchio, zi, zi, cubic_set, n_table, table_x, table_y);// initial concentration for whole halo

	double cf_pinocchio = c_pseudo(cosm, ci_pinocchio, zi, zf, cubic_set, n_table, table_x, table_y); //final conc for whole halo
	double Mf = Mi_pinocchio * mu(cf_pinocchio) / mu(ci_pinocchio); // eq. 7 of Diemer 2013

	return Mf;
}

double c_change_mass_def(Cosmology cosm, double ci, double Delta_i, double Delta_f, double zi, double zf, double *cubic_set, int n_table, double table_x[], double table_y[])
{
	//xDelta in COLOSSUS
	//Change from delta_vir_2lpt (35.4) to delta_pinocchio (15.28)
	//
	double rho_i = rho_m(cosm, zi);
	double rho_f = rho_m(cosm, zf);
	double time3 = clock();
	double y;
	y = (mu(ci)/pow(ci,3))*(Delta_f*rho_f/(Delta_i*rho_i));
	double ypval;
	double yppval;
	double c2 = spline_cubic_val(n_table, table_y, table_x, cubic_set, y, &ypval, &yppval);
	double end = clock();
	return c2;
}

double Delta_fof(double conc)
{
	double n_c = 0.652960;
	double b = 0.2;
	double Delta = (3 * n_c * pow(b,-3) * mu(conc) * pow(1+conc,2))/pow(conc,2) - 1;
	return Delta;
}

double M_pinocchio2vir(Cosmology cosm, double M_pinocchio, double z, double *cubic_set, int n_table, double table_x[], double table_y[], double mass_bins[], double z_bins[], double **conc_arr)
{
	double M_vir;
	double Delta35 = 35;
	double x = omega_m(cosm,z) - 1; 
	double Delta_vir_crit = (18 * pow(M_PI,2)) + (82.0 * x) - (39.0 * pow(x,2)); // Bryan & Norman
	double Delta_vir_BN = Delta_vir_crit / omega_m(cosm,z);
//	printf("Delta_vir_BN: %lf\n", Delta_vir_BN);
	double logM_pinocchio = log10(M_pinocchio);
	
	double c35 = conc_ishiyama(logM_pinocchio, z, mass_bins, z_bins, conc_arr); // conc for this Delta
	double b = 0.2; //linking length
	double Delta_pinocchio = 9 / (2 * M_PI * pow(b,3));
	//printf("Delta_pinocchio: %lf\n", Delta_pinocchio);
	double c_pinocchio = c_change_mass_def(cosm, c35, Delta35, Delta_pinocchio, z, z, cubic_set, n_table, table_x, table_y);// initial concentration for whole halo
	double c_vir = c_change_mass_def(cosm, c_pinocchio, Delta_pinocchio, Delta_vir_BN, z, z, cubic_set, n_table, table_x, table_y);// initial concentration for whole halo


	
	/*
	double r_pinocchio = pow(M_pinocchio * 3 / (4 * M_PI * Delta_pinocchio* rho_m(cosm,z)), 1.0/3.0);
//	printf("r_pinocchio: %lf\n", log10(r_pinocchio));
	double r_s = r_pinocchio / c_pinocchio;
	double rho_s = M_pinocchio / (4 * M_PI * pow(r_s,3) * mu(c_pinocchio)); // eq. 3 of Diemer+2013

	double r_vir = c_vir * r_s;
//	printf("r_vir: %lf\n", log10(r_vir));
	double rho_vir = Delta_vir_BN * rho_m(cosm, z);
	double M_vir = 4 * M_PI * pow(r_vir,3) * rho_vir/3; 
	printf("M_vir: %lf\n", log10(M_vir));
	*/

	//printf("M_pinocchio: %lf\n", log10(M_pinocchio));
	M_vir = M_pinocchio * mu(c_vir) / mu(c_pinocchio); // this is for static density profile i.e. dR = 0 
	//printf("M_vir: %lf\n", log10(M_vir));
//	else M_vir = M_pinocchio;
	return M_vir;
}

double M_dyn(double masses[], double redshifts[], int index, double z1, int equals)
{
	// M(t - tdyn)

	double Md;
	if (equals == 1)
	{
		// if exactly equal to a merger redshift
		double logMd = masses[index];
		Md = pow(10,logMd);
	}
	else if (redshifts[0] < z1)
	{
		Md = masses[0]; // this will underestimate sfr
		//z1 = redshifts[0];
	}
	else
	{
		double z_a = redshifts[index];
		double z_b = redshifts[index+1];

		double m_a = masses[index];
		double m_b = masses[index+1];

		double m = (m_b - m_a)/(log10(1+z_b) - log10(1+z_a));
		double c = m_a - log10(1+z_a)*m;

		double y = m*log10(1+z1) + c;
		Md = y;
	}
	return Md;
}

double tau(Cosmology cosm, double M_star_msun_h, double z)
{
	double tdyn = t_dyn(cosm, z);
	double tau_0 = 4.282;
	double tau_s = 0.363;

	double tau = tdyn * tau_0 * pow((M_star_msun_h/cosm.hubble)/1E10, -tau_s);

	return tau; 
}

double mhalo2mstar(Cosmology cosm, double M_halo_msun_h, double z)
{

	double M_halo = M_halo_msun_h/cosm.hubble; //[Msun]

	double M10 = 11.590;
	double M11 = 1.195;
	double N10 = 0.0351;
	double N11 = -0.0247;
	double beta10 = 1.376;
	double beta11 = -0.826;
	double gamma10 = 0.608;
	double gamma11 = 0.329;

	double a = 1/(z+1);

	double logM1 = M10 + M11*(1-a); 
	double N = N10 + N11*(1-a);
	double beta = beta10 + beta11*(1-a);
	double gamma = gamma10 + gamma11*(1-a);

	double M1 = pow(10,logM1);

	double m_star = M_halo_msun_h * 2 * N / (pow(M_halo/M1,-beta) + pow(M_halo/M1, gamma));
	return m_star;
}

double sfr2mstar(Cosmology cosm, double sfr, double z)
{
	// Tomczak 2016 eq.2
	double s0 = 0.195 + 1.157*z - 0.143*z*z;
	double logM0 = 9.244 + 0.753*z - 0.090*z*z;
	double gamma = 1.118;

	//printf("sfr: %lf\n", sfr);
	double Mstar = pow(pow(10,s0 - log10(sfr)) - 1, -1/gamma)*pow(10,logM0);
	//Mstar *= cosm.hubble;

	return Mstar; // [Msun]


}

double conc_ishiyama(double logM, double z, double mass_bins[], double z_bins[], double **conc_arr)
{

	double logMmin = mass_bins[0];
	double zmin = z_bins[0];
	double dlogM = mass_bins[1] - mass_bins[0];
	double dz = z_bins[1] - z_bins[0];

	int ibin_z = (int) (z - zmin)/dz;
	int ibin_m = (int) (logM - logMmin)/dlogM;

	double conc_z0 = lin_interp(mass_bins, conc_arr[ibin_z], logM, ibin_m);
	double conc_z1 = lin_interp(mass_bins, conc_arr[ibin_z+1], logM, ibin_m);

	double m = (conc_z1 - conc_z0) / (z_bins[ibin_z+1] - z_bins[ibin_z]);
	double c = conc_z0 - m*z_bins[ibin_z];
	double conc = m*z + c;

	return conc;
}
