#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <sys/stat.h>
#include <time.h>

#include "constant.h"
#include "para.h"
#include "proto.h"

#define NBIN_L 100
#define NLINES 1

int main(int argc, char *argv[])
{
	time_t now;
	time(&now);
	printf("%s\n", ctime(&now));


	if (argc != 3)
	{
		fprintf(stderr, "ERROR: incorrect number of arguments to %s\n", argv[0]);
		fprintf(stderr, "	Input: ");
		for(int i = 0; i < argc; i++) fprintf(stderr, "%s ", argv[i]);
		fprintf(stderr, "\n");
		fprintf(stderr, "	Usage: %s lambda_obs boxsize\n", argv[0]);

		exit(1);
	}

	char *linenames[NLINES] = {"Ha", "OIII", "OII"};
	double lambda_rest_Ha = 6563 * Angstrom;
	double lambda_rest_O3 = 5007 * Angstrom;
	double lambda_rest_O2 = 3727 * Angstrom;
	double lambda_rest[NLINES] = {lambda_rest_Ha, lambda_rest_O3, lambda_rest_O2};


	double lambda_obs = atof(argv[1]) * 10000 * Angstrom;
	double nu_obs = cspeed/lambda_obs;

	double redshifts[NLINES];
	for (int iline = 0; iline < NLINES; iline++)
	{
		redshifts[iline] = lambda_obs/lambda_rest[iline] - 1;
	}

	double extinction[NLINES] = {1.0, 1.35, 1.0};

	double R = 41; //spectral resolution


	double xylen = atof(argv[2]); //Mpc


	double L_max = 3;
	double L_min = -4; //42.1;
	printf("max logL = %lf\n", L_max);
	printf("min logL = %lf\n", L_min);
	double L_mins[NLINES] = {41.8, 42.1,41.9};
	double L_bin[NBIN_L] = {};
	double dL = (L_max - L_min) / NBIN_L;
	printf("%lf\n", dL);
	double dLs[NLINES];
	for (int iline = 0; iline < NLINES; iline++)
	{
		dLs[iline] = (L_max - L_mins[iline])/NBIN_L;
	}
	for (int ibin = 0; ibin < NBIN_L; ibin++) 
	{
		//L_bin[ibin] = L_min + (ibin+0.5)*dL;
		L_bin[ibin] = 0.0;
		//printf("%lf\n", L_bin[ibin]);
	}
	int L[NLINES][NBIN_L] = {};


	for (int iline = 0; iline < NLINES; iline++)
	{
		double z;
		z = redshifts[iline];
		z = round(z*10)/10; //round to 1dp

		double deltaz = (1 + z) / spectral_resolution;
		double zstart = z - deltaz/2;

		FILE *fp;
		char fname[128];
		sprintf(fname, "mass-sfr.txt");
		fp = fopen(fname, "r");
		if (fp==NULL)
		{
			fprintf(stderr, "cannot open file: %s\n", fname);
			exit(1);
		}
	
		FILE *fp_sat;
		char fname_sat[128];
		sprintf(fname_sat, "mass-sfr-sat.txt");
		fp_sat = fopen(fname_sat, "r");
		if (fp_sat==NULL)
		{
			fprintf(stderr, "cannot open file: %s\n", fname_sat);
			exit(1);
		}

		int arrlen = 100000;
		const int STEPSIZE = 100000;


		double** L_lines_all;
		L_lines_all = (double**)malloc(NLINES * sizeof(double*));
		for (int line = 0; line < NLINES; line++)
		{
			L_lines_all[line] = (double*) malloc(arrlen*sizeof(double));
		}


		double** L_lines;
		L_lines = (double**)malloc(NLINES * sizeof(double*));
		for (int line = 0; line < NLINES; line++)
		{
			L_lines[line] = (double*) malloc(arrlen*sizeof(double));
		}

		char str[256];
		int i = 0;
		int j = 0;
		
		int k = 0;

		int count_cent = 0;
		while(fgets(str, 256, fp) != NULL)
		{
			if (i == arrlen)
			{
				arrlen += STEPSIZE;
				for (int line = 0; line < NLINES; line++)
				{
					L_lines_all[line] = realloc(L_lines_all[line], arrlen * sizeof(double));
				}

				for (int line = 0; line < NLINES; line++)
				{
					L_lines[line] = realloc(L_lines[line], arrlen * sizeof(double));
				}

				if (!L_lines_all[0])
				{
					fprintf(stderr, "Can't realloc\n");
					exit(1);
				}
			}

			if (str[0] == '#') continue;
			else
			{
				sscanf(str, "%*lf %lf ", &L_lines_all[0][i] );
				if (isnan(L_lines_all[iline][i]) == 0 && isinf(L_lines_all[iline][i]) == 0)
				{
				L_lines[iline][j] = L_lines_all[iline][i];
				j++;
			
				i++;
				count_cent++;
				}
			}
			k++;
		}
		fclose(fp);
		printf("%d\n", k);
	
		printf("number of central subhaloes in box: %d\n", count_cent);	
		time(&now);
		printf("%s\n", ctime(&now));

		int count_sat = 0;
		int big_count = 0;
		int small_count = 0;
		
		
		while(fgets(str, 256, fp_sat) != NULL)
		{
			if (i == arrlen)
			{
				arrlen += STEPSIZE;
				for (int line = 0; line < NLINES; line++)
				{
					L_lines_all[line] = realloc(L_lines_all[line], arrlen * sizeof(double));
				}

				for (int line = 0; line < NLINES; line++)
				{
					L_lines[line] = realloc(L_lines[line], arrlen * sizeof(double));
				}

				if (!L_lines_all[0])
				{
					fprintf(stderr, "Can't realloc\n");
					exit(1);
				}
			}

			if (str[0] == '#') continue;
			else
			{
				sscanf(str, "%*s %lf ", &L_lines_all[0][i] );
				L_lines[iline][j] = L_lines_all[iline][i];
				j++;
			
				i++;
				count_sat++;
			}
			k++;
		}
		
		
		fclose(fp_sat);
		printf("%d\n", k);
	
		printf("number of satellite subhaloes in box: %d\n", count_sat);	
		int nhalo = j;
		printf("total number of subhalos: %d\n", nhalo);


		for (int j = 0; j < nhalo; j++)
		{
			double old = L_bin[NBIN_L-1];	
			int ibin = (int) ((L_lines[iline][j] - L_min) / dL);
//			printf("%lf\n", L_lines[iline][j]);
//			printf("ibin: %d\n", ibin);
			/*
			if (L_bin[NBIN_L-1] != 44.25)
			{
	printf("***%lf\n", L_bin[NBIN_L-1]);
	printf("%d\n", j);
	break;
		}
			*/
				double new = L_bin[NBIN_L-1];
			if (old != new) printf("old\n");

			if (ibin >= NBIN_L)
			{
				big_count++;
				printf("%d\n", ibin);
			printf("%lf\n", L_lines[iline][j]);
			}
			else if (ibin < 0)
			{
				small_count++;
			}
			else
			{
				L[iline][ibin] += 1;
				L_bin[ibin] += L_lines[iline][j];
			}
		}
		printf("cp10\n");


		for (int ibin = 0; ibin < NBIN_L; ibin++)
		{
			L_bin[ibin] /= L[iline][ibin];
			printf("%d\n", L[iline][ibin]);
		}
		printf("%d galaxies have logL < logL_min= %lf\n", small_count, L_min);
		printf("%d galaxies have logL > logL_max= %lf\n", big_count, L_max);

	}



	char *linenames_with_noise[NLINES+1];
	for (int line = 0; line < NLINES; line++)
	{
		linenames_with_noise[line] = linenames[line];
	}
	linenames_with_noise[NLINES] = "noise";



	FILE *fp1;
	char fname1[128];
	sprintf(fname1, "./luminosity_function_emerge.txt");
	fp1 = fopen(fname1, "w");

	fprintf(fp1, "# log L [erg/s]");
	fprintf(fp1, " %-3s (z=%.1f) ", linenames[0], redshifts[0]);

	for (int iline = 1; iline < NLINES; iline++)
	{
		fprintf(fp1, " %-3s (z=%.1f) ", linenames[iline], redshifts[iline]);
	}
	//fprintf(fp1, "  %s", linenames_with_noise[NLINES]);
	fprintf(fp1, "\n");

	for (int ibin = 0; ibin < NBIN_L; ibin++)
	{
		fprintf(fp1, "%-15lf ", L_bin[ibin]);
//		printf("%lf ", L_bin[ibin]);
		for (int iline = 0; iline < NLINES; iline++)
		{
			fprintf(fp1, "%-13d ", L[iline][ibin]);
		}
		fprintf(fp1, "\n");
	}

	fclose(fp1);	

	
	time(&now);
	printf("%s\n", ctime(&now));

}
