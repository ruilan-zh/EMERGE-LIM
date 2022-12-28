#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>

#include "proto.h"
#include "constant.h"
#include "para.h"

#define NLINES 1

int main(int argc, char *argv[])
{
	time_t now;
	time(&now);
	printf("%s\n", ctime(&now));

	Cosmology cosm;
	cosm.omega = 0.3153;
	cosm.lambda = 0.6847;
	cosm.omegab = 0.0493;
	cosm.hubble = 0.6736;

	printf("cp0\n");
	if (argc != 5)
	{
		fprintf(stderr, "ERROR: incorrect number of arguments to %s\n", argv[0]);
		fprintf(stderr, "	Input: ");
		for(int i = 0; i < argc; i++) fprintf(stderr, "%s ", argv[i]);
		fprintf(stderr, "\n");
		fprintf(stderr, "	Usage: %s dir pixlen boxsize Nmeshz\n", argv[0]);

		exit(1);
	}
	printf("cp1\n");
	char *dir;
	dir = argv[1];

	double z0 = 1.47;
	double d_A = angular_distance(cosm, z0); // [cm]

	double pixlen = atof(argv[2]); // [arcsec]
	pixlen = pixlen/(60*60);//degrees
	double pixlen_rad = pixlen * 2 * M_PI / 360;
	double pixlen_cm = pixlen_rad * (d_A * (1+z0));
	double pixlen_mpc = pixlen_cm / (mpc);
	double pixlen_mpc_h = pixlen_mpc * cosm.hubble; 

	double boxsize_mpc_h = atof(argv[3]); //[Mpc/h]
//	double boxsize_mpc_h = boxsize_mpc * cosm.hubble; //[Mpc/h]
	double boxsize_mpc = boxsize_mpc_h /cosm.hubble; // [Mpc]

	double boxsize_cm = boxsize_mpc * mpc; // [cm]
	double boxsize_rad = boxsize_cm / (d_A * (1+z0)); // [deg]
	double boxsize_deg = boxsize_rad / (2 * M_PI / 360);	

	int Nmeshxy = boxsize_mpc_h/pixlen_mpc_h;
	int Nmeshz = atoi(argv[4]);
	if (Nmeshz == 1) ;
	else Nmeshz = Nmeshxy;
	printf("Nmeshz: %d\n", Nmeshz);
	int Nmeshtot = Nmeshxy*Nmeshxy*Nmeshz;

	printf("cp2\n");
	char *linenames_with_noise[NLINES+1] = {"Ha", "OIII", "OII", "noise"};
	linenames_with_noise[-1] = "noise";
	double lambda_rest_Ha = 6563 * Angstrom;
	double lambda_rest_O3 = 5007 * Angstrom;
	double lambda_rest_O2 = 3727 * Angstrom;
	double lambda_rest[NLINES] = {lambda_rest_Ha, lambda_rest_O3, lambda_rest_O2};

	double lambda_obs = 1.5 * 10000 * Angstrom;
	double redshifts[NLINES];
	for (int iline = 0; iline < NLINES; iline++)
	{
		redshifts[iline] = lambda_obs/lambda_rest[iline] - 1;
	}

	FILE *fp;
	char fname[128];
	sprintf(fname, "%s/intensity.txt", dir);

	fp = fopen(fname, "r");
	if (fp==NULL)
	{
		fprintf(stderr, "cannot open file: %s\n", fname);
		exit(1);
	}

	printf("cp3\n");
	double** nu_I_lines;
	nu_I_lines = (double**)malloc((NLINES+1) * sizeof(double*));
	for (int line = 0; line < NLINES+1; line++)
	{
		nu_I_lines[line] = (double*) malloc(Nmeshtot*sizeof(double));
	}

	printf("cp4\n");
	char str[256];
	int i = 0;

	time(&now);
	printf("%s\n", ctime(&now));

	printf("Reading %s\n", fname);
	fflush(stdout);

	double max_nu_I = 0;
	double min_nu_I = 0;
	while(fgets(str, 256, fp) != NULL)
	{
		if (str[0] == '#') continue;
		else
		{
			sscanf(str, "%lf %lf %lf %lf ", &nu_I_lines[0][i], &nu_I_lines[1][i], &nu_I_lines[2][i], &nu_I_lines[3][i]);
			if (nu_I_lines[0][i] > max_nu_I) max_nu_I = nu_I_lines[0][i];
			//if (nu_I_lines[NLINES][i] < min_nu_I) min_nu_I = nu_I_lines[NLINES][i];
			//if (nu_I_lines[1][i] > max_nu_I) max_nu_I = nu_I_lines[1][i];
		
			i++;
		}
	}
	fclose(fp);

	printf("cp5\n");

	double** nu_I_lines_2d;
	nu_I_lines_2d = (double**)malloc((NLINES+1) * sizeof(double*));
	for (int line = 0; line < NLINES+1; line++)
	{
		nu_I_lines_2d[line] = (double*) malloc(Nmeshxy*Nmeshxy*sizeof(double));
	}

	printf("cp6\n");	
	time(&now);
	printf("%s\n", ctime(&now));

	printf("Computing projected map\n");
	fflush(stdout);

	int Nslice = Nmeshz;
	double max_nu_I_2d[NLINES+1] = {0.0, 0.0, 0.0, 0.0};
	for (int line = 0; line < NLINES+1; line++)
	{
		max_nu_I_2d[line] = intensity_proj_map(linenames_with_noise[line], nu_I_lines_2d[line], nu_I_lines[line], Nmeshxy, Nslice, dir);
	}

	printf("max 2d: %lf\n", max_nu_I_2d[0]);
	printf("max 2d: %lf\n", max_nu_I_2d[1]);
	printf("max 3d: %lf\n", max_nu_I);

	//power spectrum//
	
	char fnameps2d[128], fnameps3d[128];
    sprintf(fnameps2d,  "%s/2dps.txt", dir);
	sprintf(fnameps3d,  "%s/3dps.txt", dir);


	/*
	printf("cubesize: %lf Mpc/h\n", slice_depth);
	double ratio = slice_depth/boxsize;

	printf("cubesize/boxsize: %lf\n", ratio);
	printf("Nmeshxy for box: %d\n", Nmeshxy);
	double Nmesh_cube1 = Nmeshxy*ratio;
	printf("Nmeshxy for cube: %lf\n", Nmesh_cube1);
	int Nmesh_cube = (int)round(Nmeshxy*ratio);
	printf("Nmeshxy for cube: %d\n", Nmesh_cube);
	*/

	/*
	double** nu_I_lines_3d;
	nu_I_lines_3d = (double**)malloc((NLINES+1) * sizeof(double*));
	for (int line = 0; line < NLINES+1; line++)
	{
		nu_I_lines_3d[line] = (double*) malloc(Nmesh_cube*Nmesh_cube*sizeof(double));
	}
	

	for (int iline = 0; iline < NLINES+1; iline++)
	{
		for (int iy = 0; iy < Nmesh_cube; iy++)
		{
			for (int ix = 0; ix < Nmesh_cube; ix++)
			{
				nu_I_lines_3d[iline][ix + iy*Nmesh_cube] = nu_I_lines[iline][ix + iy*Nmeshxy];
			}
		}
	}
	*/

	int Npix[3] = {Nmeshxy, Nmeshxy, Nmeshz};
//	int Npix_cube[3] = {Nmesh_cube, Nmesh_cube, Nmeshz};
//	int Nmesh_cube_tot = Nmesh_cube*Nmesh_cube*Nmeshz;

	printf("cp15.5\n");

	time(&now);
	printf("%s\n", ctime(&now));

	printf("Computing 2D power spectrum\n");
	fflush(stdout);
	ps_2d(fnameps2d, Npix, nu_I_lines_2d, linenames_with_noise, boxsize_mpc_h, NLINES+1);
	ps_3d(fnameps3d, Npix, nu_I_lines, linenames_with_noise, boxsize_mpc_h, NLINES+1);	
	time(&now);
	printf("%s\n", ctime(&now));

	//probability distribution function//
	
	char fnamepdf_2d[128], fnamepdf_3d[128];
	sprintf(fnamepdf_2d, "%s/pdf_2d.txt", dir);
	sprintf(fnamepdf_3d, "%s/pdf_3d.txt", dir);
	printf("cp16.01\n");

	double ps_shot_2d[NLINES+1] = {};
	double ps_shot_3d[NLINES+1] = {};

	printf("Computing PDF\n");
	fflush(stdout);

	pdf_nd(2, fnamepdf_2d, nu_I_lines_2d, linenames_with_noise, max_nu_I_2d[0], Npix, NLINES+1);
//	pdf_with_shot_nd(2, fnamepdf_2d, nu_I_lines_2d, linenames_with_noise, max_nu_I_2d[1], Npix, NLINES+1, ps_shot_2d);
//	pdf_with_shot_nd(3, fnamepdf_3d, nu_I_lines, linenames_with_noise, max_nu_I, Npix, NLINES+1, ps_shot_3d);
	printf("cp18\n");
	time(&now);
	printf("%s\n", ctime(&now));

	double pixsize = pixlen*pixlen/pow(180/M_PI,2); //convert to steradians
	printf("cp19\n");


	printf("cp20\n");
//	double V_voxel = pixsize * dnu(lambda_obs, 1, z, deltaz);
	

//	printf("pixsize: %lf\n", pixsize);
//	printf("dnu: %lf\n", dnu(lambda_obs, 1, z, deltaz));

	printf("cp21\n");


	char fnameshot2d[128], fnameshot3d[128];
    sprintf(fnameshot2d,  "%s/shot2d.txt", dir);
	sprintf(fnameshot3d,  "%s/shot3d.txt", dir);



//	ps_2d_with_shot(fnameps2d, fnameshot2d, Npix, nu_I_lines_2d, linenames_with_noise, boxsize, NLINES+1, ps_shot_2d, V_voxel);	


	printf("cp22\n");
//	ps_3d_with_shot(fnameps3d, fnameshot3d, Npix, nu_I_lines, linenames_with_noise, boxsize, NLINES+1, ps_shot_3d);	


	//pdf noise linear//
	printf("%lf\n", min_nu_I);
	printf("%lf\n", max_nu_I);
	char fnamepdf_linear[128];
	sprintf(fnamepdf_linear, "%s/pdf_linear.txt", dir);
	pdf_linear(fnamepdf_linear, nu_I_lines, linenames_with_noise, min_nu_I, max_nu_I, Nmeshtot, NLINES+1);
	
	time(&now);
	printf("%s\n", ctime(&now));		
}
