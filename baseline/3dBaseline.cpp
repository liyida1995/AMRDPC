using namespace std;
#include <cfloat>
#include <cmath>
#include <iostream>     // std::cout   
#include <fstream>      // std::ifstream
#include <string>		// std::string
#include <sstream>		// std::stringstream
#include <vector>
#include <map>
#include <algorithm>
#include "sz.h"
#include "rw.h"
#include <chrono>
#include <bits/stdc++.h>

double *AMR_compress(double *oriData, char* cfgFile, size_t blksize, size_t length4d, double eb, char* zipFilePath, size_t* comSize)
{	
    int status = 0;
    printf("cfgFile=%s\n", cfgFile); 
    status = SZ_Init(cfgFile);
    size_t outSize;
	unsigned char* bytes = SZ_compress_args(SZ_DOUBLE, oriData, &outSize, ABS, eb, 1, 1, 0, length4d, blksize, blksize, blksize);
    writeByteData(bytes, outSize, zipFilePath, &status);
	if(status!=SZ_SCES)
    {
	printf("Error: file %s cannot be written!\n", zipFilePath);
	exit(0);
    }
	std::cout << "compressedFile=" << zipFilePath << std::endl;
	*comSize += outSize;
	double *deData = (double *)SZ_decompress(SZ_DOUBLE, bytes, outSize, 0, length4d, blksize, blksize, blksize);
 	size_t i;
    double Max, Min, diffMax, err, maxpw_relerr = 0, relerr;
    Max = oriData[0];
    Min = oriData[0];
    diffMax = fabs(deData[0] - oriData[0]);

	if (length4d == 0){
		length4d = 1;
	}
    for (i = 0; i < length4d*blksize*blksize*blksize; i++)
    {
    	if (Max < oriData[i]) Max = oriData[i];
    	if (Min > oriData[i]) Min = oriData[i];
		err = fabs(deData[i] - oriData[i]);
    	if (diffMax < err)
    		diffMax = err;
        if(oriData[i]!=0)
        {
			relerr = err/fabs(oriData[i]);
			if(maxpw_relerr<relerr)
				maxpw_relerr = relerr;
        }
    }

    //printf ("Max absolute error = %.20G\n", diffMax);
    //printf ("Max relative error = %.20G\n", diffMax/(Max-Min));
    //printf ("Max pw_relative err = %.20G\n", maxpw_relerr);

    //====================PSNR============================
    double vrange = Max - Min;
    double mse, temp;
    size_t length = length4d * blksize * blksize * blksize;
    for(i = 0; i < length; i++){
        temp += (oriData[i] - deData[i]) * (oriData[i] - deData[i]);
    }
    mse = 20 * log10(vrange) - 10 * log10(temp / length);
    printf("The value of PSNR: %lf\n", mse);

   //====================SSIM=============================
   double ori_mean = 0;
   double de_mean = 0;
   double ori_sum = 0;
   double de_sum = 0;
   for(i = 0; i < length; ++i){
       ori_sum += oriData[i];
       de_sum += deData[i];
   }

   ori_mean = ori_sum / length;
   de_mean = de_sum / length;

   double ori_temp = 0;
   double de_temp = 0;
   double cov_temp = 0;
   for(i = 0; i < length; ++i){
        ori_temp += (oriData[i] - ori_mean) * (oriData[i] - ori_mean);
        de_temp += (deData[i] - de_mean) * (deData[i] - de_mean);
        cov_temp += (oriData[i]- ori_mean) * (deData[i] - de_mean);
   }

   ori_temp = ori_temp / length;
   de_temp = de_temp / length;
   double cov = 0;
   double ori_var = 0;
   double de_var = 0;
   cov = cov_temp / (length - 1);
   ori_var = ori_temp;
   de_var = de_temp;
   double ssim = 0;
   ssim = ((2 * ori_mean * de_mean) * 2 * cov) / ((ori_mean * ori_mean + de_mean * de_mean) * (ori_var + de_var));
   printf("The value of SSIM: %lf\n", ssim);


   //===================NRMSE=============================
   double nrmse;
   nrmse = (sqrt(temp / length)) / vrange;
   printf("The value of NRMSE: %lf\n", nrmse);


    if(status!=SZ_SCES)
    {
	printf("Error: file %s cannot be written!\n");
	free(oriData);
	exit(0);
    }
    free(oriData);
    free(bytes);
    printf("done\n");
    
    SZ_Finalize();
    return deData;
}


int main(int argc, char *argv[]) {
    char inName[640];
	sprintf(inName, "%s", argv[1]);
    std::ifstream f;
    f.open(inName, std::ifstream::binary);
	if (f.fail()) {
		std::cout    << "Error opening file" << std::endl;
		return 0;
	}

    f.seekg(0, std::ios::end);
	size_t file_len = f.tellg();
	f.seekg(0, std::ios::beg);

	size_t num_cells = file_len / (6 * sizeof(double)); // changes depending on # of fields
	std::cout    << "Number of cells in the file: " << num_cells << std::endl;


	double* xx = new double[num_cells];
	double* yy = new double[num_cells];
	double* zz = new double[num_cells];
    double* oriVol = new double[num_cells];
	double* vol = new double[num_cells];
	double* rho = new double[num_cells];
	double* temp = new double[num_cells];

	double rho_c = 2.77536627e11;
	double omegab = 0.045;
	double hubble = 0.7;
	double rhomean = rho_c * omegab * hubble * hubble;

    double maxVol = -DBL_MAX;
    double minVol = DBL_MAX;
	double maxRho = -DBL_MAX;
    double minRho = DBL_MAX;
	double max[3] = { -FLT_MAX,-FLT_MAX,-FLT_MAX };
	double range = 0;

	for (size_t i = 0; i < num_cells; ++i) 
	{
		f.read( reinterpret_cast<char*>(&xx[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&yy[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&zz[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&oriVol[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&rho[i]) , sizeof(double));
		if (rho[i] > maxRho)
			maxRho = rho[i];
		if (rho[i] < minRho)
			minRho = rho[i];
        // std::cout    << rho[i]<< std::endl;
		f.read( reinterpret_cast<char*>(&temp[i]) , sizeof(double));
        vol[i] = std::cbrt(oriVol[i]);
        if (vol[i] > maxVol)
			maxVol = vol[i];
		if (vol[i] < minVol)
			minVol = vol[i];

		if (xx[i] > max[1])
			max[1] = xx[i];

		if (yy[i] > max[2])
			max[2] = yy[i];

		if (zz[i] > max[0])
			max[0] = zz[i];
	}
	range = maxRho - minRho;
	std::cout << "range: " << range << std::endl;

	size_t grid = std::max(std::max(max[1], max[2]), max[3]);
	grid = grid/minVol;
	grid = pow(2, ceil(log2(grid)));

	std::cout << "output dims: " << grid << std::endl;

    // size_t intdims[3]; intdims[0] = 511; intdims[1] = 511; intdims[2] = 511;

	double * densgrid = new double[grid*grid*grid]{0.0};
	bool *mark  = new bool[grid*grid*grid]{0};
    double * degrid = new double[grid*grid*grid]{0.0};
    double * diffgrid = new double[grid*grid*grid]{0.0};

	//to 3d
	for (size_t i = 0; i < num_cells; ++i)
	{
		size_t xi, yi, zi, cpos, mx, my, mz;
		// Debug: Make sure 'round' is working right
		double dxi = (xx[i] / minVol);
		if (round(dxi) != dxi)
		{
			std::cout.precision(10);
			if( abs(dxi - round(dxi)) > 0.5)
				cout << "WARNING: Potential mapping error: " << dxi << " vs " << round(dxi) << std::endl;
		}

		xi = round(xx[i] / minVol);
		yi = round(yy[i] / minVol);
		zi = round(zz[i] / minVol);
		
		if (vol[i] == minVol)
		{
			cpos = xi + (yi * (grid)) + (zi * (grid) * (grid));
			densgrid[cpos] = rho[i];
            degrid[cpos] = rho[i];
		}
		else if (vol[i] == maxVol){
        mx = xi + (int)round(vol[i] / minVol);
		my = yi + (int)round(vol[i] / minVol);
		mz = zi + (int)round(vol[i] / minVol);
        for(size_t a = xi; a < mx; a++)
			for(size_t b = yi; b < my; b++)
				for (size_t c = zi; c < mz; c++)
                {
                    cpos = a + (b * (grid)) + (c * (grid) * (grid));
                    if (densgrid[cpos] == 0)
                    {
                        densgrid[cpos] = rho[i];
                        degrid[cpos] = rho[i];
                    }
                }
        }
        else
        {
            mx = xi + (int)round(vol[i] / minVol);
            my = yi + (int)round(vol[i] / minVol);
            mz = zi + (int)round(vol[i] / minVol);
            for (size_t a = xi; a < mx; a++)
				for (size_t b = yi; b < my; b++)
					for (size_t c = zi; c < mz; c++)
					{
                        cpos = a + (b * (grid)) + (c * (grid) * (grid));
                        densgrid[cpos] = rho[i];
                        degrid[cpos] = rho[i];
                    }
        }
        
	}


	double eb = atof(argv[3]);
    char *cfg = argv[2];
    char path[650];
	char temp_path[660];
	sprintf(path, "%s.3dbs", inName);

	size_t outSize = 0;
	degrid = AMR_compress(degrid, cfg, grid, 0, eb, path, &outSize);
	std::cout << "3dbs compression ratio : " << double(num_cells) * 8 / outSize << std::endl;
	std::cout << "3dbs bitrate : " << double(outSize) * 8 / num_cells << std::endl;

	std::cout << "3dBaseline compress finish!" << std::endl;
        std::cout << "==================================" << std::endl;

    delete [] degrid;
    delete [] densgrid;
    delete [] diffgrid;

    return 0;
}

