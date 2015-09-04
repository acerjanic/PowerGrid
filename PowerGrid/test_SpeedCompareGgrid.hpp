//
//  test_SpeedCompareGgrid.h
//  PowerGrid
//
//  Created by Joe holtrop on 4/7/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__SpeedCompare__Ggrid
#define __PowerGrid__SpeedCompare__Ggrid

#include "PowerGrid.h"

using namespace arma;

template<typename T1, typename T2>
void cast_data(T1 in, T2 out)
{
	for (int ii = 0; ii<in.n_rows; ii++) {
		out(ii) = in(ii);
	}

};

template<typename T1>
int test_SpeedCompareGgrid(string dataPath, uword Nx, uword Ny, uword Nz, uword L, uword niter, uword nc, uword nshots,
		T1 beta)
{
	string testPath = dataPath;
	typedef complex <T1> CxT1;
	//uword Nx =64;
	//uword Ny =64;
	//uword Nz =16;

	//Setup image space coordinates/trajectory
	Cube <T1> ix(Nx, Ny, Nz);
	Cube <T1> iy(Nx, Ny, Nz);
	Cube <T1> iz(Nx, Ny, Nz);
	Col <T1> FM;
	Col <CxT1> SMap;

	ix.zeros();
	iy.zeros();
	iz.zeros();

	//generate the image space coordinates of the voxels we want to reconstruct
	// after vectorizing ix and iy the image coordinates must match the Field and SENSe map image coordinates
	for (uword ii = 0; ii<Ny; ii++) { //y
		for (uword jj = 0; jj<Nx; jj++) { //x
			for (uword kk = 0; kk<Nz; kk++) { //z

				ix(ii, jj, kk) = ((T1) jj-(T1) Nx/2.0)/((T1) Nx);
				iy(ii, jj, kk) = ((T1) ii-(T1) Ny/2.0)/((T1) Ny);
				iz(ii, jj, kk) = ((T1) kk-(T1) Nz/2.0)/((T1) Nz);
			}
		}
	}

	Col <T1> kx;
	loadmat(testPath+"kx.mat", "kx", &kx);
	//Col <double> kx;
	//kx.copy_size(kx_tmp);
	//cast_data(kx_tmp,kx);
	Col <T1> ky;
	loadmat(testPath+"ky.mat", "ky", &ky);
	//Col <double> ky;
	//ky.copy_size(ky);
	//cast_data(ky_tmp,ky);
	Col <T1> kz;
	loadmat(testPath+"kz.mat", "kz", &kz);
	//Col <double> kz;
	//kz.copy_size(kz);
	//cast_data(kz_tmp,kz);
	//Col <T1> kz = kz_tmp;


	uword nro;
	nro = kx.n_elem;

	Col <T1> tvec;
	loadmat(testPath+"t.mat", "t", &tvec);
	//savemat(testPath+"twrite.mat", "twrite", tvec);
	//Col <T1> tvec;
	//tvec.copy_size(tvec_tmp);
	//cast_data(tvec_tmp,tvec);

	//Col <T1> FM;
	loadmat(testPath+"FM.mat", "FM", &FM);
	//FM.copy_size(FM_tmp);
	//cast_data(FM_tmp,FM);



	// Fourier transfrom operator
	//cout << "Initializing Ggrid" << endl;
	Ggrid<T1> Gg(nro, 2.0, Nx, Ny, Nz, kx, ky, kz, vectorise(ix), vectorise(iy), vectorise(iz));

	// Field correction operation

	uword type = 1; // 2 for min max time seg and 1 for Hanning
	//uword L = 4;
	//cout << "Initializing FieldCorrection" << endl;
	FieldCorrection<T1, Ggrid<T1>> A(Gg, vectorise(FM), vectorise(tvec), nro, Nx*Ny*Nz, L, type, nshots);

	//cout << "Initializing Gdft" << endl;
	//Gdft<T1,T2> Gd(nro,Nx*Ny*Nz,kx,ky,kz,vectorise(ix),vectorise(iy),vectorise(iz),vectorise(FM),vectorise(tvec));

	//uword nc = 4;
	//Col <CxT1> SMap;
	loadmat(testPath+"SMap.mat", "SMap", &SMap);
	//loadmat(testPath+"FM.mat", "FM", &FM_tmp);
	//Col <CxT1> SMap;
	//SMap.copy_size(SMap_tmp);
	//cast_data(SMap_tmp,SMap);

	//cout << "Iniitalizing SENSE gdft" << endl;
	//SENSE<cx_double, Gdft<T1,T2>> Sd(Gd,SMap,nro,Nx*Ny*Nz,nc);

	// Sense operation
	//cout << "Iniitalizing SENSE Ggrid" << endl;
	SENSE<T1, FieldCorrection<T1, Ggrid<T1>>> Sg(A, SMap, nro, Nx*Ny*Nz, nc);

	//cout << "loading data" << endl;
	Col <CxT1> data;
	loadmat(testPath+"data.mat", "data", &data);
	//Col <CxT1> data;
	//data.copy_size(data_tmp);
	//cast_data(data_tmp,data);
	savemat(testPath+"datawrite.mat", "data", data);

	// Variables needed for the recon: Penalty object, num of iterations
	//ucube ReconMask(Nx,Ny,Nz);
	//ReconMask.ones();

	//cout << "Iniitalizing QuadPenalty" << endl;
	QuadPenalty<T1> R(Nx, Ny, Nz, beta);
	//cout << "QuadPenalty setup successfull" << endl;

	//uword niter = 10;
	Col <CxT1> xinit(Nx*Ny*Nz); // initial estimate of x
	xinit.zeros();
	Col <T1> W;
	//W = eye<sp_mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is it right?
	W.ones(nro*nc);

	//Col<CxT1> x_t;
	//cout << "heading into solve_pwls_pcg" << endl;
	//x_t = solve_pwls_pcg<T1, SENSE<cx_double, FieldCorrection<T1, T2, Ggrid<T1,T2>>>,QuadPenalty<T1>>(xinit, S, W, data, R, niter);
	//x_t = Sg/data;
	//savemat(testPath+"test_3D.mat","img",x_t);
/*
    cout << "Runing adjoint with ggrid" << endl;
    Col<T1> x_Sg_adjoint;
    x_Sg_adjoint = Sg/data;
    savemat(testPath+"test_adjoint_Sggrid.mat","img",x_Sg_adjoint);

    cout << "Runing forward with ggrid" << endl;
    Col<T1> x_Sg_forward;
    x_Sg_forward = Sg*x_Sg_adjoint;
    savemat(testPath+"test_forward_Sggrid.mat","img",x_Sg_forward);

    cout << "Runing adjoint with gdft" << endl;
    Col<T1> x_Sd_adjoint;
    x_Sd_adjoint = Sd/data;
    savemat(testPath+"test_adjoint_Sdft.mat","img",x_Sd_adjoint);

    cout << "Runing forward with gdft" << endl;
    Col<T1> x_Sd_forward;
    x_Sd_forward = Sd*x_Sd_adjoint;
    savemat(testPath+"test_forward_Sdft.mat","img",x_Sd_forward);
*/
	//cout << "Runing pwls with ggrid" << endl;
	Col <CxT1> test_pwls;
	test_pwls = solve_pwls_pcg<T1, SENSE<T1, FieldCorrection<T1, Ggrid<T1>>>, QuadPenalty<T1>>(xinit, Sg, W, data, R,
			niter);
	savemat(testPath+"test_pwls.mat", "img", test_pwls);
/*
    cout << "Runing pwls with ggrid" << endl;
    Col<T1> test_pwls;
    test_pwls = solve_pwls_pcg<T1,  Ggrid<T1,T2>,QuadPenalty<T1>>(xinit, Gg, W, data, R, niter);
    savemat(testPath+"test_pwls.mat","img",test_pwls);
    */

	return 0;

}

#endif
