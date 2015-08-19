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

template<typename T1>
int test_SpeedCompareGgrid(string dataPath, uword Nx, uword Ny, uword Nz, uword L, uword niter, uword nc, uword nshots)
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

	//generate the image space cordinates of the voxels we want to reconstruct
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
	Col <T1> ky;
	loadmat(testPath+"ky.mat", "ky", &ky);
	Col <T1> kz;
	loadmat(testPath+"kz.mat", "kz", &kz);

	uword nro;
	nro = kx.n_elem;

	Col <T1> tvec;
	loadmat(testPath+"t.mat", "t", &tvec);

	loadmat(testPath+"FM.mat", "FM", &FM);
	//FM.zeros();



	// Fourier transfrom operator
	cout << "Initializing Ggrid" << endl;
	Ggrid<T1> Gg(nro, 2.0, Nx, Ny, Nz, kx, ky, kz, vectorise(ix), vectorise(iy), vectorise(iz));

	// Field correction operation

	uword type = 2; // 2 for min max time seg and 1 for Hanning
	//uword L = 4;
	cout << "Initializing FieldCorrection" << endl;
	FieldCorrection<T1, Ggrid<T1>> A(Gg, vectorise(FM), vectorise(tvec), nro, Nx*Ny*Nz, L, type, nshots);

	//cout << "Initializing Gdft" << endl;
	//Gdft<T1,T2> Gd(nro,Nx*Ny*Nz,kx,ky,kz,vectorise(ix),vectorise(iy),vectorise(iz),vectorise(FM),vectorise(tvec));

	//uword nc = 4;
	loadmat(testPath+"SMap.mat", "SMap", &SMap);

	//cout << "Iniitalizing SENSE gdft" << endl;
	//SENSE<cx_double, Gdft<T1,T2>> Sd(Gd,SMap,nro,Nx*Ny*Nz,nc);

	// Sense operation
	cout << "Iniitalizing SENSE Ggrid" << endl;
	SENSE<T1, FieldCorrection<T1, Ggrid<T1>>> Sg(A, SMap, nro, Nx*Ny*Nz, nc);

	cout << "loading data" << endl;
	Col <CxT1> data;
	loadmat(testPath+"data.mat", "data", &data);

	// Variables needed for the recon: Penalty object, num of iterations
	//ucube ReconMask(Nx,Ny,Nz);
	//ReconMask.ones();

	cout << "Iniitalizing QuadPenalty" << endl;
	TVPenalty<T1> R(Nx, Ny, Nz, 1e-5, 1e-7);
	cout << "QuadPenalty setup successfull" << endl;

	//uword niter = 10;
	Col <CxT1> xinit(Nx*Ny*Nz); // initial estimate of x
	xinit.zeros();
	Col <T1> W;
	//W = eye<sp_mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is it right?
	W = ones(nro*nc);

	//Col<T1> x_t;
	//cout << "heading into solve_pwls_pcg" << endl;
	//x_t = solve_pwls_pcg<T1, SENSE<cx_double, FieldCorrection<T1, T2, Ggrid<T1,T2>>>,QuadPenalty<T1>>(xinit, S, W, data, R, niter);
	//x_t = S/data;
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
	cout << "Runing pwls with ggrid" << endl;
	Col <CxT1> test_pwls;
	test_pwls = solve_pwls_pcg<T1, SENSE<T1, FieldCorrection<T1, Ggrid<T1>>>, TVPenalty<T1>>(xinit, Sg, W, data, R,
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
