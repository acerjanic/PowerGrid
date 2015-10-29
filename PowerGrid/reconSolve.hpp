//
// Created by acerja2 on 10/20/15.
//

#ifndef POWERGRID_RECONSOLVE_HPP
#define POWERGRID_RECONSOLVE_HPP

using namespace arma;
using namespace PowerGrid;

template<typename T1>
void initImageSpaceCoords(Col <T1> &ix, Col <T1> &iy, Col <T1> &iz, uword Nx, uword Ny, uword Nz) {
//generate the image space coordinates of the voxels we want to reconstruct
// after vectorizing ix and iy the image coordinates must match the Field and SENSE map image coordinates
    Cube <T1> ixTemp(Nx, Ny, Nz);
    Cube <T1> iyTemp(Nx, Ny, Nz);
    Cube <T1> izTemp(Nx, Ny, Nz);

    for (uword ii = 0; ii < Ny; ii++) { //y
        for (uword jj = 0; jj < Nx; jj++) { //x
            for (uword kk = 0; kk < Nz; kk++) { //z

                ixTemp(ii, jj, kk) = ((T1) jj - (T1) Nx / 2.0) / ((T1) Nx);
                iyTemp(ii, jj, kk) = ((T1) ii - (T1) Ny / 2.0) / ((T1) Ny);
                izTemp(ii, jj, kk) = ((T1) kk - (T1) Nz / 2.0) / ((T1) Nz);
            }
        }
    }

    ix = vectorise(ixTemp);
    iy = vectorise(iyTemp);
    iz = vectorise(izTemp);
}

//Template parameters are  T1: data precision (double, float, FP16 etc...), TObj: Transform Object, RObj is regularization object
template<typename T1, typename TObj, typename RObj>
Col <complex<T1>> reconSolve(Col <complex<T1>> data, TObj &Sg, RObj &R, Col <T1> kx, Col <T1> ky, Col <T1> kz, uword Nx,
                             uword Ny, uword Nz, Col <T1> tvec, uword niter) {
    typedef std::complex <T1> CxT1;

    //Col<T1> ix, iy, iz;
    //initImageSpaceCoords(ix,iy,iz,Nz,Ny,Nz);

    //Data weighting term - use ones unless we have something better to use
    Col <T1> W;
    W.ones(data.n_rows);
    Col <std::complex<T1>> xinit;
    xinit.zeros(Nx * Ny * Nz);

    Col <CxT1> imageOut;
    imageOut = solve_pwls_pcg<T1, TObj, RObj>(xinit, Sg, W, data, R, niter);

    return imageOut;
}

#endif //POWERGRID_RECONSOLVE_HPP
