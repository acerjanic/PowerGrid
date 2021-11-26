/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [gridding.cpp]

    Synopsis    [Implementation of the forward and adjoint non-uniform Fast
                                Fourier Transform (NUFFT) on CPU and GPU via
 OpenACC.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#include "gridding.h"


// Explicit Instantiations
/*
template int Gnufft<float>::gridding_adjoint_2D(unsigned int, parameters<float>, float,
    ReconstructionSample<float>*,
    const float*, const uword,
    float*) const;
template int Gnufft<double>::gridding_adjoint_2D(unsigned int, parameters<double>,
    double, ReconstructionSample<double>*,
    const double*, const uword,
    double*) const;
template int Gnufft<float>::gridding_adjoint_3D(unsigned int, parameters<float>, float,
    ReconstructionSample<float>*,
    const float*, const uword,
    float*) const;
template int Gnufft<double>::gridding_adjoint_3D(unsigned int, parameters<double>,
    double, ReconstructionSample<double>*,
    const double*, const uword,
    double*) const;

template int Gnufft<float>::gridding_forward_2D(unsigned int, parameters<float>,
    const float*, const float*,
    float beta, float*,
    const float*, const uword,
    float*) const;
template int Gnufft<double>::gridding_forward_2D(unsigned int, parameters<double>,
    const double*, const double*,
    double beta, double*,
    const double*, const uword,
    double*) const;
template int Gnufft<float>::gridding_forward_3D(unsigned int, parameters<float>,
    const float*, const float*,
    const float*, float beta,
    float*, const float*,
    const uword, float*) const;
template int Gnufft<double>::gridding_forward_3D(unsigned int, parameters<double>,
    const double*, const double*,
    const double*, double beta,
    double*, const double*,
    const uword, double*) const;
template void Gnufft<float>::computeFH_CPU_Grid(int, const float*, const float*,
    const float*,
    const float*, int, int, int,
    float gridOS,
    const float, const float, const float*,
    const uword, void*, CFTHandle*, float*,
    float*, float*, float*) const;
template void Gnufft<double>::computeFH_CPU_Grid(int, const double*, const double*,
    const double*,
    const double*, int, int, int,
    double gridOS,
    const double, const double,
    const double*, const uword, void*, CFTHandle*,
    double*, double*, double*, double*) const;
template void Gnufft<float>::computeFd_CPU_Grid(int, const float*,
    const float*, const float*, const float*,
    int, int, int, float, const float, const float,
    const float*, const uword, void*, CFTHandle*,
    float*, float*, float*, float*, float*) const;
template void Gnufft<double>::computeFd_CPU_Grid(int, const double*,
    const double*, const double*, const double*,
    int, int, int, double, const double, const double,
    const double*, const uword, void*, CFTHandle*,
    double*, double*, double*, double*, double*) const;

*/
