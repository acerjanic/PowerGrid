/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [PowerGrid.h]

    Synopsis    [Main header file for the project.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/

/// @file PowerGrid.h
/// @brief Main project umbrella header: includes all encoding operators, regularizers,
///        solvers, ISMRMRD support, and (optionally) MPI-distributed variants.

#ifndef PowerGrid_PowerGrid_h
#define PowerGrid_PowerGrid_h
//#define ARMA_NO_DEBUG // Disable this comment only for release.

#include "PGIncludes.h"

// Headers for ISMRMRD Support
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/dataset.h"
#include "ismrmrd/version.h"
#include "IO/acqTracking.h"
//namespace PowerGrid {
#include "Solvers/Robject.h"
#include "Solvers/TVPenalty.h"

#include "Operators/Gdft.h"
#include "Operators/Gfft.h"
#include "Operators/Gnufft.h"

#include "Operators/pcSENSE.h"
#include "Operators/pcSenseTimeSeg.h"
#include "Solvers/solve_pwls_pcg.hpp"

#include "Operators/SENSE.h"
#include "Gridding/TimeSegmentation.h"
#include "FFT/fftGPU.h"
#include "IO/fftshift.hpp"
#include "FFT/ftCpu.h"
#include "Gridding/gridding.h"
#include "Gridding/griddingSupport.h"

#include "Solvers/reconSolve.h"

#ifdef PowerGridMPI

#include "MPI/mpipcSENSE.h"
#include "MPI/mpipcSENSETimeSeg.h"
#endif // PowerGridMPI

//}

#endif
