
//
//  PowerGrid.h
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/2/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_PowerGrid_h
#define PowerGrid_PowerGrid_h
//#define ARMA_NO_DEBUG // Disable this comment only for release.
#include <iostream>
#include <cmath>
#include <string>

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

//Support Headers for making it easier to work with Armadillo and Matio.
//#include "../Support/CeempleComplex.h"
#include "../Support/ArmaExtensions/arma_extensions.h"

#include "../Support/CeempleArmadillo.h"
#include "../Support/CeempleMatio.h"

#include "impatientTypes.h"
#include "config.hxx"
#include <memory>

//Headers for ISMRMRD Support
//#include "ismrmrd/ismrmrd.h"
//#include "ismrmrd/xml.h"
//#include "ismrmrd/dataset.h"
//#include "ismrmrd/version.h"



//
namespace PowerGrid {

#include "fftshift.hpp"
#include "Gfft.hpp"
#include "ftCpu.hpp"
#include "Gdft.hpp"
#include "fftGPU.hpp"
#include "griddingSupport.hpp"
#include "gridding.hpp"
#include "Ggrid.hpp"
#include "SENSE.hpp"
#include "Robject.hpp"
#include "QuadPenalty.hpp"
#include "TVPenalty.hpp"
#include "solve_pwls_pcg.hpp"
#include "reconSolve.hpp"
#include "test_gdft.hpp"
#include "test_ggrid.hpp"
#include "FieldCorrection.hpp"
#include "DWICGMC.hpp"
#include "DWICGMCDFT.hpp"
#include "mpiDWICGMCDFT.hpp"
#include "test_FieldCorrection.hpp"
#include "test_3D.hpp"
#include "test_DWI.hpp"
#include "test_DWIDft.hpp"
#include "test_SpeedCompareGgrid.hpp"
#include "test_SpeedCompareGdft.hpp"
#include "test_pwls_pcg.hpp"
#include "reconfMRIGdft.hpp"
#include "reconfMRIGgrid.hpp"

}

#endif
