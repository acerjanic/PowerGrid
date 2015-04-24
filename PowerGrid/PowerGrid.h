//
//  PowerGrid.h
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/2/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_PowerGrid_h
#define PowerGrid_PowerGrid_h
#define ARMA_NO_DEBUG
#include <iostream>
#include <cmath>
#include "armadillo"
#include "../Support/CeempleComplex.h"
#include "../Support/CeempleArmadillo.h"
#include "../Support/CeempleMatio.h"
#include "test_pwls_pcg.h"
#include "impatientTypes.h"
#include "fftw3.h" //for fft library
#include "config.hxx"
#include <memory>

// I throw everything into namespace arma to use the Armadillo objects and functions without arma:: everywhere
// Maybe I should reevaluate this later...
namespace arma {
#include "op_circshift_bones.hpp"
#include "op_circshift_meat.hpp"
#include "fn_circshift.hpp"
#include "fftshift.hpp"
#include "Gfft.hpp"
#include "SENSE.hpp"
#include "QuadPenalty.hpp"
#include "pwls_pcg1.hpp"
#include "ftCpu.hpp"
#include "Gdft.hpp"
#include "test_gdft.hpp" 
#include "griddingSupport.hpp"
#include "gridding.hpp"
#include "Ggrid.hpp"
#include "test_ggrid.hpp"
#include "FieldCorrection.hpp"
#include "test_FieldCorrection.hpp"
#include "test_3D.hpp"
#include "test_SpeedCompareGgrid.hpp"
#include "test_SpeedCompareGdft.hpp"

}

#endif
