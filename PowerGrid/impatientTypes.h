//
//  impatientTypes.h
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/8/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_impatientTypes_h
#define PowerGrid_impatientTypes_h
template<typename T1>
struct parameters{
    int numSamples;
    int imageSize[3];
    int gridSize[3];
    T1 gridOS;
    T1 kernelWidth;
    int binsize;
    int useLUT;
    int sync;
};

template <typename T1>
struct ReconstructionSample{
    T1 real;
    T1 imag;
    T1 kX;
    T1 kY;
    T1 kZ;
    T1 sdc;
    T1 t;
    T1 dummy;
};

#endif
