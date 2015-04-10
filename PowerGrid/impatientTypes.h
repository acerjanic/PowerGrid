//
//  impatientTypes.h
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/8/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_impatientTypes_h
#define PowerGrid_impatientTypes_h
typedef struct{
    int numSamples;
    int imageSize[3];
    int gridSize[3];
    float gridOS;
    float kernelWidth;
    int binsize;
    int useLUT;
    int sync;
}parameters;

template <T1>
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
