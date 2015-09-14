//
//  PowerGridDWI3D.cpp
//  PowerGrid
//
//  Created by Joe Holtrop on 3/12/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#include "PowerGrid.h" //Project headers.
#include "../Support/CeempleMatio.h" //Headers for using savemat and loadmat


using namespace arma; //Armdillo stuff is in the arma namespace
using namespace std; //complex type comes from the STL


int main(int argc, char** argv)
{
    uword Nx,Ny,Nz,Niter = 1,NL = 1,Ncoils;
    uword startIndex, endIndex;

    string testPath,configPath;
    if (argc > 1) {
        testPath = std::string(argv[1]);
        configPath = testPath+"config.xml";
    } else {
      cout << "Enter a path to find test files." << endl;
      return -1;
    }

    if (argc > 2) {
        startIndex = atoi(argv[2]);
        endIndex = atoi(argv[3]);
        cout << "Start Index." << startIndex << "End Index" << endIndex << endl;

        try
        {
            auto_ptr<PowerGridConfig_t> cfg(PowerGridConfig(configPath.c_str()));

            Nx = cfg->Nx();
            Ny = cfg->Ny();
            Nz = cfg->Nz();
            NL = cfg->Ntimeseg();
            Niter = cfg->Niter();
            Ncoils = cfg->Ncoils();

        }
        catch (const xml_schema::exception& e)
        {
            cerr << e << endl;
            return 1;
        }

        int test = reconfMRIGgrid<double>(testPath, Nx,Ny,Nz,NL,Niter,Ncoils,startIndex,endIndex);
    } else {

        try
        {
            auto_ptr<PowerGridConfig_t> cfg(PowerGridConfig(configPath.c_str()));

            Nx = cfg->Nx();
            Ny = cfg->Ny();
            Nz = cfg->Nz();
            NL = cfg->Ntimeseg();
            Niter = cfg->Niter();
            Ncoils = cfg->Ncoils();

        }
        catch (const xml_schema::exception& e)
        {
            cerr << e << endl;
            return 1;
        }

        //int test = test_SpeedCompareGgrid<double>(testPath, Nx,Ny,Nz,NL,Niter,Ncoils);
		int test = test_DWI<double>(testPath, Nx,Ny,Nz,NL,Niter,Ncoils);

    }

    return 0;
}

