//
//  PowerGridIsmrmrd.cpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 8/31/2015.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#include "PowerGrid.h" //Project headers.
#include "../Support/CeempleMatio.h" //Headers for using savemat and loadmat
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/dataset.h"
#include "ismrmrd/version.h"

#include <boost/program_options.hpp>

using namespace arma; //Armadillo stuff is in the arma namespace
using namespace std; //complex type comes from the STL
using namespace ISMRMRD;

namespace po = boost::program_options;

int main(int argc, char** argv)
{

	po::options_description desc("Allowed options");
	desc.add_options()
			("help,h", "produce help message");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}
	/*
	try
	{

		ISMRMRD::Dataset d(testpath.c_str(),"dataset", false);
		std::string xml;
		d.readHeader(xml);
		ISMRMRD::IsmrmrdHeader hdr;
		ISMRMRD::deserialize(xml.c_str(),hdr); 

	}
	catch (const xml_schema::exception& e)
	{
		cerr << e << endl;

		return 1;
	}
	int test = test_SpeedCompareGgrid<double>(testPath, Nx, Ny, Nz, NL, Niter, Ncoils, Nshots, Beta);
	*/

	return 0;
}

