//
//  SENSE.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/4/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_SENSE_hpp
#define PowerGrid_SENSE_hpp

// We are using two template types at the moment. One for the type of data to be processed (ie Col<cx_double>) and one for the type of G object (ie Gfft<Col<cx_double>>
template<typename T1, typename Tobj>
class SENSE {
    typedef complex <T1> CxT1;
public:
    SENSE();

    //Class variables go here
    uword n1 = 0; //Data size
    uword n2 = 0; //Image size
    uword nc = 0; //number of coils
    Tobj* G_obj;
    Mat <CxT1> SMap; //dimensions Image size b (n1 by number of coils (nc)


    //Class constructor
    SENSE(Tobj& G, Col <CxT1> SENSEmap, uword a, uword b, uword c)
    {
	    n1 = a;
	    n2 = b;
	    nc = c;
	    G_obj = &G;
	    SMap = reshape(SENSEmap, n2, nc);
    }

    //Overloaded operators go here

    //Forward transformation is *
    // d is the vector of data of type T1, note it is const, so we don't modify it directly rather return another vector of type T1
    Col <CxT1> operator*(const Col <CxT1>& d) const
    {

	    Mat <CxT1> outData = zeros<Mat<CxT1 >> (this->n1, this->nc);
	    //Col<CxT1> temp;
	    //In SENSE we store coil data using the columns of the data matrix, and we weight the data by the coil sensitivies from the SENSE map
	    for (unsigned int ii = 0; ii<this->nc; ii++) {
		    //temp = d%(this->SMap.col(ii));
		    outData.col(ii) = (*this->G_obj)*(d%(this->SMap.col(ii)));

	    }
	    Col <CxT1> out = vectorise(outData);
	    //equivalent to returning col(output) in MATLAB with IRT
	    return out;

    }

    //For the adjoint operation, we have to weight the adjoint transform of the coil data by the SENSE map.
    Col <CxT1> operator/(const Col <CxT1>& d) const
    {

	    Mat <CxT1> inData = reshape(d, this->n1, this->nc);

	    Col <CxT1> outData = zeros<Col<CxT1 >> (this->n2);

	    for (unsigned int ii = 0; ii<this->nc; ii++) {

		    outData += conj(this->SMap.col(ii))%((*this->G_obj)/inData.col(ii));
		    // temp = ((*this->G_obj)/data);
		    // savemat("/shared/mrfil-data/data/PowerGridTest/64_64_16_4coils/coil_img_" + std::to_string(ii) + ".mat","img",temp);

	    }

	    //equivalent to returning col(output) in MATLAB with IRT
	    return vectorise(outData);


    }
};

#endif
