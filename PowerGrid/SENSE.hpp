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
template <typename T1,typename Tobj>
class SENSE {
public:
    SENSE();
    
    //Class variables go here
    uword n1 = 0; //Data size
    uword n2 = 0; //Image size
    uword nc = 0; //number of coils
    Tobj *G_obj;
    Mat<T1> SMap; //dimensions Image size b (n1 by number of coils (nc)

    
    //Class constructor
    SENSE(Tobj &G, Col<T1> SENSEmap, uword a, uword b,uword c ) {
      n1 = a;
      n2 = b;
      nc = c;
      G_obj = &G;
      SMap = reshape(SENSEmap,n2,nc);
    }
    
    //Overloaded operators go here
    
    //Forward transformation is *
    // d is the vector of data of type T1, note it is const, so we don't modify it directly rather return another vector of type T1
    Col<T1> operator*(const Col<T1>& d) const {

      Mat<T1> outData = zeros<Mat<T1>>(this->n2,this->nc);

      for (unsigned int ii=0; ii < this->nc; ii++) {
        Col<T1> data = d%(this->SMap.col(ii));

        outData.col(ii) = (*this->G_obj)*(data);

      }

      //equivalent to returning col(output) in MATLAB with IRT
      return vectorise(outData);
    }
    
    Col<T1> operator/(const Col<T1>& d) const {

      Mat<T1> inData = reshape(d,this->n2,this->nc);

      Col<T1> outData = zeros<Col<T1>>(this->n2);

      for (unsigned int ii=0; ii < this->nc; ii++) {
        Col<T1> data = inData.col(ii);
        outData += this->SMap.col(ii)%((*this->G_obj)/data);

      }

      //equivalent to returning col(output) in MATLAB with IRT
      return vectorise(outData);

        
    }
};

#endif
