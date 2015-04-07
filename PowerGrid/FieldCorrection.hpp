//
//  FiledCorrection.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/4/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//
// This using field correction by time segmentation
// The data is correctied to time 0 with reference to the time vector passed
//
//

#ifndef PowerGrid_FiledCorrection_hpp
#define PowerGrid_FiledCorrection_hpp
// We are using two template types at the moment. One for the type of data to be processed (ie Col<cx_double>) and one for the type of G object (ie Gfft<Col<cx_double>>
// T1 is the data type for complex, T2 is the data type for real data
template <typename T1,typename T2, typename Tobj>
class FiledCorrection {
public:
  FiledCorrection();
    
    //Class variables go here
    uword n1 = 0; //Data size
    uword n2 = 0; //Image size
    uword L = 0; //number of time segments
    T2 tau;
    Tobj *obj;
    T2 fieldMap; //filed map in radians per second
    Col<T2> timeVec;
    Mat<T2> AA;

    
    //Class constructor
    FiledCorrection(Tobj &G, Col<T2> map_in, Col<T2> timeVec_in, uword a, uword b,uword c ) {
      n1 = a; //Data size
      n2 = b;//Image size
      L = c; //number of time segments
      obj = &G;
      filedMap = map_in;
      timeVec = timeVec_in;
      AA = zeros<mat>(n1,L); //time segments weights
      T2 rangt = timeVec.max()-timeVec.min();
      tau = (rangt+datum::eps)/(L-1);
      //Hanning interpolator
      for (unsigned int ii=0; ii < L; ii++) {
        for (unsigned int jj=0; jj < n1; jj++) {
          if (abs(timeVec(jj)-((ii-1)*tau))<=tau){
            AA(jj,ii) = 0.5 + 0.5*cos((datum::pi)*(timeVec(jj)-((ii-1)*tau))/tau);
          }
        }
      }

    }
    
    //Overloaded operators go here
    
    //Forward transformation is *
    // d is the vector of data of type T1, note it is const, so we don't modify it directly rather return another vector of type T1
    Col<T1> operator*(const Col<T1>& d) const {


      Col<T1> outData = zeros<Mat<T1>>(this->n2);

      //loop through time segments
      for (unsigned int ii=0; ii < this->L; ii++) {

        Col<T1> Wo = exp(-sqrt(-1)*this->fieldMap*((ll-1)*tau));
        outData += this->AA.col(ii)%(this->obj*(Wo%d));

      }

      return outData;
    }
    
    Col<T1> operator/(const Col<T1>& d) const {

      Col<T1> outData = zeros<Mat<T1>>(this->n1);

      for (unsigned int ii=0; ii < this->L; ii++) {

        Col<T1> Wo = exp(-sqrt(-1)*this->fieldMap*((ll-1)*tau));

        outData +=  Wo%(this->obj)/(AA.col(ii)%d);

      }

      return outData;

    }
};

#endif