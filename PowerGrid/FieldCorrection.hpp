//
//  FieldCorrection.hpp
//  PowerGrid
//
//  Created by Joe Holtrop on 4/10/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//
// This using field correction by time segmentation
// The data is corrected to time 0 with reference to the time vector passed
//
//

#ifndef PowerGrid_FieldCorrection_hpp
#define PowerGrid_FieldCorrection_hpp
// We are using two template types at the moment. One for the type of data to be processed (ie Col<cx_double>) and one for the type of G object (ie Gfft<Col<cx_double>>
// T1 is the data type for complex, T2 is the data type for real data
template <typename T1,typename T2, typename Tobj>
class FieldCorrection {
public:
  FieldCorrection();
    
    //Class variables go here
    uword n1 = 0; //Data size
    uword n2 = 0; //Image size
    uword L = 0; //number of time segments
    T2 tau;		//time segment length
    T2 T_min;   // minimum time in the time vector (i.e. TE for spiral out)
    Tobj *obj;
    Col<T2> fieldMap; //Field map (in radians per second)
    Col<T2> timeVec;  //timing vector of when each data point was collected relative to the echo time (in seconds)
    Mat<T2> AA;		//interpolator coefficients for the different time segments
    T1 i = T1(0.,1.);

    //Class constructor
    FieldCorrection(Tobj &G, Col<T2> map_in, Col<T2> timeVec_in, uword a, uword b,uword c ) {

      n1 = a; //Data size
      n2 = b;//Image size
      L = c; //number of time segments
      obj = &G;
      fieldMap = map_in;
      
      AA.set_size(n1,L); //time segments weights
      timeVec = timeVec_in;
      T_min = timeVec.min();
      T2 rangt = timeVec.max()-T_min;
      tau = (rangt+datum::eps)/(L-1);
      timeVec = timeVec - T_min;

      //Hanning interpolator
      if (L > 1){
    	  tau = (rangt+datum::eps)/(L-1);
		  for (unsigned int ii=0; ii < L; ii++) {
			for (unsigned int jj=0; jj < n1; jj++) {
			  if ((fabs(timeVec(jj)-((ii)*tau)))<=tau){
				AA(jj,ii) = 0.5 + 0.5*std::cos((datum::pi)*(timeVec(jj)-((ii)*tau))/tau);
			  } else {
                AA(jj,ii) = 0.0;
              }
			}
		  }
      }
      else { //no time segmentation needed
    	  tau = 0;
    	  AA.ones();
      }


    }
    
    //Overloaded operators go here
    
    //Forward transformation is *
    // d is the vector of data of type T1, note it is const, so we don't modify it directly rather return another vector of type T1
    Col<T1> operator*(const Col<T1>& d) const {

    //output is the size of the kspace data
      Col<T1> outData = zeros<Mat<T1>>(this->n1);
      Col<T1> Wo;
      
      //loop through time segments
      for (unsigned int ii=0; ii < this->L; ii++) {

    	//apply a phase to each time segment
		Wo = exp(-i*(this->fieldMap)*((ii)*tau+T_min));

		//perform multiplication by the object and sum up the time segments
		outData += (this->AA.col(ii))%((*this->obj)*(Wo%d));


      }

      return outData;
    }
    


    Col<T1> operator/(const Col<T1>& d) const {

       //output is the size of the image
      Col<T1> outData = zeros<Mat<T1>>(this->n2);
      Col<T1> Wo;
      //loop through the time segemtns
      for (unsigned int ii=0; ii < this->L; ii++) {

    	//create the phase map for the Lth time segment
        Wo = exp(i*(this->fieldMap)*((ii)*tau+T_min));

        //perform adjoint operation by the object and sum up the time segments
        outData +=  Wo%((*this->obj)/(AA.col(ii)%d));

      }

      return outData;

    }
};

#endif
