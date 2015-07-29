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
    uword type = 0; //type of time segmentation
    T2 tau;		//time segment length
    T2 T_min;   // minimum time in the time vector (i.e. TE for spiral out)
    Tobj *obj;
    Col<T2> fieldMap; //Field map (in radians per second)
    Col<T2> timeVec;  //timing vector of when each data point was collected relative to the echo time (in seconds)
    Mat<T1> AA;		//interpolator coefficients for the different time segments AA for T2 before but is T1 now - GCN
    T1 i = T1(0.,1.);

    FieldCorrection(Tobj &G, Col<T2> map_in, Col<T2> timeVec_in, uword a, uword b,uword c) {
      FieldCorrection(G, map_in, timeVec_in, a, b,c, 1 );

    }

    //Class constructor
    FieldCorrection(Tobj &G, Col<T2> map_in, Col<T2> timeVec_in, uword a, uword b,uword c, uword interptype ) {

      n1 = a; //Data size
      n2 = b;//Image size
      L = c; //number of time segments
      type = interptype; // type of time segmentation performed
      obj = &G;
      fieldMap = map_in;
      
      AA.set_size(n1,L+1); //time segments weights
      timeVec = timeVec_in;
      T_min = timeVec.min();
      T2 rangt = timeVec.max()-T_min;
      tau = (rangt+datum::eps)/(L); // it was L-1 before
      timeVec = timeVec - T_min;

      savemat("/vagrant/AAinit.mat","AAinit",AA);


      if (type == 1) {// Hanning interpolator
        cout << "Hanning interpolation" << endl;
        //tau = (rangt+datum::eps)/(L-1);
		  for (unsigned int ii=0; ii < L+1; ii++) {
			for (unsigned int jj=0; jj < n1; jj++) {
			  if ((fabs(timeVec(jj)-((ii)*tau)))<=tau){
              	AA(jj,ii) = 0.5 + 0.5*std::cos((datum::pi)*(timeVec(jj)-((ii)*tau))/tau);
			  } else {
                AA(jj,ii) = 0.0;
              }
			}
          }
      }

      else if (type == 2) { // Min-max interpolator: Exact LS interpolator

        cout << "Min Max time segmentation" << endl;
        //cout << "nro" << n2<< endl;
        //cout << "L" << L << endl;
        //cout << "ndat" << n1<< endl;
        Mat <T1> Ltp;
        Ltp.ones(1, L);
        Col <T1> ggtp;
        ggtp.ones(n2, 1);
        Mat <T1> gg;
        gg = exp(i * fieldMap * tau)*Ltp;
        Mat <T1> iGTGGT;
        iGTGGT.set_size(L+1,n2);
        savemat("/vagrant/gg.mat","ggc",gg);
        Mat <T1> gl;
        gl.zeros(n2, L);
        savemat("/vagrant/fieldMap.mat","fieldMap",fieldMap);
        savemat("/vagrant/timeVec.mat","timeVecc",timeVec);


        //cout << "M0" << endl;
        for (unsigned int ii = 0; ii < L; ii++) {
          for (unsigned int jj = 0; jj < n2; jj++) {
           gl(jj,ii) = pow(gg(jj,ii),double (ii+1));
           }
        }
        savemat("/vagrant/gl.mat","glc",gl);

        Mat <T1> G;
        G.set_size(n2,L+1);

        for (unsigned int jj = 0; jj<L+1; jj++){
           if (jj == 0) {
              G.col(jj) = ggtp;
            }
            else {
              G.col(jj) = gl.col(jj-1);
            }
        }
        savemat("/vagrant/G.mat","Gc",G);
        cout << "LLL" << endl;

        //G = joint_rows(ggtp, gl);
        Col <T1> glsum;
        Mat <T1> GTG;
        GTG.zeros(L + 1, L + 1);
        GTG.diag(0) += n2;
        glsum = sum(gl.t(),1);
        savemat("/vagrant/glsum.mat","glsumc",glsum);

        for (unsigned int ii = 0; ii < L ; ii++) {
          Mat <T1> GTGtp;
          GTGtp.zeros(L + 1, L + 1);
          GTGtp.diag(- double (ii+1)) += glsum(ii);
          GTGtp.diag( double (ii+1)) += std::conj(glsum(ii));
          GTG = GTG + GTGtp;
        }

        savemat("/vagrant/GTG.mat","GTGc",GTG);

        T2 rcn = 1/cond(GTG);
        if (rcn > 10*2e-16) { //condition number of GTG
         iGTGGT = inv(GTG)*G.t();
        }
        else {
         iGTGGT = pinv(GTG)*G.t(); // pseudo inverse
        }

        savemat("/vagrant/iGTGGT.mat","iGTGGTc",iGTGGT);
        savemat("/vagrant/timeVec.mat","timeVecc",timeVec);

        Mat <T1> iGTGGTtp;
        Mat <T1> ftp;
        Col<T1> res;

        cout << "Start filling the AA matrix" << endl;
        /*for (unsigned int jj = 0; jj < L+1; jj++) {
         for (unsigned int ii = 0; ii < n1; ii++) {
           iGTGGTtp =iGTGGT.row(jj);
           ftp = exp(i*fieldMap*timeVec(ii));
           res = iGTGGTtp*ftp;
           AA(ii,jj) = std::conj(res);
         }
          cout << "Loop L" << jj << endl;
        }*/



          for (unsigned int ii = 0; ii < n1; ii++) {
            ftp = exp(i*fieldMap*timeVec(ii));
            //savemat("/vagrant/ftp.mat","ftpc",ftp);
            //cout << "L" << endl;

            res = iGTGGT*ftp;
            //cout << "N" << endl;
            //savemat("/vagrant/res.mat","res",res);

            AA.row(ii) = conj(res.t());
          }



        savemat("/vagrant/AA.mat","AAc",AA);

      }


      else { //no time segmentation needed
    	  tau = 0;
    	  AA.ones();
      }

      savemat("/vagrant/AA.mat","AAc",AA);
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
        cout << "Entering time segmentation loop" << endl;
    	//apply a phase to each time segment
		Wo = exp(-i*(this->fieldMap)*((ii)*tau+T_min));

		//perform multiplication by the object and sum up the time segments
		outData += (this->AA.col(ii))%((*this->obj)*(Wo%d));


      }
      savemat("/vagrant/outForward.mat","outForward",outData);
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
      savemat("/vagrant/outBackward.mat","outBackward",outData);

      return outData;

    }
/*
protected:
  Col<T1> int_tim_seg(uword L) {

  }
*/


};

#endif
