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
template<typename T1, typename Tobj>
class FieldCorrection {
    typedef complex <T1> CxT1;
public:
    FieldCorrection();

    //Class variables go here
    uword n1; //Data size
    uword n2; //Image size
    uword L; //number of time segments
    uword type; //type of time segmentation
    uword Nshots; //Number of shots, used to reduce complexity of calculating interpolator
    T1 tau;        //time segment length
    T1 T_min;   // minimum time in the time vector (i.e. TE for spiral out)
    Tobj* obj;
    Col <T1> fieldMap; //Field map (in radians per second)
    Col <T1> timeVec;  //timing vector of when each data point was collected relative to the echo time (in seconds)
    Mat <CxT1> AA;        //interpolator coefficients for the different time segments
    CxT1 i = CxT1(0., 1.);

    //Class constructor
    FieldCorrection(Tobj& G, Col <T1> map_in, Col <T1> timeVec_in, uword a, uword b, uword c, uword interptype = 1,
		    uword shots = 1)
    {
		cout << "Entering Class constructor" << endl;
	    n1 = a; //Data size
	    n2 = b;//Image size
	    L = c; //number of time segments
	    type = interptype; // type of time segmentation performed
	    Nshots = shots; // number of shots
	    obj = &G;
	    fieldMap = map_in;
    	cout << "N1 = " << n1 << endl;
		cout << "N2 = " << n2 << endl;
		cout << "L = " << L << endl;


		AA.set_size(n1, L); //time segments weights
	    timeVec = timeVec_in;
	    T_min =timeVec.min();
	    T1 rangt = timeVec.max()-T_min;
		tau = (rangt + datum::eps) / (L - 1); // it was L-1 before
	    timeVec = timeVec-T_min;

	    uword NOneShot = n1/Nshots;
	    if (L==1) {
		    tau = 0;
		    AA.ones();
	    }
	    else {
			Mat <CxT1> tempAA(NOneShot, L);
		    if (type==1) {// Hanning interpolator
			    cout << "Hanning interpolation" << endl;
			    //tau = (rangt+datum::eps)/(L-1);
			    for (unsigned int ii = 0; ii<L; ii++) {
				    for (unsigned int jj = 0; jj<NOneShot; jj++) {
					    if ((std::abs(timeVec(jj)-((ii)*tau)))<=tau) {
						    tempAA(jj, ii) = 0.5+0.5*std::cos((datum::pi)*(timeVec(jj)-((ii)*tau))/tau);
					    }
					    else {
						    tempAA(jj, ii) = 0.0;
					    }
				    }
			    }
			    AA = repmat(tempAA, Nshots, 1);
		    }
		    else if (type==2) { // Min-max interpolator: Exact LS interpolator

			    cout << "Min Max time segmentation" << endl;

			    Mat <CxT1> Ltp;
			    Ltp.ones(1, L);
			    Col <CxT1> ggtp;
			    ggtp.ones(n2, 1);
			    Mat <CxT1> gg;
			    gg = exp(i*fieldMap*tau)*Ltp;
			    Mat <CxT1> iGTGGT;
			    iGTGGT.set_size(L+1, n2);
			    Mat <CxT1> gl;
			    gl.zeros(n2, L);


			    for (unsigned int ii = 0; ii<L; ii++) {
				    for (unsigned int jj = 0; jj<n2; jj++) {
					    gl(jj, ii) = pow(gg(jj, ii), (T1) (ii+1));
				    }
			    }

			    Mat <CxT1> G;
			    G.set_size(n2, L);

			    for (unsigned int jj = 0; jj<L; jj++) {
				    if (jj==0) {
					    G.col(jj) = ggtp;
				    }
				    else {
					    G.col(jj) = gl.col(jj-1);
				    }
			    }

			    //G = joint_rows(ggtp, gl);
			    Col <CxT1> glsum;
			    Mat <CxT1> GTG;
			    GTG.zeros(L, L);
			    GTG.diag(0) += n2;
			    glsum = sum(gl.t(), 1);
				Mat <CxT1> GTGtp(L, L);
				for (unsigned int ii = 0; ii < (L - 1); ii++) {
					GTGtp.zeros();
				    GTGtp.diag(-(T1) (ii+1)) += glsum(ii);
				    GTGtp.diag((T1) (ii+1)) += std::conj(glsum(ii));
				    GTG = GTG+GTGtp;
			    }

			    //savemat("/vagrant/GTG.mat","GTGc",GTG);

			    T1 rcn = 1/cond(GTG);
			    if (rcn>10*2e-16) { //condition number of GTG
				    iGTGGT = inv(GTG)*G.t();

			    }
			    else {
				    iGTGGT = pinv(GTG)*G.t(); // pseudo inverse
			    }


			    Mat <CxT1> iGTGGTtp;
			    Mat <CxT1> ftp;
			    Col <CxT1> res, temp;

			    //cout << "Start filling the AA matrix" << endl;
			    /*
				for (unsigned int jj = 0; jj < L+1; jj++) {
				 for (unsigned int ii = 0; ii < NOneShot; ii++) {
				   iGTGGTtp =iGTGGT.row(jj);
				   ftp = exp(i*fieldMap*timeVec(ii));
				   res = as_scalar(iGTGGTtp*ftp);
				   AA(ii,jj) = conj(res);
				 }
				  cout << "Loop L" << jj << endl;
				}
				*/


			    for (unsigned int ii = 0; ii<NOneShot; ii++) {
				    ftp = exp(i*fieldMap*timeVec(ii));
				    res = iGTGGT*ftp;
				    tempAA.row(ii) = res.t();
			    }
			    AA = repmat(tempAA, Nshots, 1);
		    }
	    }
		cout << "Exiting class constructor." << endl;
    }




//Overloaded operators go here

//Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it directly rather return another vector of type T1
    Col <CxT1> operator*(const Col <CxT1>& d) const
    {
	    Tobj* G = this->obj;
	    //output is the size of the kspace data
	    Col <CxT1> outData = zeros<Col<CxT1 >> (this->n1);
		//cout << "OutData size = " << this->n1 << endl;
	    Col <CxT1> Wo;

	    //loop through time segments
	    for (unsigned int ii = 0; ii<this->L; ii++) {
		    //cout << "Entering time segmentation loop" << endl;
		    //apply a phase to each time segment
		    Wo = exp(-i*(this->fieldMap)*((ii)*this->tau+this->T_min));

		    //perform multiplication by the object and sum up the time segments
		    outData += (this->AA.col(ii))%(*G*(Wo%d));


	    }
	    return outData;
    }

    Col <CxT1> operator/(const Col <CxT1>& d) const
    {

	    //output is the size of the image
	    Col <CxT1> outData = zeros<Col<CxT1 >> (this->n2);
	    Col <CxT1> Wo;
	    //loop through the time segemtns
	    for (unsigned int ii = 0; ii<this->L; ii++) {

		    //create the phase map for the Lth time segment
		    Wo = exp(i*(this->fieldMap)*((ii)*this->tau+this->T_min));

		    //perform adjoint operation by the object and sum up the time segments
		    outData += Wo%((*this->obj)/(AA.col(ii)%d));

	    }

	    return outData;

    }
/*
protected:
  Col<T1> int_tim_seg(uword this->L) {

  }
*/


};

#endif
