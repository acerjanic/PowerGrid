//
//  DWI_CGMC.hpp
//  PowerGrid
//
//  Created by Joe Holtrop on 4/22/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_DWICGMC_hpp
#define PowerGrid_DWICGMC_hpp

 template<typename T1>
class DWICGMC{
typedef complex<T1> CxT1;
public:
    DWICGMC();
    
    //Class variables go here
    uword Nd = 0; //Data size  (the size of one gdft or ggrid object, length of a single shot)
    uword Ni = 0; //Image size
    uword Nc = 0; //number of coils
    uword Ns = 0; //number of shots
    Mat<CxT1> SMap; //coil sensitivity, dimensions Image size(n1) by number of coils (nc)
    Mat<T1> PMap; //shot phase, dimensions Image size(n1) by number of shots. in radians.
    Col<T1> FMap; //Fieldmap
    Mat<T1> Kx; //kspace coordinates in x direction
    Mat<T1> Ky; //kspace coordinates in y direction
    Mat<T1> Kz; //kspace coordinates in z direction
    Col<T1> Tvec; //timing vector for a single shot (all shots assumed to have same timing vector)
    uword Nx;
    uword Ny;
    uword Nz;
    Col<T1> Ix;
    Col<T1> Iy;
    Col<T1> Iz;
    CxT1 i = CxT1(0.,1.);


    
    //Class constructor
    DWICGMC( Col<T1> kx, Col<T1> ky,  Col<T1> kz,  uword nx, uword ny, uword nz,uword nc, Col<T1> t, Col<CxT1> SENSEmap,  Col<T1> FieldMap, Col<T1> ShotPhaseMap ) {
        Ni = nx*ny*nz;
        Nc = nc;
        Ns = ShotPhaseMap.n_elem/Ni;
        Nd = kx.n_elem/Ns;
        cout <<"Nd = "<< Nd << endl;
        cout << "Ns = " << Ns << endl;
        cout <<"Nc = "<< Nc << endl;
        cout << "Ni = " << Ni << endl;
        SMap = reshape(SENSEmap,Ni,Nc);
        PMap = reshape(ShotPhaseMap,Ni,Ns);
        FMap = FieldMap;
        Kx = reshape(kx,Nd,Ns);
        Ky = reshape(ky,Nd,Ns);
        Kz = reshape(kz,Nd,Ns);
        Nx = nx;
        Ny = ny;
        Nz = nz;
        Tvec = t;

        Cube<T1> ix;
        ix.zeros(Nx,Ny,Nz);
        Cube<T1> iy;
        iy.zeros(Nx,Ny,Nz);
        Cube<T1> iz;
        iz.zeros(Nx,Ny,Nz);

        //generate the image space cordinates of the voxels we want to reconstruct
        // after vectorizing ix and iy the image coordinates must match the Field and SENSe map image coordinates
        for(uword ii = 0; ii < Ny; ii++) { //y
            for (uword jj = 0; jj < Nx; jj++) { //x
                for (uword kk = 0; kk < Nz; kk++) { //z
                    ix(ii, jj, kk) = ((T2) jj - (T2) Nx / 2.0) / ((T2) Nx);
                    iy(ii, jj, kk) = ((T2) ii - (T2) Ny / 2.0) / ((T2) Ny);
                    iz(ii, jj, kk) = ((T2) kk - (T2) Nz / 2.0) / ((T2) Nz);
                }
            }
        }

        Ix = vectorise(ix);
        Iy = vectorise(iy);
        Iz = vectorise(iz);

    }
    
    //Overloaded operators go here
    
    //Forward transformation is *
    // d is the vector of data of type T1, note it is const, so we don't modify it directly rather return another vector of type T1
    Col<CxT1> operator*(const Col<CxT1>& d) const {

        Mat<CxT1> outData = zeros<Mat<T1>>(Nd,Ns*Nc);
        cout <<"Nd = " << Nd << endl;
        cout <<"Ns = " << Ns << endl;
        cout <<"Nc = " << Nc << endl;
        cout <<"Ns = " << Ns << endl;
        //savemat("/shared/mrfil-data/data/PowerGridTest/DWI_2mm/tmp_SMap.mat","img",SMap);
        //savemat("/shared/mrfil-data/data/PowerGridTest/DWI_2mm/tmp_PMap.mat","img",PMap);
      /*  Col<T1> tmp1;
        Col<T1> tmp2;
        Col<T1> tmp3;
        Col<T1> tmp4;
        Col<T2> tmp5;
        Col<T2> tmp6;
        Col<T2> tmp7;
        tmp4.zeros(Nd);*/
        for (unsigned int jj=0; jj < Ns; jj++) {
			//Use grid or DFT?
            Gdft<T1> G(Nd, Ni, Kx.col(jj), Ky.col(jj), Kz.col(jj), Ix, Iy, Iz,FMap, Tvec);
            //Ggrid<T1,T2> G(Nd,2.0,Nx,Ny,Nz,Kx.col(jj),Ky.col(jj),Kz.col(jj),Ix,Iy,Iz);
            //In SENSE we store coil data using the columns of the data matrix, and we weight the data by the coil sensitivies from the SENSE map
          for (unsigned int ii=0; ii < Nc; ii++) {
             /*   tmp1 = exp(-i*(PMap.col(jj)));
                tmp2 = SMap.col(ii)%tmp1;
                tmp3 = d%tmp2;
                tmp4 = G*tmp3;*/
              outData.col(jj+ii*Ns) = G * (d % (SMap.col(ii)%exp(-i*(PMap.col(jj)))));

          }
      }
      //equivalent to returning col(output) in MATLAB with IRT
      return vectorise(outData);
    }
    //For the adjoint operation, we have to weight the adjoint transform of the coil data by the SENSE map.
    Col<CxT1> operator/(const Col<CxT1>& d) const {

      Mat<CxT1> inData = reshape(d,Nd,Ns*Nc);

      Col<CxT1> outData = zeros<Col<CxT1>>(Ni);

          for (unsigned int jj=0; jj < Ns; jj++) {
		
			//Use grid or DFT?
          Ggrid<T1> G(Nd,2.0,Nx,Ny,Nz,Kx.col(jj),Ky.col(jj),Kz.col(jj),Ix,Iy,Iz);


              for (unsigned int ii=0; ii < Nc; ii++) {

              outData += conj(SMap.col(ii)%exp(-i*(PMap.col(jj))))%((G)/inData.col(jj+ii*Ns));

          }

         // temp = ((*this->G_obj)/data);
         // savemat("/shared/mrfil-data/data/PowerGridTest/64_64_16_4coils/coil_img_" + std::to_string(ii) + ".mat","img",temp);

      }

      //equivalent to returning col(output) in MATLAB with IRT
      return vectorise(outData);

        
    }
};

#endif
