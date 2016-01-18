//
//  mpiDWICGMCDFT.hpp
//  PowerGrid
//
//  Created by Joe Holtrop on 4/22/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_mpiDWICGMCDFT_hpp
#define PowerGrid_mpiDWICGMCDFT_hpp
using namespace arma;

namespace bmpi = boost::mpi;
template<typename T1>
class mpiDWICGMCDFT {
    typedef std::complex <T1> CxT1;
public:
    mpiDWICGMCDFT();

    ~mpiDWICGMCDFT() {
        /*
        for (uword jj = 0; jj < Ns; jj++) {
            delete AObj[jj];
            //delete G[jj];
        }
        delete[] AObj;
        //delete[] G;
        */
        delete AObj;
    }

    //Class variables go here
    uword Nd = 0; //Data size  (the size of one gdft or ggrid object, length of a single shot)
    uword Ni = 0; //Image size
    uword Nc = 0; //number of coils
    uword Ns = 0; //number of shots
    Mat <CxT1> SMap; //coil sensitivity, dimensions Image size(n1) by number of coils (nc)
    Mat <T1> PMap; //shot phase, dimensions Image size(n1) by number of shots. in radians.
    Col <T1> FMap; //Fieldmap
    Mat <T1> Kx; //kspace coordinates in x direction
    Mat <T1> Ky; //kspace coordinates in y direction
    Mat <T1> Kz; //kspace coordinates in z direction
    Mat <T1> Tvec; //timing vector for a single shot (all shots assumed to have same timing vector)
    uword Nx;
    uword Ny;
    uword Nz;
    Col <T1> Ix;
    Col <T1> Iy;
    Col <T1> Iz;
    CxT1 i = CxT1(0., 1.);
    uword type = 1; // 2 for min max time seg and 1 for Hanning
    uword L = 20;
    //Gdft <T1> **AObj = NULL;
    Gdft <T1> *AObj = NULL;
    //FieldCorrection <T1, Ggrid<T1>> **AObj = NULL;
    bmpi::environment *env;
    bmpi::communicator *world;
    //Class constructor
    mpiDWICGMCDFT(Col <T1> kx, Col <T1> ky, Col <T1> kz, uword nx, uword ny, uword nz, uword nc, Col <T1> t,
                  Col <CxT1> SENSEmap, Col <T1> FieldMap, Col <T1> ShotPhaseMap, bmpi::environment &en,
                  bmpi::communicator &wor) {
        Ni = nx * ny * nz;
        Nc = nc;
        Ns = ShotPhaseMap.n_elem / Ni;
        Nd = kx.n_elem / Ns;
        std::cout << "Nd = " << Nd << std::endl;
        std::cout << "Ns = " << Ns << std::endl;
        std::cout << "Nc = " << Nc << std::endl;
        std::cout << "Ni = " << Ni << std::endl;
        SMap = reshape(SENSEmap, Ni, Nc);
        PMap = reshape(ShotPhaseMap, Ni, Ns);
        FMap = FieldMap;
        Kx = reshape(kx, Nd, Ns);
        Ky = reshape(ky, Nd, Ns);
        Kz = reshape(kz, Nd, Ns);
        Nx = nx;
        Ny = ny;
        Nz = nz;
        Tvec = reshape(t, Nd, Ns);
        env = &en;
        world = &wor;

        Cube <T1> ix;
        ix.zeros(Nx, Ny, Nz);
        Cube <T1> iy;
        iy.zeros(Nx, Ny, Nz);
        Cube <T1> iz;
        iz.zeros(Nx, Ny, Nz);

        //generate the image space coordinates of the voxels we want to reconstruct
        // after vectorizing ix and iy the image coordinates must match the Field and SENSe map image coordinates
        for (uword ii = 0; ii < Ny; ii++) { //y
            for (uword jj = 0; jj < Nx; jj++) { //x
                for (uword kk = 0; kk < Nz; kk++) { //z
                    ix(ii, jj, kk) = ((T1) jj - (T1) Nx / 2.0) / ((T1) Nx);
                    iy(ii, jj, kk) = ((T1) ii - (T1) Ny / 2.0) / ((T1) Ny);
                    iz(ii, jj, kk) = ((T1) kk - (T1) Nz / 2.0) / ((T1) Nz);
                }
            }
        }

        Ix = vectorise(ix);
        Iy = vectorise(iy);
        Iz = vectorise(iz);

        //AObj = new Gdft <T1> *[Ns];
        //AObj = new FieldCorrection <T1, Ggrid<T1>> *[Ns];
        //AObj = new Gdft <T1>;
        // Initialize the field correction and G objects we need for this reconstruction
        /*for (uword jj = 0; jj < Ns; jj++) {

            AObj[jj] = new Gdft<T1>(Nd, Nx * Ny * Nz, Kx.col(jj), Ky.col(jj), Kz.col(jj), Ix, Iy, Iz, vectorise(FMap),
                                    vectorise(Tvec.col(jj)));
            //AObj[jj] = new FieldCorrection <T1, Ggrid<T1>>(*G[jj], vectorise(FMap), vectorise(Tvec.col(jj)),
            //                                               (uword) Nd, (uword)(Nx * Ny * Nz), (uword) L, type,
            //                                               (uword) 1);

        } */
        //This gives us the rank (index) of the process in the collective MPI space.
        //The process will only calculate shot jj of the acquisition.
        uword ShotRank = world->rank();
        uword numShots = Ns / world->size();
        AObj = new Gdft<T1>(Nd, Nx * Ny * Nz, Kx.col(ShotRank), Ky.col(ShotRank), Kz.col(ShotRank), Ix, Iy, Iz,
                            vectorise(FMap), vectorise(Tvec.col(ShotRank)));

    }

    //Overloaded operators go here

    //Forward transformation is *
    // d is the vector of data of type T1, note it is const, so we don't modify it directly rather return another vector of type T1
    Col <CxT1> operator*(const Col <CxT1> &d) const {
        uword ShotRank = world->rank();
        Mat <CxT1> outData = zeros < Mat < CxT1 >> (Nd, Ns * Nc);
        Mat <CxT1> tempOutData = zeros < Mat < CxT1 >> (Nd, Nc);
        Mat <CxT1> tempOutData2 = zeros < Mat < CxT1 >> (Nd, Nc);

        //Shot loop. Each shot has it's own kspace trajectory
        //for (unsigned int jj = 0; jj < Ns; jj++) {

            //Coil loop. Each coil exists for each shot, so we need to work with these.
            for (uword ii = 0; ii < Nc; ii++) {
                tempOutData.col(ii) = (*AObj) * (d % (SMap.col(ii) % exp(-i * (PMap.col(ShotRank)))));
            }
        std::cout << "Performed forward xforms." << std::endl;
        //}
        //Now let's do some MPI stuff here.
        if (world->rank() == 0) {
            std::vector <Mat<CxT1>> OutDataGather(world->size());
            // Collect all the data into OutDataGather an std::vector collective
            std::cout << "Rank #: " << world->rank() << " reached foward xform gather" << std::endl;
            bmpi::gather < Mat < CxT1 >> (*world, tempOutData, OutDataGather, 0);
            std::cout << "Rank #: " << world->rank() << " passed forward xform gather" << std::endl;
            for(uword jj = 0; jj < Ns; jj++) {
                std::cout << "About to access element #" << jj << " in gathered data" << std::endl;
                std::cout << "Size of returned stl::vector<> = " << OutDataGather.size() << std::endl;
                tempOutData2 = OutDataGather[jj];
                for (uword ii = 0; ii < Nc; ii++) {
                    outData.col(ii + jj * Nc) = tempOutData2.col(ii);
                }
            }

        } else {
            std::cout << "Rank #: " << world->rank() << " reached foward xform gather" << std::endl;
            bmpi::gather < Mat < CxT1 >> (*world, tempOutData, 0);
        }
        std::cout << "Rank #: " << world->rank() << " reached foward xform broadcast" << std::endl;
        bmpi::broadcast < Mat < CxT1 >> (*world, outData, 0);
        std::cout << "Rank #: " << world->rank() << " passed foward xform broadcast" << std::endl;

        //equivalent to returning col(output) in MATLAB with IRT
        return vectorise(outData);
    }

    //For the adjoint operation, we have to weight the adjoint transform of the coil data by the SENSE map.
    Col <CxT1> operator/(const Col <CxT1> &d) const {
        uword ShotRank = world->rank();
        Mat <CxT1> inData = reshape(d, Nd, Ns * Nc);

        //Col <CxT1> tempOutData = zeros < Col < CxT1 >> (Ni);
        Col <CxT1> outData = zeros < Col < CxT1 >> (Ni);
        //Shot Loop. Each shot has it's own k-space trajectory
        //for (unsigned int jj = 0; jj < Ns; jj++) {

            //Coil Loop - for each shot we have a full set of coil data.
            for (unsigned int ii = 0; ii < Nc; ii++) {
                outData += conj(SMap.col(ii) % exp(-i * (PMap.col(ShotRank)))) %
                           ((*AObj) / inData.col(ShotRank + ii * Ns));
                //std::cout << "Processed shot # " << jj << " coil # " << ii << std::endl;
            }

        //}
        if (world->rank() == 0) {
            std::vector<Col<CxT1>> OutDataGather;
            // Collect all the data into OutDataGather an std::vector collective
            bmpi::gather < Col < CxT1 >> (*world, outData, OutDataGather, 0);
            for(uword jj = 1; jj < Ns; jj++) {
                //Skip the zeroth rank because that is the outData we started with on this processor.
                //Now we will manually reduce the data by adding across the elements.
                outData += OutDataGather.at(jj);
            }

        } else {
            bmpi::gather < Col < CxT1 >> (*world, outData, 0);

        }
        //Broadcast will send the data to all machines if the rank==0 and receive the broadcasted data from rank=0 to
        // all other machines in the communicator, overwriting the outData value already there.
        bmpi::broadcast < Col < CxT1 >> (*world, outData, 0);
        //equivalent to returning col(output) in MATLAB with IRT
        return vectorise(outData);

    }


};

#endif //PowerGrid_mpiDWICGMCDFT_hpp
