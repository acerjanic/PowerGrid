//
// Created by acerja2 on 10/9/15.
//

#ifndef POWERGRID_PROCESSISMRMRD_HPP
#define POWERGRID_PROCESSISMRMRD_HPP

#include "PowerGrid.h"

using namespace arma;

void openISMRMRDData(std::string inputDataFile, ISMRMRD::Dataset *&d, ISMRMRD::IsmrmrdHeader &hdr) {
    std::cout << "trying to create an ISMRMD::Dataset object" << std::endl;
    d = new ISMRMRD::Dataset(inputDataFile.c_str(), "dataset", false);
    //std::cout << "address of the  ISMRMD::Dataset object = " << d << std::endl;
    std::string xml;
    std::cout << "trying to read the header from the ISMRMD::Dataset object" << std::endl;
    d->readHeader(xml);
    std::cout << "read the header from the ISMRMD::Dataset object" << std::endl;
    //ISMRMRD::IsmrmrdHeader hdr;
    ISMRMRD::deserialize(xml.c_str(), hdr);
    std::cout << "trying to deserialze the xml header from the string" << std::endl;

}


//Write conversion from Image format to Armadillo matrix format for further use.
template<typename T1>
arma::Col<T1> convertFromNDArrayToArma(ISMRMRD::NDArray<T1> &inArray) {
    std::cout << "Converting NDArray to Arma Column" << std::endl;
    arma::Col<T1> temp;
    arma::uword numElems = inArray.getNumberOfElements();
    std::cout << "Num Elements to be converted = " << numElems << std::endl;
    const T1 *auxData = NULL;
    auxData = inArray.getDataPtr();
    temp = arma::Col<T1>(auxData, numElems);

    return temp;
}

template<typename T1>
arma::Col<T1> getISMRMRDFieldMap(ISMRMRD::Dataset *d) {
    const std::string fieldMap = "FieldMap";
    arma::Col<T1> FM;
    if (d->getNumberOfNDArrays(fieldMap) > 1) {
        //Throw error here
    }
    ISMRMRD::NDArray<T1> tempArray;
    d->readNDArray(fieldMap, 0, tempArray);
    FM = convertFromNDArrayToArma(tempArray);

    return FM;
}

template<typename T1>
arma::Col<T1> getISMRMRDSenseMap(ISMRMRD::Dataset *d) {
    const std::string senseMap = "SENSEMap";
    arma::Col<T1> sen;
    if (d->getNumberOfNDArrays(senseMap) > 1) {
        //Throw error here
    }
    ISMRMRD::NDArray<T1> tempArray;
    d->readNDArray(senseMap, 0, tempArray);
    sen = convertFromNDArrayToArma(tempArray);

    return sen;
}

template<typename T1>
void processISMRMRDInput(std::string inputDataFile, ISMRMRD::Dataset *&d, ISMRMRD::IsmrmrdHeader &hdr,
                         arma::Col<T1> &FM, arma::Col<std::complex<T1>> &sen) {
    std::cout << "About to open ISMRMRD file for input" << std::endl;
    openISMRMRDData(inputDataFile, d, hdr);
    std::cout << "Opened ISMRMRD file for input " << std::endl;
    std::cout << "About to get the SENSE map" << std::endl;
    sen = getISMRMRDSenseMap<std::complex<T1>>(d);
    std::cout << "About to get the Field map" << std::endl;
    FM = getISMRMRDFieldMap<T1>(d);

    return;
}

template<typename T1>
void getISMRMRDAcqData(ISMRMRD::Dataset *d, uword Nacq, Col<std::complex<T1>> &data, Col<T1> &kx, Col<T1> &ky,
                       Col<T1> &kz, Col<T1> &tvec) {
    ISMRMRD::Acquisition acq;
    d->readAcquisition(Nacq, acq);
    uword nro = acq.number_of_samples();
    uword nc = acq.active_channels();

    Mat<std::complex<T1>> dataTemp(nro, nc);
    for (uword jj = 0; jj < nc; jj++) {
        for (uword kk = 0; kk < nro; kk++) {
            dataTemp(kk, jj) = acq.data(kk, jj);
        }
    }

    //Deal with trajectories
    for (uword ii = 0; ii < acq.number_of_samples(); ii++) {
        kx(ii) = acq.traj(0, ii);
        ky(ii) = acq.traj(1, ii);
        kz(ii) = acq.traj(2, ii);
        tvec(ii) = acq.traj(3, ii);
    }
    data = vectorise(dataTemp);
    return;
}

template<typename T1>
ISMRMRD::Acquisition getISMRMRDAcq(ISMRMRD::Dataset *d, uword Nacq) {
    ISMRMRD::Acquisition acq;
    d->readAcquisition(Nacq, acq);
    uword nro = acq.number_of_samples();
    uword nc = acq.active_channels();

    Mat<std::complex<T1>> dataTemp(nro, nc);
    for (uword jj = 0; jj < nc; jj++) {
        for (uword kk = 0; kk < nro; kk++) {
            dataTemp(kk, jj) = acq.data(kk, jj);
        }
    }

    //Deal with trajectories
    for (uword ii = 0; ii < acq.number_of_samples(); ii++) {
        kx(ii) = acq.traj(0, ii);
        ky(ii) = acq.traj(1, ii);
        kz(ii) = acq.traj(2, ii);
        tvec(ii) = acq.traj(3, ii);
    }
    data = vectorise(dataTemp);
    return acq;
}

template<typename T1>
void writeISMRMRDImageData(ISMRMRD::Dataset *d, Col<std::complex<T1>> &image, uword Nx, uword Ny, uword Nz) {
    ISMRMRD::Image<std::complex<T1>> img_out(Nx, Ny, Nz, 1);
    savemat("testImage.mat", "img", image);
    for (int ii = 0; ii < Ny; ii++) {
        for (int jj = 0; jj < Nx; jj++) {
            for (int kk = 0; kk < Nz; kk++) {
                img_out(jj, ii, kk) = image(ii + jj * Ny + kk * Ny * Nx);
            }
        }
    }

    //Let's set some header details
    img_out.setImageType(ISMRMRD::ISMRMRD_IMTYPE_COMPLEX);

    //Write out the image
    d->appendImage("image", img_out);

}

template<typename T1>
void getCompleteISMRMRDAcqData(ISMRMRD::Dataset *d, uword NSlice, uword NSet, uword NRep, Col<std::complex<T1>> &data,
                               Col<T1> &kx, Col<T1> &ky, Col<T1> &kz, Col<T1> &tvec) {

    //Initialization
    Mat<std::complex<T1>> dataTemp;
    Mat<std::complex<T1>> acqTemp;
    Col < T1 > kxTemp, kyTemp, kzTemp, tvecTemp;
    uword numAcq = d->getNumberOfAcquisitions();
    bool firstData = true;
    ISMRMRD::Acquisition acq;
    std::cout << "Num of acquisitions in dataset = " << numAcq << std::endl;

    for (uword acqIndx = 0; acqIndx < numAcq; acqIndx++) {
        std::cout << "Scanning through acquisition # " << acqIndx << std::endl;
        // Scanning through the file.
        d->readAcquisition(acqIndx, acq);
        uword nro = acq.number_of_samples();
        uword nc = acq.active_channels();
        ISMRMRD::EncodingCounters encIdx = acq.idx();

        if ((encIdx.set == NSet) && (encIdx.repetition == NRep) && (encIdx.slice == NSlice)) {
            std::cout << "Grabbing acq index #" << acqIndx << std::endl;
            acqTemp.zeros(nro, nc);
            kxTemp.set_size(nro);
            kyTemp.set_size(nro);
            kzTemp.set_size(nro);
            tvecTemp.set_size(nro);

            for (uword jj = 0; jj < nc; jj++) {
                for (uword kk = 0; kk < nro; kk++) {
                    acqTemp(kk, jj) = static_cast<std::complex<T1>>(acq.data(kk, jj));
                }
            }
            std::cout << "Get Number of Trajectoriy entries = " << acq.getNumberOfTrajElements() << std::endl;

            //Deal with trajectories
            for (uword ii = 0; ii < nro; ii++) {
                kxTemp(ii) = static_cast<T1>(acq.traj(0, ii));
                kyTemp(ii) = static_cast<T1>(acq.traj(1, ii));
                kzTemp(ii) = static_cast<T1>(acq.traj(2, ii));
                tvecTemp(ii) = static_cast<T1>(acq.traj(3, ii));
            }

            //Append Data points to vectors
            if (firstData) {
                firstData = false; //Sentinel
                dataTemp = acqTemp;
                kx = kxTemp;
                ky = kyTemp;
                kz = kzTemp;
            } else {
                dataTemp = join_vert(dataTemp, acqTemp);
                kx = join_vert(kx, kxTemp);
                ky = join_vert(ky, kyTemp);
                kz = join_vert(kz, kzTemp);
                tvec = join_vert(tvec, tvecTemp);
            }

        }
    }

    //Vectorise coils from matrix to column vector
    data = vectorise(dataTemp);
    savemat("kxOut.mat", "kx", kx);
    savemat("kyOut.mat", "ky", ky);
    savemat("kzOut.mat", "kz", kz);
    savemat("tvecOut.mat", "tvec", tvec);
    savemat("dataOut.mat", "dataOut", data);

    return;
}

#endif //POWERGRID_PROCESSISMRMRD_HPP
