//
// Created by acerja2 on 10/9/15.
//

#ifndef POWERGRID_PROCESSISMRMRD_HPP
#define POWERGRID_PROCESSISMRMRD_HPP

#include <armadillo>
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/dataset.h"
#include "ismrmrd/version.h"


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

#endif //POWERGRID_PROCESSISMRMRD_HPP
