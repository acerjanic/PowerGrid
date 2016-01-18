//
// Created by Alex Cerjanic on 1/10/16.
//

#ifndef POWERGRID_ARMASERIALIZATION_HPP
#define POWERGRID_ARMASERIALIZATION_HPP

#include <boost/serialization/complex.hpp>
/*
BOOST_SERIALIZATION_SPLIT_FREE(arma::Col<double>)
BOOST_SERIALIZATION_SPLIT_FREE(arma::Mat<double>)
BOOST_SERIALIZATION_SPLIT_FREE(arma::Col<float>)
BOOST_SERIALIZATION_SPLIT_FREE(arma::Mat<float>)
BOOST_SERIALIZATION_SPLIT_FREE(arma::Col<std::complex<double>>)
BOOST_SERIALIZATION_SPLIT_FREE(arma::Mat<std::complex<double>>)
BOOST_SERIALIZATION_SPLIT_FREE(arma::Col<std::complex<float>>)
BOOST_SERIALIZATION_SPLIT_FREE(arma::Mat<std::complex<float>>)
*/
namespace boost {
    namespace serialization {

        //Serialization code for arma::Col<double>

        template<class Archive>
        void serialize(Archive &ar, arma::Col<double> &g, const unsigned int version)
        {   //Defer the serialization to seperate load and save functions
            boost::serialization::split_free(ar, g, version);
        }

        template<class Archive>
        void save(Archive &ar, const arma::Col<double> &g, const unsigned int version)
        {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            //Start with the length of the vector
            ar & g.n_rows;
            //Now read out the data. This is probably a slow way to do it.
            for (uword ii = 0; ii < g.n_elem; ii++) {
                ar << g(ii);
            }
        }

        template<class Archive>
        void load(Archive &ar, arma::Col<double> &g, const unsigned int version)
        {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            unsigned int nRows = 0;
            double tempT;
            //Start with the length of the vector
            ar & nRows;
            //Make sure the object is set to the right size.
            g.set_size(nRows);
            //Now read the elements out in the order they were read in. Probably a slow way to do this.
            for (uword ii = 0; ii < g.n_elem; ii++) {
                ar & tempT;
                g(ii) = tempT;
            }
        }

        //Serialization code for arma::Col<float>

        template<class Archive>
        void serialize(Archive &ar, arma::Col<float> &g,
                       const unsigned int version) {   //Defer the serialization to seperate load and save functions
            boost::serialization::split_free(ar, g, version);
        }

        template<class Archive>
        void save(Archive &ar, const arma::Col<float> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            //Start with the length of the vector
            ar & g.n_rows;
            //Now read out the data. This is probably a slow way to do it.
            for (uword ii = 0; ii < g.n_elem; ii++) {
                ar << g(ii);
            }
        }

        template<class Archive>
        void load(Archive &ar, arma::Col<float> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            unsigned int nRows = 0;
            float tempT;
            //Start with the length of the vector
            ar & nRows;
            //Make sure the object is set to the right size.
            g.set_size(nRows);
            //Now read the elements out in the order they were read in. Probably a slow way to do this.
            for (uword ii = 0; ii < g.n_elem; ii++) {
                ar & tempT;
                g(ii) = tempT;
            }
        }

        //Serialization code for arma::Col<std::complex<float>>

        template<class Archive>
        void serialize(Archive &ar, arma::Col<std::complex<float>> &g,
                       const unsigned int version) {   //Defer the serialization to seperate load and save functions
            boost::serialization::split_free(ar, g, version);
        }

        template<class Archive>
        void save(Archive &ar, const arma::Col<std::complex<float>> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            //Start with the length of the vector
            ar & g.n_rows;
            //Now read out the data. This is probably a slow way to do it.
            for (uword ii = 0; ii < g.n_elem; ii++) {
                ar & g(ii);
            }
        }

        template<class Archive>
        void load(Archive &ar, arma::Col<std::complex<float>> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            unsigned int nRows = 0;
            std::complex<float> tempT;
            //Start with the length of the vector
            ar & nRows;
            //Make sure the object is set to the right size.
            g.set_size(nRows);
            //Now read the elements out in the order they were read in. Probably a slow way to do this.
            for (uword ii = 0; ii < g.n_elem; ii++) {
                ar & tempT;
                g(ii) = tempT;
            }
        }

        //Serialization code for arma::Col<std::complex<double>>

        template<class Archive>
        void serialize(Archive &ar, arma::Col<std::complex<double>> &g,
                       const unsigned int version) {   //Defer the serialization to seperate load and save functions
            boost::serialization::split_free(ar, g, version);
        }

        template<class Archive>
        void save(Archive &ar, const arma::Col<std::complex<double>> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            //Start with the length of the vector
            ar & g.n_rows;
            //Now read out the data. This is probably a slow way to do it.
            for (uword ii = 0; ii < g.n_elem; ii++) {
                ar & g(ii);
            }
        }

        template<class Archive>
        void load(Archive &ar, arma::Col<std::complex<double>> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            unsigned int nRows = 0;
            std::complex<double> tempT;
            //Start with the length of the vector
            ar & nRows;
            //Make sure the object is set to the right size.
            g.set_size(nRows);
            //Now read the elements out in the order they were read in. Probably a slow way to do this.
            for (uword ii = 0; ii < g.n_elem; ii++) {
                ar & tempT;
                g(ii) = tempT;
            }
        }

        //Serialization code for arma::Mat<double>

        template<class Archive>
        void serialize(Archive &ar, arma::Mat<double> &g,
                       const unsigned int version) {   //Defer the serialization to seperate load and save functions
            boost::serialization::split_free(ar, g, version);
        }

        template<class Archive>
        void save(Archive &ar, const arma::Mat<double> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            //Start with the length of the vector
            ar & g.n_rows;
            ar & g.n_cols;
            //Now read out the data. This is probably a slow way to do it.
            for (uword ii = 0; ii < g.n_rows; ii++) {
                for (uword jj = 0; jj < g.n_cols; jj++) {
                    ar << g(ii, jj);
                }
            }
        }

        template<class Archive>
        void load(Archive &ar, arma::Mat<double> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            unsigned int nRows = 0;
            unsigned int nCols = 0;
            double tempT;
            //Start with the length of the vector
            ar & nRows;
            ar & nCols;
            //Make sure the object is set to the right size.
            g.set_size(nRows, nCols);
            //Now read the elements out in the order they were read in. Probably a slow way to do this.
            for (uword ii = 0; ii < nRows; ii++) {
                for (uword jj = 0; jj < nCols; jj++) {
                    ar & tempT;
                    g(ii, jj) = tempT;
                }
            }
        }

        //Serialization code for arma::Mat<float>

        template<class Archive>
        void serialize(Archive &ar, arma::Mat<float> &g,
                       const unsigned int version) {   //Defer the serialization to seperate load and save functions
            boost::serialization::split_free(ar, g, version);
        }

        template<class Archive>
        void save(Archive &ar, const arma::Mat<float> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            //Start with the length of the vector
            ar & g.n_rows;
            ar & g.n_cols;
            //Now read out the data. This is probably a slow way to do it.
            for (uword ii = 0; ii < g.n_rows; ii++) {
                for (uword jj = 0; jj < g.n_cols; jj++) {
                    ar & g(ii, jj);
                }
            }
        }

        template<class Archive>
        void load(Archive &ar, arma::Mat<float> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            unsigned int nRows = 0;
            unsigned int nCols = 0;
            float tempT;
            //Start with the length of the vector
            ar & nRows;
            ar & nCols;
            //Make sure the object is set to the right size.
            g.set_size(nRows, nCols);
            //Now read the elements out in the order they were read in. Probably a slow way to do this.
            for (uword ii = 0; ii < nRows; ii++) {
                for (uword jj = 0; jj < nCols; jj++) {
                    ar & tempT;
                    g(ii, jj) = tempT;
                }
            }
        }

        //Serialization code for arma::Mat<std::complex<double>>

        template<class Archive>
        void serialize(Archive &ar, arma::Mat<std::complex<double>> &g,
                       const unsigned int version) {   //Defer the serialization to seperate load and save functions
            boost::serialization::split_free(ar, g, version);
        }

        template<class Archive>
        void save(Archive &ar, const arma::Mat<std::complex<double>> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            //Start with the length of the vector
            ar & g.n_rows;
            ar & g.n_cols;
            //Now read out the data. This is probably a slow way to do it.
            for (uword ii = 0; ii < g.n_rows; ii++) {
                for (uword jj = 0; jj < g.n_cols; jj++) {
                    ar << g(ii, jj);
                }
            }
        }

        template<class Archive>
        void load(Archive &ar, arma::Mat<std::complex<double>> &g, const unsigned int version) {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            unsigned int nRows = 0;
            unsigned int nCols = 0;
            std::complex<double> tempT;
            //Start with the length of the vector
            ar & nRows;
            ar & nCols;
            //Make sure the object is set to the right size.
            g.set_size(nRows, nCols);
            //Now read the elements out in the order they were read in. Probably a slow way to do this.
            for (uword ii = 0; ii < nRows; ii++) {
                for (uword jj = 0; jj < nCols; jj++) {
                    ar & tempT;
                    g(ii, jj) = tempT;
                }
            }

        }

        //Serialization code for arma::Mat<std::complex<double>>

        template<class Archive>
        void serialize(Archive &ar, arma::Mat<std::complex<float>> &g, const unsigned int version)
        {   //Defer the serialization to seperate load and save functions
            boost::serialization::split_free(ar, g, version);
        }

        template<class Archive>
        void save(Archive &ar, const arma::Mat<const std::complex<float>> &g, const unsigned int version)
        {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            //Start with the length of the vector
            ar & g.n_rows;
            ar & g.n_cols;
            //Now read out the data. This is probably a slow way to do it.
            for (uword ii = 0; ii < g.n_rows; ii++) {
                for (uword jj = 0; jj < g.n_cols; jj++) {
                    ar << g(ii, jj);
                }
            }
        }

        template<class Archive>
        void load(Archive &ar, arma::Mat<std::complex<float>> &g, const unsigned int version)
        {
            /* We need to save the data starting with the metadata to reconstruct an empty arma object.
             * Since this process isn't reversible (ie. we can't just repeat each of these steps in reverse using the &
             * operator, we have to break out into save and load. */
            unsigned int nRows = 0;
            unsigned int nCols = 0;
            std::complex<float> tempT;
            //Start with the length of the vector
            ar & nRows;
            ar & nCols;
            //Make sure the object is set to the right size.
            g.set_size(nRows,nCols);
            //Now read the elements out in the order they were read in. Probably a slow way to do this.
            for (uword ii = 0; ii < nRows; ii++) {
                for (uword jj = 0; jj < nCols; jj++) {
                    ar >> tempT;
                    g(ii, jj) = tempT;
                }
            }

        }

    } // namespace serialization
} // namespace boost


#endif //POWERGRID_ARMASERIALIZATION_HPP
