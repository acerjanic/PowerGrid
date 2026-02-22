/// @file acqTracking.h
/// @brief Acquisition parameter tracker for ISMRMRD datasets.

//
// Created by acerja2 on 10/22/17.
//

#ifndef POWERGRID_ACQTRACKING_H
#define POWERGRID_ACQTRACKING_H

#include "Core/PowerGrid.h"

/// @brief Tracks ISMRMRD acquisition counters to determine dataset dimensions.
///
/// Scans an ISMRMRD dataset to find the maximum index values for each
/// acquisition dimension (shots, partitions, slices, repetitions, averages,
/// echoes, phases) and stores a lookup array mapping each acquisition to its
/// position in the full dataset.
class acqTracking {

public:
	~acqTracking() {};

	/// @brief Construct an acquisition tracker by scanning all acquisitions in a dataset.
	///
	/// @param d    Pointer to the open ISMRMRD dataset.
	/// @param hdr  Reference to the parsed ISMRMRD header.
	acqTracking(ISMRMRD::Dataset *d, ISMRMRD::IsmrmrdHeader &hdr);

	/// @brief Maximum shot (interleave) index + 1.
	int NShotMax;
	/// @brief Maximum partition (3-D encode) index + 1.
	int NParMax;
	/// @brief Maximum slice index + 1.
	int NSliceMax;
	/// @brief Maximum repetition index + 1.
	int NRepMax;
	/// @brief Maximum average index + 1.
	int NAvgMax;
	/// @brief Maximum echo index + 1.
	int NEchoMax;
	/// @brief Maximum phase index + 1.
	int NPhaseMax;

	/// @brief Lookup array mapping each acquisition to its sequential index.
	ISMRMRD::NDArray<int> acqArray;
private:
	ISMRMRD::Dataset *d;



};

#endif //POWERGRID_ACQTRACKING_H
