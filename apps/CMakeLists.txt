set(POWERGRID_ISMRMRD
	PowerGridIsmrmrd.cpp
)
add_executable(PowerGridIsmrmrd ${POWERGRID_ISMRMRD})
target_link_libraries(PowerGridIsmrmrd PGCommon)

set(POWERGRID_PCSENSE
	PowerGridPcSense.cpp
)
add_executable(PowerGridPcSense ${POWERGRID_PCSENSE})
target_link_libraries(PowerGridPcSense PGCommon)

set(POWERGRID_PCSENSE_TIMESEG
	PowerGridPcSenseTimeSeg.cpp
)
add_executable(PowerGridPcSenseTimeSeg ${POWERGRID_PCSENSE_TIMESEG})
target_link_libraries(PowerGridPcSenseTimeSeg PGCommon)

add_executable(GenSyntheticISMRMRD GenSyntheticISMRMRD.cpp)
target_include_directories(GenSyntheticISMRMRD PRIVATE
    ${CMAKE_SOURCE_DIR}/PowerGrid/Tests)
target_link_libraries(GenSyntheticISMRMRD PGCommon)

add_subdirectory(lowRank)
add_subdirectory(reconSupport)

if (MPISupport)
    add_subdirectory(MPI)
endif (MPISupport)

install(TARGETS PowerGridIsmrmrd DESTINATION bin)
install(TARGETS PowerGridPcSense DESTINATION bin)
install(TARGETS PowerGridPcSenseTimeSeg DESTINATION bin)
install(TARGETS GenSyntheticISMRMRD DESTINATION bin)