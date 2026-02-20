#Define files to be compiled
set(PG_FILES
	${PROJECT_SOURCE_DIR}/PowerGrid/Robject.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/QuadPenalty.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/TVPenalty.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/fftCPU.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/fftGPU.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/ftCpu.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/ftCpuWithGrads.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gdft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/GdftR2.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gfft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/TimeSegmentation.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/griddingSupport.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/gridding.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gnufft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/SENSE.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/pcSENSE.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/pcSenseTimeSeg.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/reconSolve.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.cpp
	#${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.h
	${PROJECT_SOURCE_DIR}/PowerGrid/directRecon.cpp
)

if(METAL_COMPUTE)
	list(APPEND PG_FILES
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalGridding.mm
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalVectorOps.mm
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/fftAccelerate.cpp)
	set_source_files_properties(
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalGridding.mm
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalVectorOps.mm
		PROPERTIES COMPILE_FLAGS "-fobjc-arc")
endif()

add_library(PGCommon SHARED ${PG_FILES})
if(METAL_COMPUTE AND TARGET MetalShaders)
	add_dependencies(PGCommon MetalShaders)
endif()
if(METAL_COMPUTE AND TARGET MetalVecOpsShaders)
	add_dependencies(PGCommon MetalVecOpsShaders)
endif()
target_link_libraries(PGCommon ${LIBS})
set_target_properties(PGCommon PROPERTIES LINK_DEPENDS_NO_SHARED true)

install(TARGETS PGCommon LIBRARY DESTINATION lib)
