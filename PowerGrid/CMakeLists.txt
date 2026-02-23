#Define files to be compiled
set(PG_FILES
	# Solvers
	${PROJECT_SOURCE_DIR}/PowerGrid/Solvers/Robject.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Solvers/QuadPenalty.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Solvers/TVPenalty.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Solvers/reconSolve.cpp
	# FFT
	${PROJECT_SOURCE_DIR}/PowerGrid/FFT/fftCPU.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/FFT/fftGPU.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/FFT/ftCpu.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/FFT/ftCpuWithGrads.cpp
	# Operators
	${PROJECT_SOURCE_DIR}/PowerGrid/Operators/Gdft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Operators/GdftR2.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Operators/Gfft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Operators/Gnufft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Operators/SENSE.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Operators/pcSENSE.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Operators/pcSenseTimeSeg.cpp
	# Gridding
	${PROJECT_SOURCE_DIR}/PowerGrid/Gridding/TimeSegmentation.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gridding/griddingSupport.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gridding/gridding.cpp
	# IO
	${PROJECT_SOURCE_DIR}/PowerGrid/IO/acqTracking.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/IO/directRecon.cpp
)

if(METAL_COMPUTE)
	list(APPEND PG_FILES
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalGridding.mm
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalDFT.mm
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalNufftPipeline.mm
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalVectorOps.mm
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/fftAccelerate.cpp)
	set_source_files_properties(
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalGridding.mm
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalDFT.mm
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalNufftPipeline.mm
		${PROJECT_SOURCE_DIR}/PowerGrid/Metal/MetalVectorOps.mm
		PROPERTIES COMPILE_FLAGS "-fobjc-arc")
endif()

add_library(PGCommon SHARED ${PG_FILES})
target_include_directories(PGCommon PUBLIC ${PROJECT_SOURCE_DIR}/PowerGrid)
if(METAL_COMPUTE AND TARGET MetalShaders)
	add_dependencies(PGCommon MetalShaders)
endif()
target_link_libraries(PGCommon ${LIBS})
set_target_properties(PGCommon PROPERTIES LINK_DEPENDS_NO_SHARED true)

install(TARGETS PGCommon LIBRARY DESTINATION lib)
