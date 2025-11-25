set( LIBNAME engine )

set( SOURCES
        cl_EN_Parameters.cpp
        cl_EN_State.cpp
        cl_EN_StepanoffChart.cpp
        cl_EN_Pump.cpp
        cl_EN_Turbine.cpp
        cl_EN_Analysis.cpp
        cl_EN_PumpArguments.cpp
        cl_EN_PumpUserLibrary.cpp
        cl_EN_Gene.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
include_directories( ${BELFEM_NONFREE_SOURCE_DIR}/combustion )
include_directories( ${BELFEM_NONFREE_SOURCE_DIR}/engine )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        mesh
        sparse
        gastables
        gasmodels
        combustion
        engine )

set( EXECNAME pump )
set( MAIN     pump.cpp )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

if( USE_EXAMPLES )
    set( EXECNAME engine )
    set( MAIN     main.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME pumptest )
    set( MAIN     pumptest.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME turbinetest )
    set( MAIN     turbinetest.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME turbinednatest )
    set( MAIN     turbinednatest.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

endif()