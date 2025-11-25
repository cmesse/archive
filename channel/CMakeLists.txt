set( LIBNAME channel )

set( SOURCES
        CH_Enums.cpp
        cl_CH_Geometry.cpp
        cl_CH_GeometryCombustor.cpp
        cl_CH_GeometryCylinderCombustor.cpp
        cl_CH_GeometryNozzle.cpp
        cl_CH_ChannelODE.cpp
        cl_CH_Factory.cpp
        cl_CH_Wall.cpp
        cl_CH_Boundarylayer.cpp
        cl_CH_Segment.cpp
        cl_CH_Element.cpp
        cl_Channel.cpp
        cl_IsotropicChannel.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/bezier )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/ode )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
include_directories( ${BELFEM_NONFREE_SOURCE_DIR}/boundarylayer )
include_directories( ${BELFEM_NONFREE_SOURCE_DIR}/combustion )
include_directories( ${BELFEM_NONFREE_SOURCE_DIR}/channel )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        sparse
        bezier
        gastables
        gasmodels
        mesh
        integration
        boundarylayer
        combustion
        ode
        interpolation
        iwg
        postproc
        integration
        channel )

if( USE_EXAMPLES )
    set( EXECNAME itlr_combustor )
    set( MAIN     itlr_combustor.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME combustor )
    set( MAIN combustor.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

endif()