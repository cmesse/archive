set( LIBNAME boundarylayer )


set( SOURCES
        cl_BL_State.cpp
        fn_BL_Spalding.cpp
        fn_BL_g_plus.cpp
        fn_BL_ReferenceTempertaure.cpp
        fn_BL_Kays_Crawford.cpp
        fn_BL_cf_flatplate_inc_laminar.cpp
        fn_BL_cf_flatplate_inc_turbulent.cpp
        fn_BL_virtual_origin.cpp
        fn_BL_Eckert.cpp
        fn_BL_Eckert_Xref.cpp
        fn_BL_VanDirest.cpp
        fn_BL_VanDriest_Xref.cpp
        fn_BL_Moody.cpp
        cl_BL_Panel.cpp
        cl_BL_StagnationPoint.cpp
        cl_BL_StreamLine.cpp
        )


include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
include_directories( ${BELFEM_NONFREE_SOURCE_DIR}/boundarylayer )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )

set ( LIBLIST
sparse
gastables
gasmodels
mesh
integration
interpolation
postproc
boundarylayer )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

if( USE_EXAMPLES )
set( EXECNAME wall )
set( MAIN     main.cpp )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()