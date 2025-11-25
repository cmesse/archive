set( LIBNAME manta )

set( SOURCES
        cl_Manta.cpp
        cl_MantaTables.cpp
        cl_MantaSurface.cpp
        fn_create_manta_mesh.cpp
)

include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )


include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_NONFREE_SOURCE_DIR}/satellite )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        manta
        sparse
        mesh
        integration
        interpolation
        materials
        postproc
        kernel
        iwg
        kernel
)

set( EXECNAME mysat )
set( MAIN mysat.cpp)


include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
