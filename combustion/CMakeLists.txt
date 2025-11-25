set( LIBNAME combustion )

set( SOURCES
        CN_Enums.cpp
        fn_CN_arrhenius.cpp
        fn_CN_molars_from_string.cpp
        cl_CN_Entry.cpp
        cl_CN_Chemkin.cpp
        cl_CN_Scheme.cpp
        cl_CN_Reaction.cpp
        cl_CN_Reaction_Arrhenius.cpp
        cl_CN_Reaction_Duplicate.cpp
        cl_CN_Reaction_Lindemann.cpp
        cl_CN_Reaction_Troe.cpp
        cl_CN_ReactionFactory.cpp
        cl_CN_Injector.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
include_directories( ${BELFEM_NONFREE_SOURCE_DIR}/combustion )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        sparse
        gastables
        gasmodels
        combustion )

if( USE_EXAMPLES )
    set( EXECNAME combustion )
    set( MAIN     main.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()