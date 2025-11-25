//
// Created by christian on 9/11/24.
//

#include <iostream>
#include "commtools.hpp"

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_SpMatrix.hpp"
#include "fn_rcond.hpp"

#include "cl_MantaSurface.hpp"
#include "cl_MantaTables.hpp"
#include "cl_Manta.hpp"
#include "fn_create_manta_mesh.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 3 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    std::cout << "hello world" << std::endl ;


    Mesh * tMesh = create_manta_mesh( "manta.hdf5") ;

    std::cout << std::endl ;
    std::cout << std::endl ;
    MantaTables tTables( *tMesh );

    if( comm_rank() == 0 )
    {
        std::cout << "load database ... " << std::endl ;
        tTables.load_database( "manta.hdf5" );

        std::cout << "create matrix ... " << std::endl ;

        std::cout << "done" << std::endl ;

        real & tTime = tMesh->time_stamp() ;
        tTime = 0 ;
        tTables.interpolate_geometry_info( tTime );
        tTables.compute_environment( tTime );

        tTables.solve_infrared( tTime );
        tTables.solve_solar( tTime );
        tMesh->save( "manta.exo");

    }

    delete tMesh ;

    // close communicator
    return gComm.finalize();
}
