//
// Created by Christian Messe on 29.04.21.
//
#include <iostream>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "assert.hpp"

#include "cl_Gas.hpp"

#include "cl_EN_PumpArguments.hpp"
#include "cl_EN_Pump.hpp"
#include "cl_EN_PumpUserLibrary.hpp"

using namespace belfem;
using namespace engine ;

Communicator gComm;
Logger       gLog( 3 );

//------------------------------------------------------------------------------

/**
 * Make sure that this executable is not called in parallel
 * @return
 */
void
check_single_core()
{
    BELFEM_ERROR( comm_size() == 1,
                 "This program can't be executed in parallel mode" );
}

//------------------------------------------------------------------------------

void
print_usage()
{
    // todo: print proper usage
    std::cout << "Usage: pump -l $library -f $fluid -s $symbol" << std::endl ;
}

//------------------------------------------------------------------------------

void
print_help()
{
    // todo: print proper helptag
    print_usage();
}

//------------------------------------------------------------------------------

void
run( const string & aLibraryPath,
     const string & aSymbolName,
     const HelmholtzModel & aFluid )
{
    std::cout << "Library : " << aLibraryPath << std::endl ;

    // check sanity
    BELFEM_ERROR( aLibraryPath.length() > 0, "No userfile specified" );

    BELFEM_ERROR( aFluid != HelmholtzModel::UNDEFINED , "No fluid specified" );

    // create fuel object
    Gas tFluid( aFluid );

    // create the pump object
    Pump tPump( tFluid );

    // create the user library
    Library tLibrary( aLibraryPath );

    // load the function from the library
    BELFEM_PUMP_USER_FUNCTION tFunction = tLibrary.load_function( aSymbolName );

    // use the settings on the pump
    tFunction( tPump );

    // compute the data
    tPump.compute();

    // print the result
    tPump.print();
}

//------------------------------------------------------------------------------

int
main( int    argc,
      char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    // do not run this in single core
    check_single_core();

    PumpArguments tArgs(  argc, argv );

    std::cout << "Hallo Welt" << std::endl;

    switch ( tArgs.state() )
    {
        case( RunState::PrintUsage ) :
        {
            print_usage();
            break ;
        }
        case( RunState::PrintHelp ) :
        {
            print_help();
            break;
        }
        case ( RunState::Compute ) :
        {
            run( tArgs.library_path(), tArgs.symbol_name(), tArgs.fluid() );
            break ;
        }
        default:
        {
            BELFEM_ERROR( false, "Something went wrong");
        }
    }

    // close communicator
    return gComm.finalize();
}