//
// Created by Christian Messe on 19.01.21.
//

//
// Created by Christian Messe on 31.08.20.
//

#include <iostream>
#include "dlfcn.h"


#include "typedefs.hpp"
#include "constants.hpp"

#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"

#include "cl_EN_Parameters.hpp"
#include "cl_EN_Analysis.hpp"
#include "cl_EN_State.hpp"
#define private public
#define protected public
#include "cl_EN_Pump.hpp"
#undef private
#undef protected


using namespace belfem;
using namespace engine ;

real myfun( const real & aP, const real & aOF )
{
    Parameters tParams;

    tParams.set_chamber_pressure( aP );
    tParams.set_mixture_ratio( aOF );

    Analysis tGasGenerator( tParams );
    tGasGenerator.compute_injector( tParams.mixture_ratio(), tParams.chamber_pressure() );
    tGasGenerator.compute_total( 600.0 ) ;
    return tGasGenerator.total()->T() - 1000.0 ;
}

Communicator gComm;
Logger       gLog( 3 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    print_banner();


    /*Parameters tParams;

    tParams.set_chamber_pressure( 90e5 );
    tParams.set_mixture_ratio( 0.179 );

    tParams.set_fuel_and_oxidizer_temperatures( 114.0,  95.0 ) ;
    tParams.set_fuel_and_oxidizer( Fuel::LCH4, Oxidizer::LOX );

    Analysis tGasGenerator( tParams );
    tGasGenerator.compute_injector( tParams.mixture_ratio(), tParams.chamber_pressure() );
    tGasGenerator.compute_total( 600.0 ) ;
    tGasGenerator.total()->Tt() - 800.0 ;

    // tGasGenerator.total()->print() ;
    // tGasGenerator.combgas()->print() ;

    real tT = 200.0 ;

    while ( tT <= 1000.0 )
    {
        std::cout << tT << " " << tGasGenerator.combgas()->cp( tT, tParams.chamber_pressure()) << std::endl ;
        tT += 25.0;
    }*/

    /*real tP = 65e5 ;

    real x0 = 0.1 ;
    real x1 = 0.4 ;

    real f0 = myfun( tP, x0 );
    while( abs( x0 - x1 ) > 0.0000001 )
    {
        real x = 0.5 * ( x0 + x1 );
        real f = myfun( tP, x );
        if ( f*f0 > 0 )
        {
            x0 = x;
            f0 = f;
        }
        else
        {
            x1 = x ;
        }
        std::cout << x << " " << f << std::endl ;
    }
    exit( 0 ); */

    // crate a gas object
    Gas tLCH4( HelmholtzModel::Methane ) ;
    Gas tLOX( HelmholtzModel::Oxygen ) ;


    // create a pump object
    Pump tMethanePump( tLCH4 );
    Pump tOxygenPump( tLOX );

    real tNRPM = 20000 ;

    real tDotMf =  49.63 ;
    real tDotMo = 147.11 ;
    real tDw    = 0.04 ;
    real tDN    = 1.4 * tDw ;

     // entry conditions
    tMethanePump.set_entry( 111.0, 3.61e5 );
    tMethanePump.set_mass_flux( tDotMf );
    //tMethanePump.set_deltap( 117e5 ); // 77e5
    tMethanePump.set_pt2( 115.5e5 );

    tMethanePump.set_nrpm( tNRPM );
    tMethanePump.set_s2( 5e-3 ); // 1.2

    //tMethanePump.set_D1a( 0.113 ); // 0.01
    //tMethanePump.set_D2a( 0.200 );
    //tMethanePump.set_psi( 0.6 );

    tMethanePump.set_DN( tDN );  //tMethanePump.set_D1a( 0.1 );

    tMethanePump.set_z1( 4 );
    //tMethanePump.set_z2( 15 );

    tMethanePump.set_beta2( 18.0 );
    tMethanePump.compute() ;
    tMethanePump.print() ;

    // exit( 0 );

    tOxygenPump.set_entry( 93.0, 3.84e5 );
    tOxygenPump.set_mass_flux( tDotMo );
    // tOxygenPump.set_deltap( 117e5); // 58e5
    tOxygenPump.set_pt2( 115.5e5 );

    tOxygenPump.set_nrpm( tNRPM );
    tOxygenPump.set_s2( 5e-3 );
    tOxygenPump.set_DN( tDN ); // 0.007
    //tOxygenPump.set_psi( 0.4 );
    tOxygenPump.set_D2a( 0.150 );
    tOxygenPump.set_D1a( 0.132 );

    //tOxygenPump.set_z1( 3 ) ;
    //tOxygenPump.set_z2( 9 ) ; // 8

    tOxygenPump.set_beta2( 16.0 );
    tOxygenPump.compute() ;
    tOxygenPump.print() ;


    // close communicator
    return gComm.finalize();
}