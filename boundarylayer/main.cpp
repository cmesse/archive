//
// Created by Christian Messe on 29.11.19.
//

#include <iostream>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

#include "cl_Gas.hpp"
#include "cl_BL_State.hpp"
#include "fn_BL_Eckert.hpp"
#include "fn_BL_VanDriest.hpp"

using namespace belfem;
using namespace boundarylayer;
Communicator gComm;
Logger       gLog( 4 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    Gas tAir;

    // reference conditions for temperature and pressure @ 30 km
    real tT0 = 226.51;
    real tP0 = 1196.7;
    real tU0 = tAir.c( tT0, tP0 ) * 8.0 ;

    // perform shock
    // inclinded plate
    real tAlpha = 60.0 * constant::pi / 180.0 ;

    // surface temperature
    real tTw = 600.0 ;

    // shock angle
    real tBeta ;

    // conditions after shock
    real tT1;
    real tP1;
    real tU1;

    tAir.shock(tT0, tP0, tU0, tAlpha, tT1, tP1, tU1, tBeta );


    // create state
    State tState( tAir );

    // set state parameters
    tState.compute( tT1, tP1, tU1 );

    // set wall temperature
    tState.set_wall_temperature( tTw );

    // compute wall state
    tState.compute_wall_state() ;

    real tX = 1.00 ;

    while ( tX < 20 )
    {
        // Eckert solution ( may be 20% higher than van driest, see Meador & Smart )
        eckert( tState, tX, true );

        real tDotQ1 = tState.dot_q() *1E-6;

        vandriest( tState, tX );
        real tDotQ2 = tState.dot_q()  * 1e-6;

        std::cout << tX << " " << tDotQ1 << " " << tDotQ2 << std::endl;

        tX += 0.5 ;
    }

    return gComm.finalize();
}