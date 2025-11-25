//
// Created by Christian Messe on 31.08.20.
//

#include <iostream>

#include "typedefs.hpp"
#include "constants.hpp"

#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"

#include "cl_EN_Parameters.hpp"
#include "cl_EN_Analysis.hpp"
#include "cl_EN_State.hpp"
#include "fn_linspace.hpp"

using namespace belfem;
using namespace engine ;

Communicator gComm;
Logger       gLog( 3 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    print_banner();

    // create parameter object
    Parameters tParams ;

    // chamber pressure in Pa
    tParams.set_chamber_pressure( 90e5 );

    // mixture ratio
    tParams.set_mixture_ratio( 3.38489426 );

    tParams.set_throat_diameter( 1.0 );

    // fuel
    tParams.set_fuel_and_oxidizer( Fuel::LNG, Oxidizer::LOX );

    // if you don't set any temperatures, the default settings are
    // LH2 :  20 K
    // LOX :  90 K
    // LCH4: 110 K
    tParams.set_fuel_and_oxidizer_temperatures( 115., 95. );

    index_t tN = 22 ;
    Vector< real > tPc = linspace( 20.0, 230.0, tN );


    index_t tM = 31 ;
    Vector< real > tOf = linspace( 0.1, 0.4, tM );

    /*index_t tM = 20 ;
    Vector< real > tPe = linspace( 0.01, 0.2, tM ); */

    // create the engine
    Analysis tEngine( tParams );


    tEngine.fuel()->print();

    /*for( real pc : tPc )
    {
        for( real of : tOf )
        {
            real p = pc * 1e5 ;
            tEngine.compute_injector(of, p );
            tEngine.compute_total( 600.0 );
            real T = tEngine.total()->T() ;

            std::cout << pc << ", " << of << ", " << tEngine.total()->T() << " " << tEngine.total()->rho() << " " << tEngine.combgas()->R( T, p ) << " "
                << tEngine.combgas()->cp( 200.0, p ) << ", " << tEngine.combgas()->cp( 400.0, p ) << ", " << tEngine.combgas()->cp( 600.0, p )
                << tEngine.combgas()->cp( 800.0, p )  << tEngine.combgas()->cp( 1000.0, p ) << std::endl ;
        }
    } */

    // ambient pressure in Pa
    tParams.set_exit_pressure( 0.6e5 );

    // alternatively: set expansion ratio
    //tParams.set_expansion_ratio( 18.088 );


    // perform the analysis, and return the ISP ( uses parameter value if no argument is passed )


    real tOF = tEngine.find_best_mixture( 2.0, 6.0, IspMode::Sealevel );
    tEngine.run( tOF );

    std::cout << "OF " << tOF << std::endl ;
    // run engine with passed parameter
    //tEngine.run( tOF );

    // print the results
    tEngine.injector()->print();
    tEngine.total()->print();
    tEngine.throat()->print() ;
    tEngine.nozzle()->print() ;
    tEngine.print_performance() ;

    // print fluid properties
    Gas * tGas = tEngine.combgas() ;
    tGas->remix( tEngine.throat()->molar_fractions() );

    /*real tT = 1000 ;
    real tP = tEngine.throat()->p() ;

    while( tT < 3800 )
    {
        std::cout << tT << " " << tGas->cp( tT, tP ) << " " << tGas->mu( tT, tP ) << " " << tGas->lambda( tT, tP ) << std::endl ;
        tT += 100 ;
    }*/

    // close communicator
    return gComm.finalize();
}