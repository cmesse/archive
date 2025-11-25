//
// Created by Christian Messe on 25.02.20.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Gas.hpp"
#include "cl_Vector.hpp"

#include "fn_GT_data_path.hpp"
#include "cl_CN_Scheme.hpp"
#include "cl_CH_GeometryCombustor.hpp"
#include "cl_CH_ChannelODE.hpp"
#include "cl_ODE_Integrator.hpp"

#include "cl_CN_Injector.hpp"

#include "en_ODE_Type.hpp"

#include "cl_HDF5.hpp"
#include "fn_linspace.hpp"

using namespace belfem;
using namespace combustion;
using namespace channel;

Communicator gComm;
Logger       gLog( 3 );

//------------------------------------------------------------------------------

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    // create geometry with injector
    GeometryCombustor tGeo;


    // create the combustion scheme
    string tPath = gastables::data_path() + "/jachimowski.inp";
    Scheme tScheme( tPath, Fuel::LH2, Oxidizer::AIR );

    // get gas from scheme
    Gas * tGas = tScheme.combgas();

    Gas tFuel( "H2" );

//------------------------------------------------------------------------------

    // tGas->print();

    // reset scheme at initial condition
    tScheme.reset_combgas_mixture();

//------------------------------------------------------------------------------

    // initial x-Position
    real tX = 0.0;

    // position of injector
    real tXinj = tGeo.injector_position() ;

//------------------------------------------------------------------------------

    // initial conditions
    real tT3 = 1100; //820.0;
    real tP3 = 0.43e5 ; //42237.0;
    real tMa3 = 2.77; //3.75;
    real tTw  = 600.0;

    // duct height
    real tDuctHeight = 2.0 * tGeo.R( 0.0 );

    // initial massflow
    real tDotM_air = tGas->c( tT3, tP3 ) * tMa3 * tGas->rho( tT3, tP3 ) * tGeo.A( 0 );

    Injector tInjector(
            tScheme,
            tXinj,
            tDuctHeight,
            0.8, 15.0 );

    tInjector.set_phi( 0.6 );

    // tInjector.print() ;
    tInjector.set_oxidizer_massflow( tDotM_air );


//------------------------------------------------------------------------------

    // Create ODE object
    ChannelODE tODE( tGeo, *tScheme.combgas(), ChannelMode::Combustor );

    tODE.set_combustion( 0.0, 0.0 );

    // crete the integrator
    ode::Integrator tIntegrator( tODE, ode::Type::RK45 );

//------------------------------------------------------------------------------

    // set surface temperature
    tODE.set_wall_temperature( tTw );

    // initial state
    Vector < real > tY0( 3 );
    tY0( 0 ) = tGas->v( tT3, tP3 );
    tY0( 1 ) = tMa3 * tGas->c( tT3, tP3 );
    tY0( 2 ) = tT3;


    // state vector
    Vector< real > tY( tY0 );

    // derivative
    Vector< real > tdYdX( 3 );

    tY.print("Y0");

    tODE.compute( tX, tY,tdYdX );

    tdYdX.print("tdYdX");

    //tGas->print();

    tIntegrator.maxtime() = tGeo.injector_position();

    tIntegrator.timestep() = 0.01;

    while( tX < tXinj )
    {
        tIntegrator.step( tX, tY );

        real tT = tY( 2 );
        real tP = tGas->p( tY( 2 ), tY( 0 ) );
        real tMa = tY( 1 ) / tGas->c( tT, tP );
        std::cout << tX << " " << tT << " " << tP << " " << tMa << " " << tGas->s( tT, tP ) << std::endl;

    }

//------------------------------------------------------------------------------

    real tA1 =  tGeo.A( tXinj - BELFEM_EPSILON );
    real tA2 =  tGeo.A( tXinj + BELFEM_EPSILON );

    tX += 0.01;
    tIntegrator.timestep() = 0.001;


    real tT1 = tY( 2 );
    real tP1 = tGas->p( tY( 2 ), tY( 0 ) );
    real tU1 = tY( 1 );

//------------------------------------------------------------------------------
// Expansion of gas after step
//------------------------------------------------------------------------------
    real tT2;
    real tP2;
    real tU2;

    std::cout << "A " << tA1 << " " << tA2 << std::endl ;


    tGas->expand( tA1, tT1, tP1, tU1, tA2, tT2, tP2, tU2 );

    tY( 0 ) = tGas->v( tT2, tP2 );
    tY( 1 ) = tU2;
    tY( 2 ) = tT2 ;

    std::cout << "expand " << tT2 << " " << tP2 << " " << tU2 << std::endl ;

//------------------------------------------------------------------------------
// mixing of gas
//------------------------------------------------------------------------------


    real tDotM_fuel = tDotM_air / tInjector.of() ;
    std::cout << "dotm air  " << tDotM_air<< std::endl ;
    std::cout << "dotm fuel " << tDotM_fuel << std::endl ;
    std::cout << "R " << tGas->R( tT2, tP2 ) << std::endl ;
    // injection conditions of fuel
    real tTf = 600.0 ;
    real tUf = 300.0 ;

    // enthalpy of air
    real tHair = tGas->h( tT2, tP2 ) + 0.5 * tU2 * tU2 ;

    // enthalpy of fuel
    real tHfuel = tFuel.h( tTf, tP2 ) + 0.5 * tUf * tUf;

    // new velocity
    tU2 = ( tDotM_air * tU2 + tDotM_fuel * tUf ) / ( tDotM_fuel + tDotM_air );

    // new total enthalpy
    real tH2t = ( tHair * tDotM_air + tHfuel * tDotM_fuel ) / ( tDotM_fuel + tDotM_air ) ;

    // inject fuel into gas
    Vector< real > tMassFractions = tGas->mass_fractions() ;
    tMassFractions *= tDotM_air ;
    tMassFractions( tScheme.inert_fuel_index() ) += tDotM_fuel ;
    tMassFractions /= ( tDotM_fuel + tDotM_air );

    // remix gas
    tGas->remix_mass( tMassFractions );

    real tH2 = tH2t - 0.5 * tU2 * tU2 ;

    // compute new temperature
    tT2 = tGas->T_from_h( tH2, tP2 );
    std::cout << "R " << tGas->R( tT2, tP2 ) << std::endl ;

    std::cout << "enthalpy " << tHair << " " << tH2 << std::endl ;
    std::cout << "mix " << tT2 << " " << tP2 << " " << tU2 / tGas->c( tT2, tP2 )<< std::endl ;

    // exit( 0 );

    // tY( 2 ) = 1200 ;
    //tY( 2 ) = 1000.0 ;
    //tY( 1 ) = 3.3 * tScheme.combgas()->c( tY( 2 ), 0.6e5 );
    //tY( 0 ) = tScheme.combgas()->v( tY( 2 ), 0.6e5 );

//------------------------------------------------------------------------------

    tIntegrator.maxtime() = tGeo.length() ;

    Cell< Vector< real > > tCell;

    Cell< real > tTimes ;

//------------------------------------------------------------------------------

    // tIntegrator.run()

    // steps for fuel injection
    uint tNumInjectionSteps = 10 ;
    Vector< real > tXi( tNumInjectionSteps );
    linspace( 1./tNumInjectionSteps, 1.0, tNumInjectionSteps, tXi );

    real tdUdX = 0 ;

    while( tX < tGeo.length() )
    {
        //real tX0  ;
        real tX1 ;
        real tNextX = tX + 0.001;

        tIntegrator.maxtime() = tNextX ;

        //std::cout << tX << "A---------" << std::endl;

        while( tX < tNextX )
        {


            // compute state
            real tT = tY( 2 );
            real tP = tGas->p( tY( 2 ), tY( 0 ) );
            real tU = tY( 1 );

            // steps
            //tX0 = tX ;
            tX1 = tX ;
            real tDeltaX = tIntegrator.timestep() ;
            real tDeltaXk = tDeltaX / tNumInjectionSteps ;

            // mix fuel
            tInjector.inject( tX );

            // copy mass fractions
            tScheme.set_Y0( tScheme.combgas()->mass_fractions() );

            real tUb = tU ;

            for( uint k=0; k<tNumInjectionSteps; ++k )
            {
                // shift velocity
                real tUa = tUb ;

                // predict new velocity
                tUb = tUa + tdUdX * tDeltaXk;

                // average velocity
                real tUm = 0.5 * ( tUa + tUb );
                // next step
                tX1 += tDeltaXk ;

                // compute the combustion and increment temperature
                tT += tScheme.compute( tT, tP, tUm, tDeltaXk );

                tGas->remix_mass( tScheme.Y(), false, false );
            }

            real tdRdX = tScheme.delta_R() / tDeltaX  ;
            //real tdWdX = tScheme.delta_w() / tDeltaX ;

            // std::cout << "comb " << tT << " " << tdRdX * tScheme.combgas()->M( tT, tP ) << " " << tdWdX << std::endl ;

            // set combustion BCs
            tODE.set_combustion( tdRdX * tScheme.combgas()->M( tT, tP ) , 0 );

            tY( 2 ) = tT ;

            // perform integration
            tIntegrator.step( tX, tY );

            tdUdX = ( tY( 1 ) - tU ) / tDeltaX;
        }

        // remix gas
        tGas->remix_mass( tScheme.Y(), true, true );

        //std::cout << tX << "B---------" << std::endl;

        real tT = tY( 2 );
        real tP = tGas->p( tY( 2 ), tY( 0 ) );
        real tMa = tY( 1 ) / tGas->c( tT, tP );

        std::cout << tX << " " << tT << " " << tP << " " << tMa << " " << tGas->s( tT, tP )  << std::endl ;
        tCell.push( tY );
        tTimes.push( tX );
    }


//------------------------------------------------------------------------------
    
    // tIntegrator.get_data( tData );

    uint tN = tCell.size() ;

    Matrix < real > tData( tN, tODE.dimension() + 1 );

    for ( uint k=0; k< tN; ++k )
    {
        tData( k, 0 ) = tTimes( k );

        for ( uint i=0; i<tODE.dimension(); ++i )
        {
            tData( k, i+1 ) = tCell( k )( i );
        }
    }

//------------------------------------------------------------------------------

    // tIntegrator.clear_data();
    tTimes.clear();
    tCell.clear();

//------------------------------------------------------------------------------
// sava data to file

    HDF5 tFile( "data.hdf5", FileMode::NEW );

    tFile.save_data( "Combustor", tData );

    tFile.close();

//------------------------------------------------------------------------------

    return gComm.finalize();
}
