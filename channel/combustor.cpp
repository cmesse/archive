//
// Created by Christian Messe on 30.08.20.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Gas.hpp"
#include "cl_Vector.hpp"


#include "cl_HDF5.hpp"
#include "fn_linspace.hpp"
#include "fn_norm.hpp"
#include "fn_sum.hpp"
#include "cl_GT_RefGas.hpp"
#include "fn_dot.hpp"
#include "cl_CN_Scheme.hpp"
#include "GT_globals.hpp"
#include "fn_GT_data_path.hpp"
using namespace belfem;
Communicator gComm;
Logger       gLog( 3 );

//------------------------------------------------------------------------------

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );
//------------------------------------------------------------------------------


    /*Cell< string > tSpecies = {
            "CH4",
            "C2H6",
            "C3H8",
            "CO2",
            "O2",
            "CO",
            "CO2",
            "COOH",
            "H",
            "H2",
            "H2O",
            "H2O2",
            "HCHO,formaldehy",
            "HCO",
            "HCOOH",
            "HO2",
            "O",
            "OH",
            "CH2",
            "CH3",
            "CH2O",
            "CH3O",
            "C2H3",
            "C2H4",
            "C2H5",
            "C2H6",
            "CH3O2",
            "O3"}; */
    Cell< string > tSpecies  = {
            "CH4",
            "C2H6",
            "C3H8",
            "CO2",
            "CH2O",
            "CO",
            "COOH",
            "H",
            "H2",
            "H2O",
            "H2O2",
            "HCHO,formaldehy",
            "HCO",
            "HCOOH",
            "O",
            "OH" };

    //string tPath = gastables::data_path() + "/zhukov_kong.inp";
    //combustion::Scheme tScheme( tPath, Fuel::LH2, Oxidizer::LOX );

    real tT = 100.0 ; // <<-- injection temperature
    real tP = 70e5 ;

    real tAt = 5491.743034e-6 ;
    real tA  = 10647.2536 ;

    real tOF = 3.2 ;
    uint tN = tSpecies.size() ;

    Vector< real > tX( tN, 0.0 );
    tX( 1 ) = 1.0 ;

    // get gas from scheme
    Gas * tGas = new Gas( tSpecies, tX );

    const uint tOxidizer = 4;

    tX( 0 ) = 0.8 / tGas->component( 0 )->M() ;
    tX( 1 ) = 0.1 / tGas->component( 1 )->M() ;
    tX( 2 ) = 0.1 / tGas->component( 2 )->M() ;
    tX( 3 ) = 0.0 / tGas->component( 3 )->M() ;

    tX( tOxidizer ) = tOF / tGas->component( tOxidizer )->M() ;
    tX /= sum( tX );

//------------------------------------------------------------------------------
// Computing the combustion temperature
//------------------------------------------------------------------------------
    tGas->remix( tX );
    real tH = tGas->h( tT, tP );

    real tOmega = 0.5 ;

    real tDeltaT = BELFEM_REAL_MAX ;
    while( std::abs( tDeltaT ) > 1e-5 )
    {
        tGas->remix_to_equilibrium( tT, tP, true, false );

        tDeltaT = tGas->T_from_h( tH, tP ) - tT ;
        tT += tOmega * tDeltaT ;
        std::cout << tT << std::endl ;
    }

    real tS = tGas->s( tT, tP );
    real tC = tGas->c( tT, tP );

    std::cout << "T: " << tT << std::endl ;
    real tK = tGas->gamma( tT, tP ) ;
    std::cout << "kappa " << tK << std::endl ;
    std::cout << "M " << tGas->M( tT, tP ) *1000 << std::endl ;
    real tRho = tGas->rho( tT, tP );

    for( uint k=0; k<tN; ++k )
    {
        std::cout << tSpecies( k ) << " " << tGas->mass_fraction( k ) << std::endl ;
    }
//------------------------------------------------------------------------------
// throat state
//------------------------------------------------------------------------------

    // initial guess at throat
    real tTt = tT / ( 1.0 + 0.5 * ( tK - 1.0 ) );

    // initial guess for pressure at throat
    real tPt = tGas->isen_p( tT, tP, tTt );

    real tWt ;
    real tHt ;

    // old temperature
    real tTt0 = BELFEM_REAL_MAX ;

    while( std::abs( tTt0 - tTt ) > 1e-4 )
    {
        // compute velocity at throat
        tWt = tGas->c( tTt, tPt );

        // correct temperature
        tHt = tH - 0.5 * tWt * tWt;

        tTt0 = tTt ;

        // fixme: T_from_h does not work for real gas
        real tDeltaT = BELFEM_REAL_MAX ;
        real tCp = tGas->cp( tTt, tPt);
        real tH0 = tGas->h( tTt, tPt );

        while( std::abs( tDeltaT ) > 1e-6 )
        {
            tDeltaT = ( tH0 + tCp * ( tTt - tTt0 ) - tHt ) / tCp;
            tTt -= 0.5 *  tDeltaT;
        }

        real tDeltaS = BELFEM_REAL_MAX;

        while ( std::abs( tDeltaS ) > 0.0001 )
        {
            tDeltaS = tGas->s( tTt, tPt ) - tS;
            tPt -= 0.5 * tDeltaS / tGas->dsdp( tTt, tPt );
        }


        // compute equilibrium composition at throat
        Vector< real > tX0 = tGas->molar_fractions() ;
        tGas->remix_to_equilibrium( tTt, tPt, true, false );
        tX0 += tGas->molar_fractions() ;
        tX0 *= 0.5 ;
        tGas->remix( tX0 );

    }
    tGas->remix_to_equilibrium( tTt, tPt );

    std::cout << "Tt " << tTt << " " << tPt << " " << tGas->M( tT, tP ) * 1000 << " " << tGas->cp( tTt, tPt ) << std::endl;

    real tDotM = tGas->c( tTt, tPt ) * tGas->rho( tTt, tPt ) * tAt ;

    std::cout << "tDotM " << tDotM << " " << tGas->c( tTt, tPt ) << " " << tGas->rho( tTt, tPt ) << std::endl ;

    for( uint k=0; k<tN; ++k )
    {
        std::cout << tSpecies( k ) << " " << tGas->mass_fraction( k ) <<std::endl ;
    }

//------------------------------------------------------------------------------
// Inlet state
//------------------------------------------------------------------------------

    // guess mach number
    tK = tGas->gamma( tTt, tPt );

    real tUa = 0.0 ;
    real tUb = tC ;

    real tFa = tRho * tA * tUa - tDotM ;
    real tFb = tRho * tA * tUb - tDotM ;
    real tF = BELFEM_REAL_MAX ;
    real tU ;

    while( std::abs( tF ) > 1e-6 )
    {

        tU = tUa - tFa * ( tUb - tUa ) / ( tFb - tFa ) ;

        tF = tRho * tA * tU - tDotM ;


        if ( tF * tFa < 0.0 )
        {
            tUa = tU;
            tFa = tF;
        }
        else
        {
            tUb = tU;
            tFb = tF;
        }
    }
    std::cout << "u0 " << tU << std::endl ;

    delete tGas ;
    return gComm.finalize();
}