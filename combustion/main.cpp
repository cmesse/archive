//
// Created by Christian Messe on 21.09.19.
//

#include <iostream>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "cl_Timer.hpp"

#include "cl_Vector.hpp"
#include "fn_linspace.hpp"
#include "fn_r2.hpp"

#include "cl_Spline.hpp"
#include "cl_SpMatrix.hpp"



#define private public
#define protected public
#include "cl_CN_Scheme.hpp"
#include "cl_CN_ReactionFactory.hpp"
#include "cl_CN_Reaction.hpp"
#include "cl_CN_Reaction_Duplicate.hpp"
#include "cl_CN_Reaction_Lindemann.hpp"
#include "cl_CN_Reaction_Troe.hpp"
#undef private
#undef protected

#include "fn_GT_data_path.hpp"
#include "cl_GT_RefGas.hpp"

using namespace belfem;
using namespace combustion;

Communicator gComm;
Logger       gLog( 3 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    string tPath = gastables::data_path() + "/jachimowski_short.inp";
    Scheme tScheme( tPath, Fuel::LH2, Oxidizer::AIR );

    real tT = 1171.64333660929;
    real tP = 40163.3592549777;
    real tMa = 2.54391153449378;
    real tU  = tScheme.combgas()->c( tT, tP ) * tMa ;

    Vector< real > tX( tScheme.combgas()->number_of_components(), 0.0 );


    tX( 0 )  = 3.1271749360219196e-2 ; // H
    tX( 1 )  = 0.14241524643246337e0 ; // H2
    tX( 2 )  = 4.4082034586038585e-2 ; // H2O
    tX( 3 )  = 2.7957329847635032e-6 ; // H2O2
    tX( 4 )  = 7.4322154527340527e-21 ; // HNO
    tX( 5 )  = 1.5296339195900822e-5  ; // HO2
    tX( 6 )  = 1.6456626281687479e-20 ; // N
    tX( 7 )  = 0.63349784386650232e0  ; // N2
    tX( 8 ) =  7.6818722639633952e-21 ; // NO
    tX( 9 ) =  5.0103363149995554e-21  ; // NO2
    tX( 10 )  = 7.1156516385142499e-3 ; // O
    tX( 11 )  = 0.13867136102297425e0  ; // O2
    tX( 12 )  = 2.9280210211072711e-3  ; // OH


    tScheme.combgas()->remix( tX );


    std::cout << "Running ... "<< std::endl;
    Timer tTimer;
    //for( uint k=0; k<1000; ++k )
    //{

        tScheme.compute( tT, tP, tU, 0.001 );
    //}
    uint tTime = tTimer.stop();

    //tScheme.mJacobi.print("J");

    std::cout << "Time " << tTime *1e-3 << std::endl;
    return 0;
}