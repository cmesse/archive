//
// Created by Christian Messe on 01.04.21.
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
#include "cl_EN_Turbine.hpp"
#include "cl_TensorMeshFactory.hpp"
#include "cl_Progressbar.hpp"
#include "cl_EN_Gene.hpp"
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

    //Gene tGene ;

    //tGene.set_values( 1.1, 2.5, 0.05 );
    //std::cout << "psi " << tGene.psi() << std::endl ;

    exit( 0 );

    Parameters tParams;


    tParams.set_fuel_and_oxidizer_temperatures( 115.7399 ,  97.4270 );

    Analysis tGasGenerator( tParams );

    real tOF = tGasGenerator.compute_gas_generator( 750.0,
                                                    112.5e5,
                                                    0.1,
                                                    0.4 );

    std::cout << "OF " << tOF << std::endl ;

    Turbine tTurbine( *tGasGenerator.combgas() );

    // reset counter
    tTurbine.set_n( 20000 );
    tTurbine.set_entry( tGasGenerator.total()->T(), tGasGenerator.total()->p() );

    // tTurbine.set_massflow( 6.454980087  );
    tTurbine.set_massflow( 6.454980087 * 2.5 / 3.5  );

    tTurbine.set_power( 3.9e6 );
    //tTurbine.set_pitch_chord_ratio( 2.0 );
    //tTurbine.set_Z2( 60 );
    //tTurbine.set_alpha2( 90.0 / 180.0 * constant::pi );

    tTurbine.set_psi( 2.513 );
    //tTurbine.set_epsilon( 0.125 );
    tTurbine.set_phi( 0.769 );
    tTurbine.set_bd( 72. / 554. );
    //tTurbine.set_b( 0.069 );
    tTurbine.compute() ;
    tTurbine.print() ;

    // exit( 0 );



    TensorMeshFactory tFactory ;

    Mesh * tMesh = tFactory.create_tensor_mesh( { 50, 50, 50 }, {0.75, 2.495, 1.25}, { 0.8, 2.525, 1.35 }, 1 ) ;
    //Mesh * tMesh = tFactory.create_tensor_mesh( { 50, 50, 50 }, {0.75, 2.25, 0.5}, { 1.0, 2.75, 1.35 }, 1 ) ;


    Vector< real > & tPsi = tMesh->create_field("psi");

    Vector< real > & tPhi0 = tMesh->create_field("phi0");
    Vector< real > & tPhi1 = tMesh->create_field("phi1");
    Vector< real > & tPhi2 = tMesh->create_field("phi2");

    Vector< real > & tMa0 = tMesh->create_field("Ma0");

    Vector< real > & tTt1 = tMesh->create_field("Tt1");
    Vector< real > & tpt1 = tMesh->create_field("pt1");
    Vector< real > & tMa1 = tMesh->create_field("Ma1");

    Vector< real > & tTt1r = tMesh->create_field("Tt1r");
    Vector< real > & tpt1r = tMesh->create_field("pt1r");
    Vector< real > & tMa1r = tMesh->create_field("Ma1r");

    Vector< real > & tGamma1 = tMesh->create_field("gamma1");

    Vector< real > & tTt2 = tMesh->create_field("Tt2");
    Vector< real > & tpt2 = tMesh->create_field("pt2");
    Vector< real > & tMa2 = tMesh->create_field("Ma2");

    Vector< real > & tTt2r = tMesh->create_field("Tt2r");
    Vector< real > & tpt2r = tMesh->create_field("pt2r");
    Vector< real > & tMa2r = tMesh->create_field("Ma2r");
    Vector< real > & tGamma2 = tMesh->create_field("gamma2");

    Vector< real > & tDm = tMesh->create_field("Dm");
    Vector< real > & tB = tMesh->create_field("b");
    Vector< real > & tBd = tMesh->create_field("bD");

    Vector< real > & tEta = tMesh->create_field("eta");
    Vector< real > & tR = tMesh->create_field("r");
    Vector< real > & tHaller = tMesh->create_field("haller");
    Vector< real > & tEpsilon = tMesh->create_field("epsilon");

    Vector< real > & tAlpha1= tMesh->create_field("alpha1");
    Vector< real > & tBeta1 = tMesh->create_field("beta1");
    Vector< real > & tBeta2  = tMesh->create_field("beta2");
    Vector< real > & tMu  = tMesh->create_field("mu");
    Vector< real > & tBladeError  = tMesh->create_field("bladerror");
    Vector< real > & tError  = tMesh->create_field("errorcode");

    Vector< real > & tZ2  = tMesh->create_field("numblades");
    Vector< real > & tPitch  = tMesh->create_field("pitch");
    Vector< real > & tPitchChordRatio  = tMesh->create_field("pitchchordratio");

    index_t tNumNodes = tMesh->number_of_nodes() ;

    Cell< mesh::Node * > & tNodes = tMesh->nodes() ;

    Progressbar tProgress( tNumNodes );
    tProgress.reset() ;

    for( index_t k=0; k<tNumNodes; ++k )
    {
        tProgress.step( k ) ;
        tTurbine.reset();
        tTurbine.set_phi( tNodes( k )->x() );
        tTurbine.set_psi( tNodes( k )->y());
        tTurbine.set_bd( tNodes( k )->z() * 0.1 );

        tTurbine.compute();

        tPsi( k )  = tTurbine.psi();

        tPhi0( k ) = tTurbine.phi0() ;
        tPhi1( k ) = tTurbine.phi1() ;
        tPhi2( k ) = tTurbine.phi2() ;
        tMu( k ) = tTurbine.phi2() / tTurbine.phi1() ;

        tBd( k ) = tTurbine.b() / tTurbine.Dm() ;
        tB( k ) = tTurbine.b();
        tDm( k ) = tTurbine.Dm();

        tMa0( k ) = tTurbine.nozzle_entry()->Ma();



        tTt1( k ) = tTurbine.turbine_entry()->Tt() ;
        tpt1( k ) = tTurbine.turbine_entry()->pt() * 1e-5 ;
        tMa1( k ) = tTurbine.turbine_entry()->Ma();

        tTt1r( k ) = tTurbine.turbine_entry_rotating()->Tt() ;
        tpt1r( k ) = tTurbine.turbine_entry_rotating()->pt() * 1e-5 ;
        tMa1r( k ) = tTurbine.turbine_entry_rotating()->Ma();

        tGamma1( k ) =
                ( tTurbine.turbine_entry()->Tt() / tTurbine.turbine_entry()->T() - 1.0 )
                * 2.0 / ( tTurbine.turbine_entry()->Ma() *  tTurbine.turbine_entry()->Ma() )
                + 1.0 ;

        tTt2( k ) = tTurbine.turbine_discharge()->Tt();
        tpt2( k ) = tTurbine.turbine_discharge()->pt() * 1e-5 ;
        tMa2( k ) = tTurbine.turbine_discharge()->Ma();

        tTt2r( k ) = tTurbine.turbine_discharge_rotating()->Tt();
        tpt2r( k ) = tTurbine.turbine_discharge_rotating()->pt() * 1e-5 ;
        tMa2r( k ) = tTurbine.turbine_discharge_rotating()->Ma();

        tGamma2( k ) =
                ( tTurbine.turbine_discharge()->Tt() / tTurbine.turbine_discharge()->T() - 1.0 )
                * 2.0 / ( tTurbine.turbine_discharge()->Ma() *  tTurbine.turbine_discharge()->Ma() )
                + 1.0 ;

        tEta( k ) = tTurbine.eta() ;
        tEpsilon( k ) = tTurbine.epsilon();
        tHaller( k ) = tTurbine.haller() ;
        tR( k ) = tTurbine.reaction();

        tAlpha1( k ) = tTurbine.alpha1() / constant::deg ;
        tBeta1( k ) = tTurbine.beta1() / constant::deg ;
        tBeta2( k ) = tTurbine.beta2() / constant::deg ;

        tBladeError( k ) = tTurbine.blade_entry_error() * 1000 ;
        tError( k ) = ( real ) tTurbine.error_code() ;
        tZ2( k ) = ( real ) tTurbine.Z2() ;
        tPitchChordRatio( k ) = tTurbine.pitch_chord_ratio() ;
        tPitch( k ) = tTurbine.pitch() ;
    }

    tProgress.finish() ;
    tMesh->save( "/tmp/mesh.exo");


    // tTurbine.print();

}
