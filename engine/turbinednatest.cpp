//
// Created by Christian Messe on 03.05.21.
//

#include "typedefs.hpp"
#include "constants.hpp"
#include "random.hpp"
#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "cl_EN_Parameters.hpp"
#include "cl_EN_Analysis.hpp"
#include "cl_EN_State.hpp"
#include "cl_EN_Turbine.hpp"
#include "cl_TensorMeshFactory.hpp"
#include "cl_Progressbar.hpp"
#include "cl_EN_Gene.hpp"

#include "cl_Cell.hpp"

using namespace belfem;
using namespace engine ;

Communicator gComm;
Logger       gLog( 3 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    print_banner();


//------------------------------------------------------------------------------



    // temperature of gas generator in K
    real tTgg = 750.0 ;

    // pressure of gas generator in pa
    real tPgg = 112.5e5 ;

    // shaft speed in RPM
    real tN = 50000 ;

    // mass flow in kg/s
    real tDotM = 0.664 ;

    // power in Watt
    real tPower = 0.536e6 ;

//------------------------------------------------------------------------------

    // number of individuals to keep
    index_t tNumKeep = 100 ;

    // number of indivuduals
    index_t tNumIndividuals = 50000 ;

    // max number of generations
    index_t tNumGenerations = 8 ;

    real tPhiMin = 0.25 ;
    real tPhiMax = 1.3 ;
    real tPsiMin = 1.75 ;
    real tPsiMax = 3.25 ;
    real tBdMin  = 0.04 ;
    real tBdMax  = 0.4 ;

//------------------------------------------------------------------------------

    // create the turbine
    Parameters tParams;

    tParams.set_fuel_and_oxidizer( Fuel::LNG, Oxidizer::LOX ) ;
    // tParams.set_fuel_and_oxidizer_temperatures( 115.7399 ,  97.4270 );

    Analysis tGasGenerator( tParams );

    std::cout << "computing gas mixture ... " << std::endl;
    real tOF = tGasGenerator.compute_gas_generator( tTgg,
                                                    tPgg,
                                                    0.1,
                                                    0.4 );

    std::cout << "OF " << tOF << std::endl ;

    // Gas tMethane("CH4");

    Turbine tTurbine( *tGasGenerator.combgas() );

    // reset counter
    tTurbine.set_n( tN );
    tTurbine.set_entry( tGasGenerator.total()->T(), tGasGenerator.total()->p() );

    tTurbine.set_massflow( tDotM );

    tTurbine.set_power( tPower );

//------------------------------------------------------------------------------

    // create original population
    Cell< Gene * > tGenes( tNumIndividuals, nullptr );
    for( uint k=0; k<tNumIndividuals; ++k )
    {
        tGenes( k ) = new Gene( tTurbine );
    }

//------------------------------------------------------------------------------

    // set values to initial data
    random_seed();
    for( Gene * tGene : tGenes )
    {
        real tPhi = tPhiMin + ( tPhiMax - tPhiMin ) * belfem::rand() * belfem::rand() * 2.0 ;
        real tPsi = tPsiMin + ( tPsiMax - tPsiMin ) * belfem::rand() * belfem::rand() * 2.0 ;
        real tBd  = tBdMin  + ( tBdMax - tBdMin ) * belfem::rand() * belfem::rand() * 2.0 ;

        tGene->set_values( tPhi, tPsi, tBd );
    }

//------------------------------------------------------------------------------


    for( index_t k=0; k<tNumGenerations; ++k )
    {
        Progressbar tProgress( tNumIndividuals );

        tProgress.reset();

        index_t tCount = 0;
        std::cout << "computing generation " << k << std::endl;
        for ( Gene * tGene : tGenes )
        {
            if ( tGene->alive() )
            {
                tGene->compute();

                if ( tGene->alive())
                {
                    ++tCount;
                }
            }
            tProgress.step();
        }
        tProgress.finish();



        sort( tGenes, opCompareFitness );

        std::cout << "alive individuals: " << tCount << " min. penalty: " << tGenes( 0 )->fitness() << std::endl;
        std::cout << " phi=" << tGenes( 0 )->phi() << " psi=" << tGenes( 0 )->psi() << " b/Dm=" << tGenes( 0 )->bd() << std::endl ;

        // number of parents
        index_t tNumParents = std::min( tCount, tNumKeep );

        tCount = tNumParents;

        // create next generation
        for ( index_t j = 0; j < tNumParents; ++j )
        {
            for ( index_t i = j; i < tNumParents; ++i )
            {
                // check limit
                if ( tCount == tNumIndividuals )
                {
                    break;
                }

                // create next generation
                tGenes( tCount++ )->inherit( tGenes( j ), tGenes( i ));
            }

            // check limit
            if ( tCount == tNumIndividuals )
            {
                break;
            }
        }

        // randomize rest
        tCount = std::min( tCount, index_t( 0.75 * ( real ) tNumIndividuals ) ) ;

        real tPhiBest = tGenes( 0 )->phi() ;
        real tPsiBest = tGenes( 0 )->psi() ;
        real tBdBest  = tGenes( 0 )->bd() ;

        for( index_t i=tCount; i<tNumIndividuals; ++i )
        {
            Gene * tGene = tGenes(  i  );

            real tPhi =  tPhiBest * ( 1.0 + 0.25 * ( belfem::rand() * belfem::rand() * 4.0  - 1.0 ) );
            real tPsi =  tPsiBest * ( 1.0 + 0.25 * ( belfem::rand() * belfem::rand() * 4.0  - 1.0 ) );
            real tBd  =  tBdBest  * ( 1.0 + 0.25 * ( belfem::rand() * belfem::rand() * 4.0  - 1.0 ) );

            tGene->set_values( tPhi, tPsi, tBd );
        }
    }

    tGenes( 0 )->compute() ;
    tTurbine.print();

//------------------------------------------------------------------------------

    // tidy up
    for( Gene * tGene : tGenes )
    {
        delete tGene ;
    }

//------------------------------------------------------------------------------
}
