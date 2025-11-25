//
// Created by Christian Messe on 02.05.21.
//
#include "assert.hpp"
#include "cl_EN_Gene.hpp"
#include "random.hpp"

namespace belfem
{
    namespace engine
    {

//------------------------------------------------------------------------------

        Gene::Gene( Turbine & aTurbine ) :
            mTurbine( aTurbine )
        {

        }

//------------------------------------------------------------------------------

        void
        Gene::reset()
        {
            mDNA.reset();
        }

//------------------------------------------------------------------------------

        void
        Gene::set_values( const real & aPhi, const real & aPsi, const real & aBD )
        {
            Vector< real > tValues( 3 );
            tValues( 0 ) = aPhi ;
            tValues( 1 ) = aPsi ;
            tValues( 2 ) = aBD ;
            mDNA.set_values( tValues );
        }

//------------------------------------------------------------------------------

        void
        Gene::set_values( const Vector< real > & aValues )
        {
            mDNA.set_values( aValues );
        }

//------------------------------------------------------------------------------

        void
        Gene::inherit( const Gene * aMom, const Gene * aDad )
        {

            // length of DNA string
            index_t tNumGenes = mDNA.data().size() ;

            // compute split index
            index_t tSplit = index_t( belfem::rand() * tNumGenes );

            // compute mutation index
            index_t tMutate = index_t( belfem::rand() * tNumGenes );


            // inherit genes from mom
            for( index_t k=0; k<tSplit; ++k )
            {
                if ( aMom->test( k ) )
                {
                    mDNA.set( k );
                }
                else
                {
                    mDNA.reset( k );
                }
            }

            // inherit genes from dad
            for( index_t k=tSplit; k<tNumGenes; ++k )
            {
                if ( aDad->test( k ) )
                {
                    mDNA.set( k );
                }
                else
                {
                    mDNA.reset( k );
                }
            }

            // mutate
            mDNA.flip( tMutate );

            // sel life flag to on
            mDNA.resurrect();
        }

//------------------------------------------------------------------------------

        void
        Gene::compute()
        {
            // get values
            real tPhi = this->phi() ;
            real tPsi = this->psi() ;
            real tBD  = this->bd() ;

            mDNA.resurrect() ;

            // check if the values make sense
            if( tPhi < 0.2 || tPhi > 1.3 )
            {
                mDNA.kill() ;
                return;
            }
            if( tPsi < 1.75 || tPsi > 3.25 )
            {
                mDNA.kill() ;
                return;
            }
            if( tBD < 0.04 || tBD > 0.4 )
            {
                mDNA.kill() ;
                return;
            }

            // check if values make sense
            if ( mDNA.alive() )
            {
                // set values
                mTurbine.set_phi( tPhi );
                mTurbine.set_psi( tPsi );
                mTurbine.set_bd( tBD );

                // compute the data
                mTurbine.compute() ;

                // kill if it makes no sense
                if ( mTurbine.error_code() != 0 )
                {
                    mDNA.kill() ;
                    return;
                }
                if( ! ( mTurbine.eta() < 1.0 ) )
                {
                    mDNA.kill() ;
                    return;
                }

                if( mTurbine.epsilon() < 0.1  )
                {
                    mDNA.kill() ;
                    return;
                }

                // compute the fitness
                real tFitness = 0.0 ;

                // just pull a tiny bit into best efficiency
                tFitness += std::pow( std::abs( mTurbine.eta() - 1.0 )*10, 2 );

                // we don't want a negative reaction
                if( mTurbine.reaction() < 0.01 )
                {
                    tFitness += std::pow( std::abs( mTurbine.reaction() -0.01 ) * 1000.0, 4 );
                }
                else if ( mTurbine.reaction() > 0.1 )
                {
                    tFitness += std::pow( std::abs( mTurbine.reaction() - 0.1 ) * 100.0, 2 );
                }

                // epsilon value must be reasonable
                if( mTurbine.epsilon() < 0.1 )
                {
                    tFitness += std::pow( std::abs( mTurbine.epsilon() - 0.1 ) * 100.0, 3 );
                }
                else if ( mTurbine.epsilon() > 1.0 )
                {
                    tFitness += std::pow( std::abs( mTurbine.epsilon() - 1.0 ) * 10.0, 3 );
                }

                // we don't want to cross the sonic area
                if( mTurbine.turbine_entry()->Ma() > 1.0 && mTurbine.turbine_discharge()->Ma() < 1.2 )
                {
                    tFitness += std::pow( std::abs( mTurbine.turbine_discharge()->Ma() - 1.2 ) * 10, 3 );
                }
                if( mTurbine.turbine_entry()->Ma() < 1.0 && mTurbine.turbine_discharge()->Ma() > 0.85 )
                {
                    tFitness += std::pow( std::abs(mTurbine.turbine_discharge()->Ma() - 0.85), 3 );
                }

                tFitness += std::pow( std::abs( mTurbine.blade_entry_error() ) * 100000, 5 );

                // std::cout << "compute " << this->phi() << " " << this->psi() << " " << this->bd() << " " << mTurbine.eta() << " " << tFitness << std::endl ;

                if( mTurbine.haller() < 0.8 )
                {
                    tFitness += std::pow( std::abs( mTurbine.haller() - 0.8 ) * 10, 2 );
                }

                mDNA.set_fitness( tFitness );
            }
        }

//------------------------------------------------------------------------------
    }
}