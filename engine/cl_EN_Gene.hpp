//
// Created by Christian Messe on 02.05.21.
//

#ifndef BELFEM_CL_EN_GENE_HPP
#define BELFEM_CL_EN_GENE_HPP

#include "typedefs.hpp"
#include "cl_DNA_old.hpp"
#include "cl_EN_Turbine.hpp"
namespace belfem
{
    namespace engine
    {
        class Gene
        {
            Turbine & mTurbine ;

            // genome
            belfem::DNA< 3 > mDNA ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

             Gene( Turbine & aTurbine ) ;

//------------------------------------------------------------------------------

            ~Gene() = default ;

//------------------------------------------------------------------------------

            /**
             * reset the bitset
             */
            void
            reset();

//------------------------------------------------------------------------------

            void
            set_values( const real & aPhi, const real & aPsi, const real & aBD );

//------------------------------------------------------------------------------

            void
            set_values( const Vector< real > & aValues );

//------------------------------------------------------------------------------

            /**
             * inherit data from parents
             */
             void
             inherit( const Gene * aMom, const Gene * aDad );

//------------------------------------------------------------------------------

            void
            compute();

//------------------------------------------------------------------------------

            /**
             * return the phi value
             */
            inline real
            phi() const
            {
                return mDNA.get_value( 0 );
            }

//------------------------------------------------------------------------------

            /**
             * return the phi value
             */
            inline real
            psi() const
            {
                return mDNA.get_value( 1 );
            }

//------------------------------------------------------------------------------

            /**
             * return the phi value
             */
            inline real
            bd() const
            {
                return mDNA.get_value( 2 );
            }

//------------------------------------------------------------------------------

            /**
             * test a genome
             */
             inline bool test( const index_t & aIndex ) const
             {
                return mDNA.test( aIndex );
             }

//------------------------------------------------------------------------------

            /**
             * return the fitness
             */
            inline const real &
            fitness() const
            {
                return mDNA.fitness() ;
            }

//------------------------------------------------------------------------------

            // return the alive state
            inline bool
            alive() const
            {
                return mDNA.alive() ;
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

        // comparision object
        struct
        {
            inline bool
            operator()( const Gene * aA, const Gene * aB )
            {
                return aA->fitness() < aB->fitness();
            }
        } opCompareFitness ;

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif


//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_EN_GENE_HPP
