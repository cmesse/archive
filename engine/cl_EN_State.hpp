//
// Created by Christian Messe on 31.08.20.
//

#ifndef BELFEM_CL_EN_STATE_HPP
#define BELFEM_CL_EN_STATE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Gas.hpp"

// geometry
#define BELFEM_ENGINE_STATE_X       0
#define BELFEM_ENGINE_STATE_A       1
#define BELFEM_ENGINE_STATE_DH      2

// thermodynamic state
#define BELFEM_ENGINE_STATE_T       3
#define BELFEM_ENGINE_STATE_P       4
#define BELFEM_ENGINE_STATE_U       5
#define BELFEM_ENGINE_STATE_RHO     6

// mixture
#define BELFEM_ENGINE_STATE_M       7
#define BELFEM_ENGINE_STATE_R       8

// caloric
#define BELFEM_ENGINE_STATE_H       9
#define BELFEM_ENGINE_STATE_CP     10
#define BELFEM_ENGINE_STATE_S      11
#define BELFEM_ENGINE_STATE_GAMMA  12
#define BELFEM_ENGINE_STATE_W      13

// total
#define BELFEM_ENGINE_STATE_TT     14
#define BELFEM_ENGINE_STATE_PT     15
#define BELFEM_ENGINE_STATE_HT     16
// transport
#define BELFEM_ENGINE_STATE_MU     17
#define BELFEM_ENGINE_STATE_LAMBDA 18

// similarity
#define BELFEM_ENGINE_STATE_PR     19
#define BELFEM_ENGINE_STATE_MA     20
#define BELFEM_ENGINE_STATE_RE     21

#define BELFEM_ENGINE_NUMSTATES    22

namespace belfem
{
    namespace engine
    {
        class Analysis ;
//------------------------------------------------------------------------------

        class State
        {
            Gas & mCombgas ;

            const string mLabel ;

            // state values according to BELFEM_ENGINE_***
            Vector< real > mValues ;

            // mass fractions of gas
            Vector< real > mMassFractions ;

            // molar fractions
            Vector< real > mMolarFractions ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            State( Gas & aCombgas, const string & aLabel, const uint aNumSpecies=0 );

//------------------------------------------------------------------------------

            ~State() = default ;

//------------------------------------------------------------------------------

            /**
             * expose the values container
             */
             Vector< real > &
             values() ;

//------------------------------------------------------------------------------

            /**
             * return the temperature value
             */
            const real &
            T() const ;

//------------------------------------------------------------------------------

            /**
             * return the pressure value
             */
            const real &
            p() const ;
//------------------------------------------------------------------------------

            /**
             * return the total temperature value
             */
            const real &
            Tt() const ;

//------------------------------------------------------------------------------

            /**
             * return the total pressure value
             */
            const real &
            pt() const ;

//------------------------------------------------------------------------------

            /**
             * total enthalpy
             */
            const real &
            ht() const ;

//------------------------------------------------------------------------------

            /**
             * return the density
             */
            const real &
            rho() const ;

//------------------------------------------------------------------------------

            /**
             * return the ratio of specific heats
             */
            const real &
            gamma() const ;

//------------------------------------------------------------------------------

            /**
             * return the enthalpy
             */
            const real &
            h() const ;

//------------------------------------------------------------------------------

            /**
             * return the entropy
             */
            const real &

            s() const ;
//------------------------------------------------------------------------------

            /**
             * return the velocity
             */
            const real &
            u() const ;

//------------------------------------------------------------------------------

            /**
             * return the Mach number
             */
            const real &
            Ma() const ;

//------------------------------------------------------------------------------

            /**
             * return the speed of sound
             */
            const real &
            w() const ;

//------------------------------------------------------------------------------

            /**
             * return the cross section
             */
            const real &
            A() const ;

//------------------------------------------------------------------------------

            /**
             * expose the mass fractions
             */
            Vector< real > &
            mass_fractions() ;

            const Vector< real > &
            mass_fractions() const ;

//------------------------------------------------------------------------------

            /**
             * expose the molar fractions
             */
            Vector< real > &
            molar_fractions() ;

            const Vector< real > &
            molar_fractions() const ;

//------------------------------------------------------------------------------

            real &
            value( const index_t aIndex );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            const real &
            value( const index_t aIndex ) const ;

//------------------------------------------------------------------------------

            void
            compute_caloric( const real & aT, const real & aP, const real aU=0.0 );

//------------------------------------------------------------------------------

            void
            compute_equilibrium(
                    const real & aT,
                    const real & aP,
                    const real & aH );

//------------------------------------------------------------------------------

            void
            print() const ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline const real &
        State::T() const
        {
            return mValues( BELFEM_ENGINE_STATE_T );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::p() const
        {
            return mValues( BELFEM_ENGINE_STATE_P );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::Tt() const
        {
            return mValues( BELFEM_ENGINE_STATE_TT );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::pt() const
        {
            return mValues( BELFEM_ENGINE_STATE_PT );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::ht() const
        {
            return mValues( BELFEM_ENGINE_STATE_HT );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::rho() const
        {
            return mValues( BELFEM_ENGINE_STATE_RHO );
        }

//------------------------------------------------------------------------------


        inline const real &
        State::gamma() const
        {
            return mValues( BELFEM_ENGINE_STATE_GAMMA );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::h() const
        {
            return mValues( BELFEM_ENGINE_STATE_H );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::s() const
        {
            return mValues( BELFEM_ENGINE_STATE_S );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::u() const
        {
            return mValues( BELFEM_ENGINE_STATE_U );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::Ma() const
        {
            return mValues( BELFEM_ENGINE_STATE_MA );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::w() const
        {
            return mValues( BELFEM_ENGINE_STATE_W );
        }

//------------------------------------------------------------------------------

        inline const real &
        State::A() const
        {
            return mValues( BELFEM_ENGINE_STATE_A );
        }
//------------------------------------------------------------------------------

        inline Vector< real > &
        State::values()
        {
            return mValues ;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        State::mass_fractions()
        {
            return mMassFractions ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        State::mass_fractions() const
        {
            return mMassFractions ;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        State::molar_fractions()
        {
            return mMolarFractions ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        State::molar_fractions() const
        {
            return mMolarFractions ;
        }

//------------------------------------------------------------------------------

        inline real &
        State::value( const index_t aIndex )
        {
            return mValues( aIndex );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline const real &
        State::value( const index_t aIndex ) const
        {
            return mValues( aIndex );
        }


    }
}
#endif //BELFEM_CL_EN_STATE_HPP
