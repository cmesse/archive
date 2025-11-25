//
// Created by Christian Messe on 24.09.19.
//
#include "fn_CN_arrhenius.hpp"
#include "constants.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        real
        arrhenius( const Vector< real > & aCoeffs, const real & aT)
        {
            return aCoeffs( 0 ) * std::pow( aT, aCoeffs( 1 ) ) *
                    std::exp( -  aCoeffs( 2 ) / ( constant::Rm_cal * aT ) );
        }

//------------------------------------------------------------------------------

        real
        darrheniusdT(
                const Vector< real > & aCoeffs,
                const real & aT,
                const real & aArrhenius )
        {
            return aArrhenius * ( aCoeffs( 2 ) + aCoeffs( 1 ) *  constant::Rm_cal * aT )
                / ( constant::Rm_cal * aT * aT );
        }

//------------------------------------------------------------------------------

        real
        arrhenius_simple( const Vector< real > & aCoeffs, const real & aT)
        {
            BELFEM_ASSERT( aCoeffs( 1 ) == 0.0, "Forbidden arrhenius function called" );

            return aCoeffs( 0 ) * std::exp( -  aCoeffs( 2 ) / ( constant::Rm_cal * aT ) );
        }

//------------------------------------------------------------------------------

        real
        darrhenius_simpledT( const Vector< real > & aCoeffs, const real & aT, const real & aArrhenius )
        {
            BELFEM_ASSERT( aCoeffs( 1 ) == 0.0, "Forbidden arrhenius function derivative called" );

            return aArrhenius * aCoeffs( 2 ) / ( constant::Rm_cal * aT * aT );
        }

//------------------------------------------------------------------------------
    }
}