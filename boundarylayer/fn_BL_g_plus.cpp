//
// Created by Christian Messe on 2019-03-30.
//
#include "fn_BL_g_plus.hpp"

namespace belfem
{
    namespace boundarylayer
    {
        // see 10.2514/6.2017-4743
//------------------------------------------------------------------------------

        real
        g_plus( const real & aKappa, const real & aPi, const real & aEta )
        {
            return aEta * aEta * ( ( 6.0 * aPi + 1.0 )
                                 - ( 4.0 * aPi + 1.0 ) * aEta ) / aKappa;
        }

//------------------------------------------------------------------------------

        real
        dg_plus_deta( const real & aKappa, const real & aPi, const real & aEta )
        {
            return aEta * ( ( 12.0 * aPi + 2.0 )
                   - ( 12.0 * aPi + 3.0 ) * aEta ) / aKappa;
        }

//------------------------------------------------------------------------------
    }
}