//
// Created by Christian Messe on 2019-03-30.
//

#ifndef BELFEM_FN_BL_G_PLUS_HPP
#define BELFEM_FN_BL_G_PLUS_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace boundarylayer
    {
        // see 10.2514/6.2017-4743
//------------------------------------------------------------------------------

        real
        g_plus( const real & aKappa, const real & aPi, const real & aEta );

//------------------------------------------------------------------------------

        real
        dg_plus_deta( const real & aKappa, const real & aPi, const real & aEta );

//------------------------------------------------------------------------------
    }
}
#endif // BELFEM_FN_BL_G_PLUS_HPP
