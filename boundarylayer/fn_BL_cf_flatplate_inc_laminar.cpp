//
// Created by Christian Messe on 15.04.20.
//

#include "fn_BL_cf_flatplate_inc_laminar.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    namespace boundarylayer
    {
        real
        cf_flatplate_inc_laminar( const real & aReX )
        {
            // number from Blasius solution
            return 0.66411468 / std::sqrt( aReX );
        }
    }

//------------------------------------------------------------------------------
}