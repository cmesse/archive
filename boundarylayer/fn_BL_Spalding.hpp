//
// Created by Christian Messe on 2019-03-30.
//

#ifndef BELFEM_FN_BL_SPALDING_HPP
#define BELFEM_FN_BL_SPALDING_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace boundarylayer
    {
        // see 10.1115/1.3641728.
//------------------------------------------------------------------------------

        real
        spalding( const real & aB,
                  const real & aKappa,
                  const real & aExp, // factor e^(-kappa * B )
                  const real & aYplus,
                  const real   aF_guess=0.0 );
//------------------------------------------------------------------------------

        real
        spalding_y(
                const real & aB,
                const real & aKappa,
                const real & aExp,   // factor e^(-kappa * B )
                const real & aF );

//------------------------------------------------------------------------------

        real
        spalding_dydf(
                const real & aB,
                const real & aKappa,
                const real & aExp,   // factor e^(-kappa * B )
                const real & aF );

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_FN_BL_SPALDING_HPP
