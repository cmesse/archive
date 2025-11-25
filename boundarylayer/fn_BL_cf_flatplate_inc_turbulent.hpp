//
// Created by Christian Messe on 15.04.20.
//

#ifndef BELFEM_FN_BL_CF_FLATPLATE_INC_TURBULENT_HPP
#define BELFEM_FN_BL_CF_FLATPLATE_INC_TURBULENT_HPP

#include "typedefs.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace boundarylayer
    {
        /**
         * see Schlichting: Grenzschichttheorie, 10th Edition, Springer 2006
         * @param aReX      incompressible Reynolds number
         * @param aBplus    Mixing integration offset
         * @param aKarman   Kármán constant
         * @param aPiWake   Coles Wake Parameter
         * @param aFe       Defect at edge of boundary layer
         * @return
         */
        real
        cf_flatplate_inc_turbulent(
                const real & aReX,
                const real aBplus=5.0,
                const real aKarman=0.41,
                const real aPiWake=0.55,
                const real aFe=3.78 );

    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_BL_CF_FLATPLATE_INC_TURBULENT_HPP
