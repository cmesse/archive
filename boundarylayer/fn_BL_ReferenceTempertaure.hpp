//
// Created by Christian Messe on 30.11.19.
//

#ifndef BELFEM_FN_BL_REFERENCETEMPERTAURE_HPP
#define BELFEM_FN_BL_REFERENCETEMPERTAURE_HPP

#include "typedefs.hpp"
#include "cl_Gas.hpp"

namespace belfem
{
    namespace boundarylayer
    {
        /**
         * computes the reference temperature
         */
         real
         reference_temperature(
                 Gas        & aGas,
                 const real & aTe,
                 const real & aPe,
                 const real & aUe,
                 const real & aTw,
                 const bool   aIsTurbulent );

    }
}
#endif //BELFEM_FN_BL_REFERENCETEMPERTAURE_HPP
