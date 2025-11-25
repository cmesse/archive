//
// Created by Christian Messe on 16.12.19.
//

#ifndef BELFEM_FN_BL_KAYS_CRAWFORD_HPP
#define BELFEM_FN_BL_KAYS_CRAWFORD_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace boundarylayer
    {
//------------------------------------------------------------------------------

        // see Kays, Crawford, Weigand:
        // Convective Heat and Mass Transfer
        // McGraw-Hill Education Ltd;  4th edition (2011)
        real
        kays_crawford(
                const real & aPr,
                const real & aMu,
                const real & aMuT,
                const real   aPrTinf=0.85  );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_BL_KAYS_CRAWFORD_HPP
