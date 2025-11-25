//
// Created by Christian Messe on 16.12.19.
//
#include <iostream>
#include "fn_BL_Kays_Crawford.hpp"

namespace belfem
{
    namespace boundarylayer
    {
//------------------------------------------------------------------------------
        real
        kays_crawford(
                const real & aPr,
                const real & aMu,
                const real & aMuT,
                const real   aPrTinf  )
        {
            // help constant
            const real tC = 0.3 * aMuT / aMu * aPr ;

            // Equation ( 12-7 ) ( note typo in book! )
            return 1.0 / ( 0.5 / aPrTinf + tC * ( 1.0/ sqrt( aPrTinf )
                   - tC * ( 1.0 - std::exp(
                   - 1.0 / ( tC * std::sqrt( aPrTinf ) ) ) ) ) );
        }

//------------------------------------------------------------------------------
    }
}