//
// Created by Christian Messe on 24.09.19.
//

#ifndef BELFEM_FN_CN_ARRHENIUS_HPP
#define BELFEM_FN_CN_ARRHENIUS_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        /**
         * Arrhenius function
         * \f$ A \, T^b \, \exp{\left( -\frac{E_a}{R_m \, T} \right)} \f$
         */
        real
        arrhenius( const Vector< real > & aCoeffs, const real & aT);

//------------------------------------------------------------------------------

        /**
         * temperature derivative of Arrhenius function
         *
         * \f$ \frac{ A \, \left( E_a + b \, R_m \, T\right)}{R_m \, T^{2-b}} \,\exp{\left( -\frac{E_a}{R_m \, T} \right)} \f$
         */
        real
        darrheniusdT(
                const Vector< real > & aCoeffs,
                const         real & aT,
                const         real & aArrhenius  );

//------------------------------------------------------------------------------

        /**
         * simplified Arrhenius funcition with b=0
         *
         * \f$ A \,\exp{\left( -\frac{E_a}{R_m \, T } \right)} \f$
         */
        real
        arrhenius_simple( const Vector< real > & aCoeffs, const real & aT );

//------------------------------------------------------------------------------
        /**
         * temperatue derivative of simplified Arrhenius funcition with b=0
         *
         * \f$ \frac{ A \, E_a}{R_m \, T^2} \,\exp{\left( -\frac{E_a}{R_m \, T } \right)} \f$
         */
        real
        darrhenius_simpledT(
                const Vector< real > & aCoeffs,
                const real & aT,
                const real & aArrhenius);

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_CN_ARRHENIUS_HPP
