//
// Created by Christian Messe on 30.11.19.
//

#include "fn_BL_ReferenceTempertaure.hpp"

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
                const bool   aIsTurbulent )
        {

            // approximate total temperature
            real tTt = aTe + 0.5 * aUe * aUe / aGas.cp( aTe, aPe );

            // initial guess for Te using Eckerts approximation
            real aTref = aTe + 0.5 * ( aTw - aTe ) + 0.22 * ( tTt - aTe );

            // relaxation
            real tOmega = 0.9;

            // Prandtl Number
            real tPrandtl;

            // recovery factor
            real tRecovery;

            // stagnation enthalpy
            real tHe = aGas.h( aTe, aPe );

            // wall enthalpy
            real tHw = aGas.h( aTw, aPe );

            // recovery enthalpy
            real tHr;

            // reference enthalpy
            real tHref = tHe;

            real tHref_old;

            real tResidual = BELFEM_INT_MAX;

            // iteration counter
            uint tCount = 0;

            real tCe;
            real tCr;
            real tCw;
            real tPow;

            // from Meador and Smart, 10.2514/1.2656
            if( aIsTurbulent )
            {
                tCe  = 0.34;
                tCr  = 0.16;
                tCw  = 0.50;
                tPow = 1.0/3.0;
            }
            else
            {
                tCe  = 0.29;
                tCr  = 0.16;
                tCw  = 0.55;
                tPow = 0.50;
            }

            while( tResidual > 1e-9 )
            {
                // compute prandtl number
                tPrandtl = aGas.Pr( aTref, aPe );

                // compute recovery factor
                tRecovery = std::pow( tPrandtl, tPow );

                // compute recovery enthalpy
                tHr = tHe + 0.5 * tRecovery * aUe * aUe;

                // shift enthalpy
                tHref_old = tHref;

                // compute reference enthalpy using values
                // from Meador and Smart, 10.2514/1.2656
                tHref = tCe * tHe + tCr * tHr + tCw * tHw;

                // compute residual
                tResidual = std::abs(
                        ( tHref - tHref_old ) / tHref_old );

                aTref *= ( 1.0 - tOmega );

                // compute temperature
                aTref += tOmega * aGas.T_from_h( tHref, aPe );

                BELFEM_ERROR( tCount++ < 1000,
                             "Too many iterations ( T: %12.6f, p: %12.6f, u: %12.6f, Tw: %12.6f, Res: %8.3e",
                             ( double ) aTe,
                             ( double ) aPe,
                             ( double ) aUe,
                             ( double ) aTw,
                             ( double ) tResidual );

                if( tCount == 50 )
                {
                    tOmega = 0.1 ;
                }
            }

            return aTref;
        }
    }
}