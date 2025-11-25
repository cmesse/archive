//
// Created by Christian Messe on 15.04.20.
//

#include "fn_BL_Eckert.hpp"
#include "fn_BL_ReferenceTempertaure.hpp"
#include "fn_BL_cf_flatplate_inc_laminar.hpp"
#include "fn_BL_cf_flatplate_inc_turbulent.hpp"
namespace belfem
{
    namespace boundarylayer
    {
        real
        eckert( State & aState, const real & aX, bool aIsTurbulent, const real aMangler )
        {
            // get the gas of the state
            Gas & tGas = aState.gas();

            // reference temperature for the state
            real aTref = reference_temperature(
                   tGas,
                   aState.T(),
                   aState.p(),
                   aState.u(),
                   aState.Tw(),
                   aIsTurbulent );

            // reference enthalpy
            // real tHref = tGas.h( aTref, aState.p() );

            // compute density
            real tRho = tGas.rho( aTref, aState.p() );

            // viscosity
            real tMu  = tGas.mu( aTref, aState.p() );

            // incompressible Reynolds number
            real tReX_inc = tRho * aState.u() * aX / tMu ;

            // Prandtl  number
            real tPr = tGas.Pr( aTref, aState.p() );

            // recovery factor
            real tRecovery ;

            // Reynolds-Colburn transformation
            real tSigma = std::pow( tPr, 2.0/3.0 );

            // incompressible friction factor
            real tCf_inc ;

            //set_loads
            if( aIsTurbulent )
            {
                tRecovery = std::pow( tPr, 1.0/3.0 );
                tCf_inc = cf_flatplate_inc_turbulent( tReX_inc ) ;
                //tCf_inc = 2.0 * 0.0296 / std::pow( tReX_inc, 0.2 );
            }
            else
            {
                tRecovery = std::sqrt( tPr );
                tCf_inc = cf_flatplate_inc_laminar( tReX_inc );
            }

            // compute actual friction factor
            real tCf = tCf_inc * aMangler * tRho / aState.rho() ;

            // shear stress
            real tTauw = 0.5 * tCf * aState.rho() * aState.u() * aState.u() ;

            // recovery enthalpy
            real tHr =  aState.h() + 0.5 * tRecovery * aState.u() * aState.u() ;

            // heat load
            real tDotQ = tTauw * ( tHr - aState.h() ) / ( tSigma * aState.u() );

            // write values in state
            aState.set_loads( tTauw, tDotQ, tHr );

            // return reference temperature
            return aTref ;
        }
    }
}