//
// Created by Christian Messe on 06.05.20.
//

#include "fn_BL_virtual_origin.hpp"
#include "assert.hpp"
#include "fn_BL_ReferenceTempertaure.hpp"
namespace belfem
{
    namespace boundarylayer
    {
        /**
         *
         * @param aState
         * @param x
         * @param aMode
         *      1: from laminar to turbulent on flat plate / cylinder
         *      2: from laminar to turbulent on cone
         *      3: from cone to cylinder, laminar
         *      4: from cone to cylinder, turbulent
         * @return
         */

        real virtual_origin( State & aState, const real & aX, uint aMode ) {
            // constants from Blasius solution
            const real Cl = 0.33205734;
            const real nl = 0.5 ;

            // Mangler factor for laminar case
            const real Mangler_l = std::sqrt( 3 );

            // constants for turbulent boundary layer, consistent with
            //  BL_flatplate_inc_turbulent
            const real Ct = 0.025077191459 ;
            const real nt = 0.186146870926 ;

            // Mangler factor for turbulent case
            const real Mangler_t = 1.176 ;

            // get gas
            Gas & tGas = aState.gas() ;

            real C0 ;
            real n0 ;
            bool turbulent0 ;

            real C ;
            real n ;
            bool turbulent ;

            // set factors
            switch ( aMode )
            {
                case ( 1 ) :
                {
                    // from laminar to turbulent on flat plate / cylinder
                    C0 = Cl;
                    n0 = nl;
                    turbulent0 = false;

                    C = Ct;
                    n = nt;
                    turbulent = true;
                    break;
                }
                case ( 2 ) :
                {
                    //  from laminar to turbulent on cone
                    C0 = Cl * Mangler_l ;
                    n0 = nl;
                    turbulent0 = false;

                    C = Ct * Mangler_t ;
                    n = nt;
                    turbulent = true;
                    break;
                }
                case ( 3 ):
                {
                    // from cone to cylinder, laminar
                    C0 = Cl * Mangler_l ;
                    n0 = nl;
                    turbulent0 = false;

                    C = Cl;
                    n = nl;
                    turbulent = false;
                    break;
                }
                case ( 4 ):
                {
                    // from cone to cylinder, turbulent
                    C0 = Ct * Mangler_t ;
                    n0 = nt;
                    turbulent0 = true;

                    C = Ct;
                    n = nt;
                    turbulent = true;
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "invalid state for virtual origin" );
                    return BELFEM_QUIET_NAN;
                }
            }

            // compute rigt hand side
            real T = reference_temperature(
                    tGas,
                    aState.T(),
                    aState.p(),
                    aState.u(),
                    aState.Tw(),
                    turbulent0 );

            const real & p = aState.p() ;
            const real & u = aState.u() ;

            real rho = tGas.rho( T, p );
            real mu = tGas.mu( T, p );

            real Re_x = aState.rho() * aState.u()  * aX / aState.mu() ;

            real delta2 = C0 * std::pow( aX, 1.0 - n0 ) / ( 1.0 - n0 )
                          * std::pow( rho * mu / ( aState.mu() * aState.rho() * Re_x ) , n0 )
                          * std::pow( rho / aState.rho(), 1.0 - 2.0 * n0 ) ;

            real RHS = delta2 ;
            real  x = aX ;

            // update reference temperature
            if( turbulent != turbulent0 )
            {
                T = reference_temperature(
                        tGas,
                        aState.T(),
                        aState.p(),
                        aState.u(),
                        aState.Tw(),
                        turbulent );

                rho = tGas.rho( T, p );
                mu = tGas.mu( T, p );
            }

            // help constants
            real K1 = std::pow( rho * mu / ( aState.mu() * aState.rho() ) , n )
                      * std::pow( rho / aState.rho(), 1.0 - 2.0 * n ) ;

            real K2 = mu / ( rho * u );

            real K3 ;

            real f = 1.0 ;
            real df = 1.0 ;

            // delta2 = C * (x^(1-n))/(1-n) * K1 * ( K2/x)^n
            uint tCount = 0;
            real tOmega ;
            real tOmega0 = 1.0 ;
            real tOmega1 ;
            while ( ( abs( f ) > 1e-7 ) && abs( df ) > 1e-7 )
            {
                K3 = std::pow( K2 / x, n );

                delta2 = C * std::pow( x, 1.0 - n ) / ( 1.0 - n ) * K1 * K3;
                f = delta2 - RHS;
                df = ( C * K1 * std::pow( x, -n ) * ( n * 2.0 - 1.0 ) * K3 ) / ( n - 1.0 );

                tOmega1 = std::abs( 0.9 * ( x - 1e-7 ) * df/f );
                tOmega = tOmega1 < tOmega0 ? tOmega1 : tOmega0 ;

                x -= tOmega * f / df ;


                BELFEM_ERROR( tCount++ < 1000,
                    "too many iterations while trying to find virtual origin" );
            }

            return x ;
        }

    }
}