//
// Created by Christian Messe on 19.04.20.
//
#include "typedefs.hpp"
#include "assert.hpp"
#include "fn_BL_VanDriest.hpp"
#include "fn_BL_Eckert.hpp"
#include "fn_BL_ReferenceTempertaure.hpp"
#include "fn_BL_cf_flatplate_inc_turbulent.hpp"
namespace belfem
{
    namespace boundarylayer
    {
//------------------------------------------------------------------------------

        void
        vandriest( State & aState, const real & aX, const real aKs, const real aMangler )
        {
            // uncertainty for friction
            real tBias = 1.0e-5 ;

            // turbulent prandtl number
            const real tPrT = 0.85 ;

            // Kármán constant
            const real tKarman = 0.41 ;

            // Defect Parameter
            const real tFe = 3.78 ;

            // parameter for law of the wall
            real tBplus = 5.0 ;

            // relaxation factor
            real tOmega = 0.9 ;

            // flag that tells if we want to compute the recovery factor
            // ( we don't do this for a cone )
            bool tComputeSigmaRecov = aMangler == 1.0 ;

            // pressure in boundary layer
            const real & tP = aState.p() ;

            // velocity
            const real & tU = aState.u() ;

            // freestream enthalpy
            const real & tH = aState.h() ;

            // enthalpy at wall
            const real & tHw = aState.hw() ;

            // help values for Reynolds-Colburn computation

            const real tXi   = aState.Pr() / tPrT - 1.0 ;
            const real tEta  = std::log( 1.0 + 5.0 / 6.0 * tXi );
            const real tZeta = std::log( 6.0 ) * std::log( ( 1.0  + 0.875 * tXi) / ( 1.0 + 0.25 * tXi ) );
            const real tBeta = constant::pi * constant::pi / 6.0 + 1.5 * ( 1.0 - tPrT ) ;

            // get pointer to gas
            Gas & tGas = aState.gas();

            // reference temperature for the state
            real tTref = eckert( aState, aX, true, aMangler );

            // initial value for recovery factor
            real tRecovery = std::pow( tGas.Pr( tTref, tP ), 1.0/3.0 ) ;

            // initial value for Reynolds-Colburn analogy
            real tSigma = tRecovery * tRecovery ;

            // recovery enthalpy
            real tHr = tH + 0.5 * tRecovery * tU * tU;

            // freestream reynolds number
            real tReX = aState.rho() * tU * aX / aState.mu() ;


            // initial value to start the iteration loop
            real tCfreset = aState.tau_w() / ( 0.5 * aState.rho() * tU * tU ) ;
            real tCf  = tCfreset;
            real tCf0 = 0.0 ;
            real tCf1 ;

            // boundary layer thickness
            real tDelta ;

            // rotta clauser parameter
            real tRottaClauser  = 0.0 ;

            real tUtau = 0.0;

            // Wake parameter
            real tPi = 0.55 ;

            real tTauW;

            bool tCompressible = aState.Ma() > 0.2 ;

            // initialize counter
            uint tCount = 0 ;

            if( tCompressible )
            {
                const cplx tOne( 1.0, 0.0 );

                // factors for Crocco-Busemann transformation
                cplx tPhi ;
                cplx tPsi ;
                cplx tChi ;
                cplx tA ;
                cplx tB ;

                // Friction conversion factor
                real tFcf;

                // Reymolds conversion factor
                real tFReX;


                real tError = 1.0 ;
                while( tError > tBias && tCount++ < 1000 )
                {
                    // catch error if heat flux is zero
                    if( std::abs( tHr - tHw ) < 1e-6 )
                    {
                        // use Eckert relations for calculation
                        tFcf  = aState.rho() / tGas.rho( tTref, tP );
                        tFReX = aState.mu() / ( tFcf * tGas.mu( tTref, tP ) );

                        // friction factor
                        tCf = cf_flatplate_inc_turbulent( tFReX * tReX, tBplus, tKarman ) * aMangler / tFcf ;

                        // compute shear stress
                        tTauW = 0.5 * tCf * aState.rho() * tU * tU ;

                        // write values into state
                        aState.set_loads( tTauW, 0.0, tHw );

                        // exit this function
                        return ;
                    }
                    else if( tCount == 100 )  // catch oscillation
                    {
                        tOmega = 0.05 ;
                        tBias = 0.01;

                        // reset values
                        tCf0 = tCfreset;
                        tCf = tCfreset ;
                        tBplus = 5.0 ;
                    }
                    else
                    {
                        // Eq. ( 21 ) doi: 10.2514/6.2017-4743
                        tPsi = aState.mu_w() * ( tHr - aState.hw() ) /
                               ( tSigma * aState.lambda_w() * aState.Tw() );

                        // Eq. ( 20 ) doi: 10.2514/6.2017-4743
                        tPhi = tPsi + ( 1.0 - aState.T() / aState.Tw() ) * tOne ;

                        tChi = std::sqrt( tPsi * tPsi + 4.0 * tPhi );

                        // Factors for Crocco-Busemann Transformation, see doi: 10.2514/6.2017-4743, Eq. ( 22 )
                        tA = std::asin( ( 2.0 * tPhi - tPsi ) / tChi );
                        tB = std::asin( tPsi / tChi );

                        //  Driest Transformation Functions, see Driest but using
                        // enthalpies instead of temperatures
                        tFcf = ( tHr / tH - 1.0 ) / std::real( std::pow( tA + tB, 2 ) );

                        // White and Christoph, Eq. ( 25 ), but using enthalpies
                        tFReX = aState.mu() / aState.mu_w() * std::sqrt( tH / ( tHw * tFcf ) );
                        //tFReX = aState.mu() / ( tFcf * aState.mu_w()  );

                        // shift value
                        tCf0 = tCf ;

                        // compute new value
                        tCf1 = tOmega * cf_flatplate_inc_turbulent( tReX * tFReX, tBplus, tKarman, tPi, tFe ) * aMangler / tFcf ;

                        // sanity check
                        if ( ! ( tCf1 > 0.0 ) )
                        {
                            // fallback to incompressible case
                            tCompressible = false ;

                            // break the loop
                            break ;
                        }

                        // relax value
                        tCf = tCf0 * ( 1.0 - tOmega ) + tCf1 * tOmega ;

                        // the mangler transformation is a simplified version of this
                        // we don't want to do both at the same time
                        if( aMangler == 1.0 )
                        {
                            // boundary layer thickness
                            tDelta = 0.37 * aX / std::pow( tFReX * tReX, 0.2 );

                            // shear velocity
                            tUtau = tU * std::sqrt( 0.5 * tCf );

                            // Rotta-Clauser parameter, Schlichting (18.85)
                            tRottaClauser = -tDelta / tUtau * aState.dudx() * tFe;

                            // correlation representing Fig. 18.2
                            tPi = ( -2.65920E-02 * tRottaClauser + 6.87076E-01 ) * tRottaClauser + 0.55;
                        }
                    }

                    // Calculate Reynolds-Colburn Analogy and Recovery Factor
                    if( tComputeSigmaRecov )
                    {
                        real tAlpha = std::sqrt( 0.5 * tCf );

                        // Driest 1954, Eq. ( 41 )
                        tSigma = tPrT * ( 1.0 + 5.0 * tAlpha
                                * ( 0.2 / tKarman * ( 1.0 - tPrT ) * tBeta + tXi + tEta ) );

                        // Driest 1954, Eq. ( 54 )
                        tRecovery = tPrT * ( 1.0 + 2.0/tKarman * tAlpha
                             * ( 1.0- tPrT ) * tBeta + 12.5 * tCf
                              *( ( tXi + 2.0 * tEta + tZeta ) ) );


                        // sanity check
                        if( ! ( tRecovery > 0.0 ) )
                        {
                            // initial value for recovery factor
                            tRecovery = std::pow( tGas.Pr( tTref, tP ), 1.0/3.0 ) ;

                            // initial value for Reynolds-Colburn analogy
                            tSigma = tRecovery * tRecovery ;

                            // deactivate auto-computation for this panel
                            tComputeSigmaRecov = false ;
                        }

                        // compute new recovery enthalpy
                        tHr = tH + 0.5 * tRecovery * tU * tU ;

                    }

                    // compute friction
                    tTauW = 0.5 * tCf * aState.rho() * aState.u() * aState.u() ;

                    if( aKs > 0.0 )
                    {
                        // compute shear velocity
                        tUtau = std::sqrt( tTauW / aState.rho_w() );

                        // Roughness Factor Schlichting ( 17.31 )
                        real tKsPlus =  aState.rho_w() * tUtau * aKs / aState.mu_w() * 0.001 ;

                        //  Offset for Wall Function, Schlichting ( 17.40 )
                        tBplus = 5.0 - std::log( 1.0 + tKsPlus / 3.4 ) / tKarman ;
                    }

                    // compute bias
                    tError = std::abs( tCf - tCf0 ) / tCf0 ;
                }
            }

            // incompressible case
            if( ! tCompressible )
            {
                // reset values
                tBplus = 5.0 ;
                tCf  = tCfreset;
                tCf0 = 0.0 ;

                // reset counter
                tCount = 0 ;

                real tError = 1.0 ;
                while( tError > tBias && tCount++ < 1000 )
                {

                    // shift value
                    tCf0 = tCf ;

                    // compute new value
                    tCf1 = tOmega * cf_flatplate_inc_turbulent( tReX, tBplus, tKarman, tPi, tFe ) * aMangler ;

                    // relax
                    tCf = tCf0 * ( 1.0 - tOmega ) + tOmega * tCf1 ;


                    // boundary layer thickness
                    tDelta = 0.37 * aX / std::pow( tReX, 0.2 );
                    if( aMangler == 1 )
                    {
                        // shear velocity
                        tUtau = tU * std::sqrt( 0.5 * tCf );

                        // Rotta-Clauser parameter, Schlichting (18.85)
                        tRottaClauser = -tDelta / tUtau * aState.dudx() * tFe;

                        // correlation representing Fig. 18.2
                        tPi = ( -2.65920E-02 * tRottaClauser + 6.87076E-01 ) * tRottaClauser + 0.55;
                    }

                    // shear stress
                    tTauW = 0.5 * tCf * aState.rho() * aState.u() * aState.u() ;

                    if( aKs > 0.0 )
                    {
                        // compute shear velocity
                        tUtau = std::sqrt( tTauW / aState.rho_w() );

                        // Roughness Factor Schlichting ( 17.31 )
                        real tKsPlus =  aState.rho_w() * tUtau * aKs / aState.mu_w() * 0.001 ;

                        //  Offset for Wall Function, Schlichting ( 17.40 )
                        tBplus = 5.0 - std::log( 1.0 + tKsPlus / 3.4 ) / tKarman ;

                    }

                    // Calculate Reynolds-Colburn Analogy and Recovery Factor
                    if( tComputeSigmaRecov )
                    {
                        real tAlpha = std::sqrt( 0.5 * tCf );

                        // Driest 1954, Eq. ( 41 )
                        tSigma = tPrT * ( 1.0 + 5.0 * tAlpha
                                                * ( 0.2 / tKarman * ( 1.0 - tPrT ) * tBeta + tXi + tEta ) );

                        // Driest 1954, Eq. ( 54 )
                        tRecovery = tPrT * ( 1.0 + 2.0/tKarman * tAlpha
                                    * ( 1.0- tPrT ) * tBeta + 12.5 * tCf
                                    *( ( tXi + 2.0 * tEta + tZeta ) ) );

                        // sanity check
                        if( ! ( tRecovery > 0.0 ) )
                        {
                            // initial value for recovery factor
                            tRecovery = std::pow( tGas.Pr( tTref, tP ), 1.0/3.0 ) ;

                            // initial value for Reynolds-Colburn analogy
                            tSigma = tRecovery * tRecovery ;

                            // deactivate auto-computation for this panel
                            tComputeSigmaRecov = false ;
                        }

                        // compute new recovery enthalpy
                        tHr = tH + 0.5 * tRecovery * tU * tU ;

                    }

                    // compute bias
                    tError = std::abs( tCf - tCf0 ) / tCf0 ;
                }
            }

            // sanity check
            BELFEM_ERROR( tCount < 1000,
                    "Too many iterations for Re_x = %E, Te = %E, pe=%E, Ma=%E, Tw=%E",
                    tReX, aState.T(), aState.p(), aState.Ma(), aState.Tw() );

            // finalize
            real tDotQ = tTauW / ( tSigma * aState.u() ) * ( tHr - tHw ) ;

            // write result into state
            aState.set_loads( tTauW, tDotQ, tHr );
        }

//------------------------------------------------------------------------------
    }
}