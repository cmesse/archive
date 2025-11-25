//
// Created by Christian Messe on 31.08.20.
//

#include "cl_EN_Analysis.hpp"
#include "assert.hpp"
#include "fn_gesv.hpp"
#include "cl_Matrix.hpp"
#include "fn_linspace.hpp"
#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "cl_GT_RefGas.hpp"
#include "cl_Timer.hpp"

namespace belfem
{
    namespace engine
    {

//------------------------------------------------------------------------------

        Analysis::Analysis( const Parameters & aParams ) :
            mParams( aParams ),
            mCombgas( aParams.create_gas() ),
            mInjector( *mCombgas, "Injector", mParams.number_of_species() ),
            mTotal( *mCombgas, "Total",mParams.number_of_species()  ),
            mThroat( *mCombgas, "Throat",mParams.number_of_species() ),
            mNozzle( *mCombgas, "Nozzle",mParams.number_of_species() )
        {
            // initialize work matrices
            mJ2.set_size( 2,2 );
            mF2.set_size( 2 );
            mJ3.set_size( 3,3 );
            mF3.set_size( 3 );
            mPivot.set_size( 3 );

            // get the cross section at the throat
            mThroat.value( BELFEM_ENGINE_STATE_A )
                = 0.25 * std::pow( mParams.throat_diameter(), 2 )
                                                    * constant::pi ;

            this->create_gasmodels() ;
        }

//------------------------------------------------------------------------------

        Analysis::~Analysis()
        {
            if( mFuel != nullptr )
            {
                delete mFuel ;
            }
            if( mOxidizer != nullptr )
            {
                delete mOxidizer ;
            }
            if( mCombgas != nullptr )
            {
                delete mCombgas ;
            }
        }

//------------------------------------------------------------------------------

        real
        Analysis::compute_gas_generator( const real & aT,
                                         const real & aP,
                                         const real aOFmin,
                                         const real aOFmax,
                                         const real aT0 )
        {
            BELFEM_ERROR( aT > aT0, "Combustion temperature is too low" );

            real tX0 = aOFmin ;
            real tF0 = this->compute_combustion_temperature( aP, tX0, aT0 ) - aT ;

            real tX1 = aOFmax ;
            real tF1 = this->compute_combustion_temperature( aP, tX1, aT0 ) - aT ;

            BELFEM_ERROR( tF0 * tF1 < 0, "Invalid OF range" );

            real aX ;
            real tF = 1.0 ;

            // initialize counter
            uint tCount = 0;

            // start regula falsi
            while( abs( tF ) > 1e-6 )
            {
                aX = tX0 - tF0 * ( tX1 - tX0 ) / ( tF1 - tF0 );
                tF = this->compute_combustion_temperature( aP, aX, aT0 ) - aT ;
                if ( tF * tF0  > 0 )
                {
                    tX0 = aX ;
                    tF0 = tF ;
                }
                else
                {
                    tX1 = aX ;
                    tF1 = tF ;
                }
                BELFEM_ERROR( tCount++ < 100, "Too many iterations");

            }
            return aX ;
        }

//------------------------------------------------------------------------------

        real
        Analysis::run( const real aOF, const IspMode aMode )
        {
            // create timer
            Timer tTimer ;

            // check if a mixture has been set, otherwise
            // take the one from the parameters
            real tOF = aOF == 0.0 ? mParams.mixture_ratio() : aOF ;


            this->compute_injector( tOF, mParams.chamber_pressure() );

            this->compute_total();

            this->compute_throat();

            if( mParams.nozzle_mode() == NozzleMode::ComputeExitPressure )
            {
                this->compute_nozzle_aconst();
            }
            else
            {
                this->compute_nozzle_pconst() ;
            }

            // take ambient pressure from params
            mPambient = mNozzle.p() ;

            this->compute_performance();


            switch( aMode )
            {
                case( IspMode::Sealevel ) :
                {
                    return  mISPsl ;
                }
                case( IspMode::OptimalExpansion ) :
                {
                    return mISPref ;
                }
                case( IspMode::Vacuum ) :
                {
                    return mISPvac ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid Mode");
                    return BELFEM_QUIET_NAN ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Analysis::create_gasmodels()
        {

            // create oxidizer
            switch ( mParams.oxidizer() )
            {
                case( Oxidizer::LOX ) :
                {
                    mOxidizer = new Gas( HelmholtzModel::Oxygen );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown oxidizer: %s",
                                 oxidizer_to_string( mParams.oxidizer() ).c_str() );
                }
            }

            // create fuel
            switch ( mParams.fuel() )
            {
                case( Fuel::LH2 ) :
                {
                    mFuel = new Gas( HelmholtzModel::ParaHydrogen );
                    break ;
                }
                case( Fuel::LCH4 ) :
                {
                    mFuel = new Gas( HelmholtzModel::Methane );
                    break ;
                }
                default:
                {
                    // temporary molar fractions
                    mFuel = new Gas( mParams.fuel_species(),
                                     mParams.fuel_molar_fractions(),
                                     GasModel::SRK );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Analysis::compute_injector( const real & aOF, const real & aP )
        {
            // link to states
            real & tT = mInjector.value( BELFEM_ENGINE_STATE_T );
            real & tP = mInjector.value( BELFEM_ENGINE_STATE_P );
            real & tH = mInjector.value( BELFEM_ENGINE_STATE_H );
            real & tS = mInjector.value( BELFEM_ENGINE_STATE_S );

            // initialize the mixture
            this->create_initial_mixture( aOF, aP );

            // compute the temperatue of the mixture
            tT = 200.0 ;
            tP = aP ;

            real tDeltaT = BELFEM_REAL_MAX ;

            uint tCount = 0 ;

            // temperature for this enthalpy
            while( std::abs( tDeltaT ) > 1e-7 )
            {
                tDeltaT = ( mCombgas->h( tT, tP ) - tH ) / mCombgas->cp( tT, tP );

                tT -= 0.9 * tDeltaT ;

                BELFEM_ERROR( tCount++ < 1000 ,
                             "Failed to compute mixture temperature" );

            }

            // entropy of mixture
            tS = mCombgas->s( tT, tP );
        }

//------------------------------------------------------------------------------

        void
        Analysis::compute_total( const real aInitialTemperatureGuess )
        {
            mTotal.compute_equilibrium( aInitialTemperatureGuess,
                                        mInjector.p(), mInjector.h() );

            mTotal.compute_caloric( mTotal.T(),mTotal.p(), 0.0 );
        }

//------------------------------------------------------------------------------

        void
        Analysis::compute_throat()
        {
            // link to states
            real & tT     = mThroat.value( BELFEM_ENGINE_STATE_T );
            real & tP     = mThroat.value( BELFEM_ENGINE_STATE_P );

            // speed of sound
            real tW ;
            real  tH0 = mTotal.h() ;
            real  tS0 = mTotal.s();

            real tOmega = 0.5 ;

            // initial guess for throat
            tT = mTotal.T() / ( 1.0 + 0.5 * ( mTotal.gamma() - 1.0 ) );

            // initial guess for pressure
            tP = mTotal.p() * std::pow( tT / mTotal.T(), mTotal.gamma()/( mTotal.gamma() - 1.0 )  ) ;


            real tH ;
            real tRes = BELFEM_REAL_MAX ;

            uint tCount = 0;

            while ( tRes > 1e-9 )
            {
                mCombgas->remix( mTotal.molar_fractions(), false, false );
                mCombgas->remix_to_equilibrium( tT, tP, true, false );

                tH = mCombgas->h( tT, tP );

                // speed of sound
                tW = mCombgas->c( tT, tP );

                // Jacobian :

                // dhdT
                mJ2( 0, 0 ) = mCombgas->cp( tT, tP ) + 0.5 * mCombgas->R( tT, tP ) * mCombgas->gamma( tT, tP );

                // dsdt
                mJ2( 1, 0 ) = mCombgas->dsdT( tT, tP );

                // dhdp
                mJ2( 0, 1 ) = 0.0;

                mJ2( 1, 1 ) = mCombgas->dsdp( tT, tP );

                // rhs
                mF2( 0 ) = tH + 0.5 * tW * tW - tH0;
                mF2( 1 ) = mCombgas->s( tT, tP ) - tS0;


                tRes = std::sqrt( std::pow( mF2( 0 ) / tH0, 2 )
                                  + std::pow( mF2( 1 ) / tS0, 2 ));

                gesv( mJ2, mF2, mPivot );

                tT -= tOmega * mF2( 0 );
                tP -= tOmega * mF2( 1 );

                BELFEM_ERROR( tCount++ < 200, "Failed to compute throat state" );
            }

            mCombgas->remix( mCombgas->molar_fractions(),
                             true,
                             true );

            tRes = BELFEM_REAL_MAX ;


            mThroat.compute_caloric( tT, tP, mCombgas->c( tT, tP ) );

            mThroat.mass_fractions()  = mCombgas->mass_fractions() ;
            mThroat.molar_fractions() = mCombgas->molar_fractions() ;
        }

//------------------------------------------------------------------------------

        void
        Analysis::create_initial_mixture( const real & aOF, const real & aP  )
        {
            // remember mixture
            mOF = aOF ;

            Vector< real > & tY = mInjector.mass_fractions() ;

            // - - - - - - - - - - - - - - - - - - - - - - - - -
            // compute initial enthalpies for the mixture
            // - - - - - - - - - - - - - - - - - - - - - - - - -

            // set the oxidizer in order to compute the oxidizer enthalpy
            const Vector< real > & tYO = mOxidizer->mass_fractions() ;
            tY.fill( 0.0 );
            index_t tCount = 0 ;
            for( index_t k : mParams.oxidizer_indices() )
            {
                tY( k ) = tYO( tCount++ );
            }
            mCombgas->remix_mass( tY, true, false );
            real tOxidizerEnthalpy = mCombgas->h( mParams.oxidizer_temperature(), aP );

            // set the fuel
            const Vector< real > & tYF = mFuel->mass_fractions() ;
            tCount = 0 ;
            tY.fill( 0.0 );
            for( index_t k : mParams.fuel_indices() )
            {
                tY( k ) = tYF( tCount++ );
            }
            mCombgas->remix_mass( tY, true, false );
            real tFuelEnthalpy = mCombgas->h( mParams.fuel_temperature(), aP );

            // add value for oxidizer
            tCount = 0.0 ;
            for( index_t k : mParams.oxidizer_indices() )
            {
                tY( k ) += tYO( tCount++ ) * aOF ;
            }

            // create the mixture
            mCombgas->remix_mass( tY, true, false );

            // save enthalpy and pressure into values list
            mInjector.value( BELFEM_ENGINE_STATE_P ) = aP ;
            mInjector.value( BELFEM_ENGINE_STATE_H ) =
                    ( tFuelEnthalpy + aOF * tOxidizerEnthalpy ) / ( 1.0 + aOF );
        }

//------------------------------------------------------------------------------

        void
        Analysis::compute_nozzle_pconst()
        {
            // relaxation factor
            real tOmega = 0.5 ;

            // mass flux
            real tDotM0 = mThroat.rho() * mThroat.u() ;

            // entropy and total enthalpy
            //real tS0 = mTotal.s() ;
            //real tH0 = mTotal.h(); //mThroat.h() + 0.5 * mThroat.u() * mThroat.u() ;

            real tS0 = mThroat.s() ;
            real tH0 = mThroat.h() + 0.5 * mThroat.u() * mThroat.u() ;

            real & tT     = mNozzle.value( BELFEM_ENGINE_STATE_T     );
            real & tP     = mNozzle.value( BELFEM_ENGINE_STATE_P     );

            // get pressure from parameters
            tP = mParams.exit_pressure() ;

            // guess value for p
            tT = mThroat.T() * std::pow( tP / mThroat.p(),
                                         ( mThroat.gamma() - 1.0 ) / mThroat.gamma() );

            // guess value for Mach number
            real tMa = std::sqrt ( 2.0 / ( mThroat.gamma() - 1.0 ) * ( mTotal.T() / tT - 1. ) );

            // initial guess for velocity
            real tU = tMa * mCombgas->c( tT, tP );

            // density
            real tRho = mCombgas->rho( tT, tP );

            // initial guess for cross section
            real tA = tDotM0 / ( tU * tRho );

            real tRes = BELFEM_REAL_MAX ;
            uint tCount = 0 ;

            real tH = mCombgas->h( tT, tP );

            while( tRes > 1e-6 )
            {
                mCombgas->remix( mThroat.molar_fractions(), false, false );
                mCombgas->remix_to_equilibrium( tT, tP, true, false );

                // compute the Jacobian
                mJ3( 0, 0 ) = mCombgas->dsdT( tT, tP );
                mJ3( 1, 0 ) = mCombgas->cp( tT, tP );
                mJ3( 2, 0 ) = -tRho * mCombgas->alpha( tT, tP ) * tU * tA ;
                mJ3( 0, 1 ) = 0.0;
                mJ3( 1, 1 ) = tU;
                mJ3( 2, 1 ) = tRho * tA;
                mJ3( 0, 2 ) = 0.0;
                mJ3( 1, 2 ) = 0.0;
                mJ3( 2, 2 ) = tRho * tU;

                // compute the right hand side
                mF3( 0 ) = mCombgas->s( tT, tP ) - tS0;
                mF3( 1 ) = tH + 0.5 * tU * tU - tH0;
                mF3( 2 ) = tRho * tA * tU - tDotM0;

                tRes = std::sqrt(   std::pow( mF3( 0 ) / tS0 , 2 )
                                  + std::pow( mF3( 1 ) / tH0 , 2 )
                                  + std::pow( mF3( 2 )/ tDotM0 ,2 ) );

                gesv( mJ3, mF3, mPivot );

                tT -= tOmega * mF3( 0 );
                tU -= tOmega * mF3( 1 );
                tA -= tOmega * mF3( 2 );

                tH = mCombgas->h( tT, tP );

                tRho = mCombgas->rho( tT, tP );

                // adative relaxation
                tOmega = 1.0 - 0.11 * std::log( ++tCount );

                BELFEM_ERROR( tCount < 500 , "Failed to converge in Nozzle computation" );

            }

            mCombgas->remix_to_equilibrium( tT, tP, true, false );

            mNozzle.compute_caloric( tT, tP , tU );

            mNozzle.value( BELFEM_ENGINE_STATE_A ) = tA * mThroat.A() ;

        }

//------------------------------------------------------------------------------

        void
        Analysis::compute_nozzle_aconst()
        {
            real & tT     = mNozzle.value( BELFEM_ENGINE_STATE_T     );
            real & tP     = mNozzle.value( BELFEM_ENGINE_STATE_P     );
            real & tU     = mNozzle.value( BELFEM_ENGINE_STATE_U     );
            real & tH     = mNozzle.value( BELFEM_ENGINE_STATE_H     );
            real & tRho   = mNozzle.value( BELFEM_ENGINE_STATE_RHO   );
            real & tA     = mNozzle.value( BELFEM_ENGINE_STATE_A     );

            // reference conditions
            real tS0 = mThroat.s() ;
            real tH0 = mThroat.h() + 0.5 * mThroat.u() * mThroat.u() ;

            // compute initial guess for exit conditions
            // see CEA Eq. 6.21 and 6.22
            real tVal ;
            if( mParams.expansion_ratio() < 2.0 )
            {
                tVal = std::log( mParams.expansion_ratio() );

                tVal = std::log( mTotal.p() / mThroat.p() )
                        + std::sqrt( tVal * ( 1.535 +
                                     3.294 * tVal ) );
            }
            else
            {
                tVal = mThroat.gamma()
                        + 1.4 * std::log( mParams.expansion_ratio() );

            }

            tP = mTotal.p() / std::exp( tVal );

            tT = mCombgas->isen_T( mThroat.T(), mThroat.p(), tP );

            tH = mCombgas->h( tT, tP );
            tU = std::sqrt( 2.0 * ( tH0 - tH ) );

            tA = mThroat.A() * mParams.expansion_ratio() ;

            tRho = mCombgas->rho( tT, tP );

            real tDotM0 = mThroat.rho() * mThroat.A() * mThroat.u() ;

            real tRes = BELFEM_REAL_MAX ;
            real tOmega = 0.5 ;

            uint tCount = 0 ;

            while( tRes > 1e-6 )
            {
                mCombgas->remix( mThroat.molar_fractions(), false, false );
                mCombgas->remix_to_equilibrium( tT, tP, true, false );

                // compute the right hand side
                mF3( 0 ) = mCombgas->s( tT, tP ) - tS0;
                mF3( 1 ) = tH + 0.5 * tU * tU - tH0;
                mF3( 2 ) = tRho * tA * tU - tDotM0;

                // dsdT
                mJ3( 0, 0 ) = mCombgas->dsdT( tT, tP );

                // dhdT
                mJ3( 1, 0 ) = mCombgas->cp( tT, tP );

                // dmdT
                mJ3( 2, 0 ) = - tRho * tA * tU * mCombgas->alpha( tT, tP );

                // dsdp
                mJ3( 0, 1 ) = mCombgas->dsdp( tT, tP );

                // dhdp
                mJ3( 1, 1 ) = 0.0 ;

                // dmdp
                mJ3( 2, 1 ) = tRho * tA * tU * mCombgas->kappa( tT, tP );

                // dsdu
                mJ3( 0, 2 ) = 0.0 ;

                // dhdu
                mJ3( 1, 2 ) = tU ;

                // dmdu
                mJ3( 2, 2 ) = tRho * tA ;

                tRes = std::sqrt(     std::pow( mF3( 0 ) / tS0 , 2 )
                                    + std::pow( mF3( 1 ) / tH0 , 2 )
                                    + std::pow( mF3( 2 )/ tDotM0 ,2 ) );

                gesv( mJ3, mF3, mPivot );

                tT -= tOmega * mF3( 0 );
                tP -= tOmega * mF3( 1 );
                tU -= tOmega * mF3( 2 );

                tRho = mCombgas->rho( tT, tP );
                tH   = mCombgas->h( tT, tP );

                tOmega = 1.0 - 0.11 * std::log( ++tCount );

                BELFEM_ERROR( tCount < 100 , "Failed to converge in Nozzle computation" );

            }

            mCombgas->remix_to_equilibrium( tT, tP, true, false );

            mNozzle.compute_caloric( tT, tP , tU );
        }

//------------------------------------------------------------------------------

        void
        Analysis::compute_performance()
        {
            // mass flow
            mDotM = mNozzle.u() * mNozzle.rho() * mNozzle.A() ;

            // thrust
            mF = mDotM * mNozzle.u() + (
                    mNozzle.p() - mPambient ) * mNozzle.A() ;

            // specific impulse
            mISPref = mF / ( mDotM * constant::g0 );

            mISPsl  = ( mDotM * mNozzle.u() + (
                      mNozzle.p() - 1.01325e5 ) * mNozzle.A() ) /
                       ( mDotM * constant::g0 );

            mISPvac = ( mDotM * mNozzle.u() + mNozzle.p() * mNozzle.A() ) /
                       ( mDotM * constant::g0 );

            // thrust coefficient
            mCF = mF / ( mTotal.p() * mThroat.A() );
        }

//------------------------------------------------------------------------------
        real
        Analysis::find_best_mixture(
                const real & aOFmin,
                const real & aOFmax,
                const IspMode aMode )
        {
            // do a stepwise search
            uint tNumSteps = 5 ;

            Vector< real > tOF( tNumSteps );
            Vector< real > tISP( tNumSteps );

            linspace( aOFmin, aOFmax, tNumSteps, tOF );

            real tISPmax = 0.0;
            uint tK = BELFEM_UINT_MAX;

            for ( uint k = 0; k < tNumSteps; ++k )
            {
                // perform analysis
                tISP( k ) = this->run( tOF( k ), aMode );

                // remember best value
                if ( tISP( k ) > tISPmax )
                {
                    tK = k;
                    tISPmax = tISP( k );
                }
            }

            // initial values

            real tA = tK > 1 ? tOF( tK-1 ) : tOF( 0 );
            real tB = tK < tNumSteps-3 ? tOF( tK+1 ) : tOF( tNumSteps-1 );

            // function values
            real tFa = ( this->run( 1.05 * tA ) - this->run( 0.95 * tA ) ) / ( 0.1 * tA );
            //real tFb = ( this->run( 1.05 * tB ) - this->run( 0.95 * tB ) ) / ( 0.1 * tB );

            real tF = BELFEM_REAL_MAX ;
            real aOF ;
            uint tCount = 4 + tNumSteps ;
            while( std::abs( tF ) > 1e-7 )
            {
                aOF = 0.5 * ( tA + tB );

                tF =  ( this->run( 1.05 * aOF, aMode ) - this->run( 0.95 * aOF, aMode ) ) / ( 0.1 * aOF );
                if ( tF * tFa > 0 )
                {
                    tA = aOF;
                    tFa = tF;
                }
                else
                {
                    tB = aOF;
                    //tFb = tF;
                }
                tCount += 2 ;
            }

            return aOF ;
        }

//------------------------------------------------------------------------------

        // help function for gas generator mode
        real
        Analysis::compute_combustion_temperature( const real & aP,
                                                  const real & aOF,
                                                  const real   aT0 )
        {
            this->compute_injector( aOF, aP );
            this->compute_total( aT0 ) ;
            return this->total()->T();
        }

//------------------------------------------------------------------------------

        void
        Analysis::print_performance()
        {
            std::fprintf( stdout, "     Ideal Performance : \n\n" );

            std::fprintf( stdout, "         Mixture Ratio  : %10.3f -\n\n",
                          mOF );

            std::fprintf( stdout, "         Mass Flow      : %10.3f kg/s\n\n",
                          mDotM );

            std::fprintf( stdout, "         Thrust @ p_amb : %10.3f kN \n\n",
                          mF * 0.001 );

            std::fprintf( stdout, "         CF  @ p_amb    : %10.3f - \n\n",
                          mCF );

            std::fprintf( stdout, "         ISP @ p_amb    : %10.3f s \n\n",
                          mISPref );

            std::fprintf( stdout, "         ISP @ SL       : %10.3f s \n\n",
                          mISPsl );

            std::fprintf( stdout, "         ISP @ VAC      : %10.3f s \n\n",
                          mISPvac );
        }

//------------------------------------------------------------------------------
    }
}