//
// Created by Christian Messe on 01.04.21.
//

#include "cl_EN_Turbine.hpp"
#include "assert.hpp"
#include "constants.hpp"
#include "fn_polyval.hpp"

namespace belfem
{
    namespace engine
    {
//------------------------------------------------------------------------------

        Turbine::Turbine( Gas & aGas ) :
            mGas( aGas ),
            mNozzleEntry( aGas, "Nozzle Entry"),
            mTurbineEntry( aGas, "Turbine Entry"),
            mTurbineDischarge( aGas, "Turbine Discharge" ),
            mTurbineEntryRotating( aGas, "Turbine Entry Rotating" ),
            mTurbineDischargeRotating( aGas, "Turbine Discharge Rotating" )
    {
            // polynomial for eta_la, David Fig. 3
            mEtaPoly = {  -2.216822e-3, 1.000617e-2, -2.272745e-2,
                          -5.187984e-3,0.0, 0.930243 } ;

            // polynomial for pitch-chord ratio, Bohl 2 Fig. 2.45, p. 99
            mPitchChordPoly = { -0.06537458, 0.4604635, - 0.9772311, 0.1874318, 1.628015, 0.0 };

            // temperature range
            mHmin = mGas.h( 200.0, 1e5 );
            mHmax = BELFEM_QUIET_NAN;

            this->reset();

        }

//------------------------------------------------------------------------------

        void
        Turbine::set_entry( const real & aTt, const real & aPt )
        {
            mEntryFlag = true ;
            mNozzleEntry.value( BELFEM_ENGINE_STATE_TT ) = aTt ;
            mNozzleEntry.value( BELFEM_ENGINE_STATE_PT ) = aPt ;
            mNozzleEntry.value( BELFEM_ENGINE_STATE_HT ) = mGas.h( aTt, aPt );
            mNozzleEntry.value( BELFEM_ENGINE_STATE_S )  = mGas.s( aTt, aPt );

            // can't get hotter than turbine entry
            mHmax = mGas.h( aTt*1.1, aPt );
        }

//------------------------------------------------------------------------------

        bool
        Turbine::compute()
        {
            // reset error
            mErrorCode = 0;

            this->check_sanity();
            this->compute_diameter();

            real tX0 = 0.2;
            real tF0;

            real tX1 ;
            real tF1 ;
            real tDeltaX = 1.0 ;

            for ( uint k=0 ; k<3; ++k )
            {
                tDeltaX *= 0.1 ;
                tX0 -= tDeltaX ;
                if ( this->compute_blade_height( tX0, tF0 ) )
                {
                    this->reset( 1 );
                    return false;
                }

                tX1 = tX0 ;
                tF1 = tF0 ;

                while ( tF0 * tF1  > 0 )
                {
                    tX0 = tX1 ;
                    tF0 = tF1 ;
                    tX1 += tDeltaX ;
                    if ( tX1 > 2.0 )
                    {
                        this->reset( 2 );
                        return false;
                    }
                    if ( this->compute_blade_height( tX1, tF1 ) )
                    {
                        this->reset( 3 );
                        return false;
                    }
                }
            }

            // "Invalid condition"
            if( tF0 * tF1 > 0 )
            {
                this->reset( 4 ) ;
                return false ;
            }

            real tF = 1.0 ;
            uint tCount = 0 ;
            real tX ;

            while( abs( tF ) > 1e-7 )
            {
                tX = tX0 - tF0 * ( tX1 - tX0 ) / ( tF1 - tF0 ) ;
                //tX = 0.5 * ( tX0 + tX1 );

                if( this->compute_blade_height( tX, tF ) )
                {
                    this->reset( 5 ) ;
                    return false ;
                }
                if ( tF * tF0 > 0 )
                {
                    tX0 = tX ;
                    tF0 = tF ;
                }
                else
                {
                    tX1 = tX ;
                    tF1 = tF ;
                }

                if( tCount++ > 100 )
                {
                    this->reset( 6 ) ;
                    return false ;
                }


            }

            // finalize
            if( this->compute_nozzle_entry() )
            {
                this->reset( 7 ) ;
                return false ;
            }
            if( this->compute_rotating_states() )
            {
                this->reset( 8 ) ;
                return false ;
            }

            this->compute_blade_entry();
            return true ;
        }

//------------------------------------------------------------------------------

        bool
        Turbine::compute_blade_height( const real & aMu, real & aDeltaB  )
        {

            if ( mPhi1Flag )
            {
                if( this->compute_turbine_entry( aMu, mPhi1, mEpsilon ) )
                {
                    return true;
                }
            }
            else if ( mEpsilonFlag )
            {

                real tX1 = 0.4 ;
                real tF1 ;

                if( this->compute_turbine_entry( aMu, tX1, tF1 ) )
                {
                    return true ;
                }
                tF1 -= mEpsilon ;

                real tX0 ;
                real tF0 = tF1 ;

                while (  tF0 * tF1 > 0 )
                {
                    tX0 = tX1 ;
                    tF0 = tF1 ;
                    tX1 += 0.1 ;

                    if( this->compute_turbine_entry( aMu, tX1, tF1 ) )
                    {
                        return false ;
                    }
                    tF1 -= mEpsilon ;

                    BELFEM_ERROR( tX1 < 2.0 , "Invalid value of admission factor" );
                }
                real tX ;
                real tF = 1 ;

                index_t tCount = 0;

                while ( abs( tF ) > 1e-6 )
                {
                    //tX = 0.5 * ( tX1 + tX0 );
                    tX = tX0 - tF0 * ( tX1 - tX0 ) / ( tF1 - tF0 );

                    if( this->compute_turbine_entry( aMu, tX, tF ) || tCount++ > 100 )
                    {
                        return false ;
                    }
                    tF -= mEpsilon ;

                    if ( tF * tF0  > 0.0 )
                    {
                        tX0 = tX ;
                        tF0 = tF ;
                    }
                    else
                    {
                        tX1 = tX ;
                        tF1 = tF ;
                    }
                }

            }
            else
            {
                BELFEM_ERROR( false, "Invalid or incomplete parameters.");
            }

            this->compute_rotor_efficiency();

            // compute and abort if error has been found
            if( this->predict_turbine_discharge() )
            {
                return true ;
            }


            this->compute_pitch_and_chord();

            real tDeltaEta = 1.0 ;

            uint tCount = 0 ;

            while( tDeltaEta > 1e-6 )
            {
                if( this->correct_turbine_discharge( tDeltaEta, 0.2 ) || tCount++ > 100 )
                {
                    return true ;
                }
            }

            aDeltaB = this->finalize_turbine_discharge();

            return false ;
        }

//------------------------------------------------------------------------------

        bool
        Turbine::compute_rotating_states()
        {
            real & Tt1  = mTurbineEntryRotating.value( BELFEM_ENGINE_STATE_TT );
            real & pt1  = mTurbineEntryRotating.value( BELFEM_ENGINE_STATE_PT );
            real & ht1  = mTurbineEntryRotating.value( BELFEM_ENGINE_STATE_HT );

            real & T1  = mTurbineEntryRotating.value( BELFEM_ENGINE_STATE_T );
            real & p1  = mTurbineEntryRotating.value( BELFEM_ENGINE_STATE_P );
            real & h1  = mTurbineEntryRotating.value( BELFEM_ENGINE_STATE_H );
            real & s1  = mTurbineEntryRotating.value( BELFEM_ENGINE_STATE_S );

            real & w1  = mTurbineEntryRotating.value( BELFEM_ENGINE_STATE_U );
            real & Ma1  = mTurbineEntryRotating.value( BELFEM_ENGINE_STATE_MA );

            real & Tt2  = mTurbineDischargeRotating.value( BELFEM_ENGINE_STATE_TT );
            real & pt2  = mTurbineDischargeRotating.value( BELFEM_ENGINE_STATE_PT );
            real & ht2  = mTurbineDischargeRotating.value( BELFEM_ENGINE_STATE_HT );

            real & T2  = mTurbineDischargeRotating.value( BELFEM_ENGINE_STATE_T );
            real & p2  = mTurbineDischargeRotating.value( BELFEM_ENGINE_STATE_P );
            real & h2  = mTurbineDischargeRotating.value( BELFEM_ENGINE_STATE_H );
            real & s2  = mTurbineDischargeRotating.value( BELFEM_ENGINE_STATE_S );

            real & w2  = mTurbineDischargeRotating.value( BELFEM_ENGINE_STATE_U );
            real & Ma2  = mTurbineDischargeRotating.value( BELFEM_ENGINE_STATE_MA );


            // copy state parameters for state 1
            T1  = mTurbineEntry.T();
            p1  = mTurbineEntry.p();
            h1  = mTurbineEntry.h();
            s1  = mTurbineEntry.s();
            w1  = mW1;

            // compute Mach number and total state
            Ma1 = w1 / mGas.c( T1, p1 );

            // catch error
            if ( std::isnan( Ma1 ) )
            {
                return true ;
            }

            mGas.total( T1, p1, w1, Tt1, pt1 );
            ht1 = h1 + 0.5 * w1 * w1 ;
            // copy state parameters for state 1
            T2  = mTurbineDischarge.T();
            p2  = mTurbineDischarge.p();
            h2  = mTurbineDischarge.h();
            s2  = mTurbineDischarge.s();
            w2  = mW2;

            // compute Mach number and total state
            Ma2 = w2 / mGas.c( T2, p2 );

            // catch error
            if ( std::isnan( Ma2 ) )
            {
                return true ;
            }

            mGas.total( T2, p2, w2, Tt2, pt2 );
            ht2 = h2 + 0.5 * w2 * w2 ;
            return false ;
        }

//------------------------------------------------------------------------------

        void
        Turbine::compute_blade_entry()
        {

            mBladeEntry = mPitch * std::sin( mBeta1 );

            mBladeEntryError = mBladeOpening * mTurbineDischarge.rho() * mW2 * mB2B1 /
                    ( mTurbineEntry.rho() * mW1 ) - mBladeEntry ;

        }

//------------------------------------------------------------------------------

        void
        Turbine::reset( const int aErrorCode )
        {
            mErrorCode = aErrorCode ;

            mYs  = BELFEM_QUIET_NAN ;
            mPsi = BELFEM_QUIET_NAN ;
            mBD  = BELFEM_QUIET_NAN ;
            mB0   = BELFEM_QUIET_NAN ;
            mB1   = BELFEM_QUIET_NAN ;
            mB2   = BELFEM_QUIET_NAN ;
            mUm  = BELFEM_QUIET_NAN ;
            mDm  = BELFEM_QUIET_NAN ;

            mPhi0 = BELFEM_QUIET_NAN ;
            mPhi1 = BELFEM_QUIET_NAN ;
            mPhi2 = BELFEM_QUIET_NAN ;

            mA0 = BELFEM_QUIET_NAN ;
            mA1 = BELFEM_QUIET_NAN ;
            mA2 = BELFEM_QUIET_NAN ;

            mCm0 = BELFEM_QUIET_NAN ;
            mCm1 = BELFEM_QUIET_NAN ;
            mCm2 = BELFEM_QUIET_NAN ;
            mC1  = BELFEM_QUIET_NAN ;
            mC2  = BELFEM_QUIET_NAN ;
            mCu2 = BELFEM_QUIET_NAN ;
            mCu1 = BELFEM_QUIET_NAN ;
            mWu1 = BELFEM_QUIET_NAN ;
            mWu2 = BELFEM_QUIET_NAN ;
            mW1 = BELFEM_QUIET_NAN ;
            mW2 = BELFEM_QUIET_NAN ;

            mAlpha1 = BELFEM_QUIET_NAN ;
            mBeta1  = BELFEM_QUIET_NAN ;
            mBeta2  = BELFEM_QUIET_NAN ;

            mDeltaHsNozzle = BELFEM_QUIET_NAN ;
            mDeltaHsRotor = BELFEM_QUIET_NAN ;
            mReaction = BELFEM_QUIET_NAN ;

            mEtaRotor = BELFEM_QUIET_NAN ;
            mEtaRotorFullAdmission = BELFEM_QUIET_NAN ;

            mEpsilon = BELFEM_QUIET_NAN ;

            mPsiFlag = false ;
            mBDflag  = false ;
            mBflag   = false ;
            mPhi1Flag = false ;

            if( ! mZ2Flag )
            {
                mZ2 = BELFEM_UINT_MAX ;
            }
            if( ! mPitchChordRatioFlag )
            {
                mPitchChordRatio = BELFEM_QUIET_NAN ;
            }

            mPitch = BELFEM_QUIET_NAN ;
            mChord = BELFEM_QUIET_NAN ;
            mAxialChord = BELFEM_QUIET_NAN ;
            mBladeOpening = BELFEM_QUIET_NAN ;
            mBladeRadius1 = BELFEM_QUIET_NAN ;
            mBladeRadius2 = BELFEM_QUIET_NAN ;
            mChordAngle = BELFEM_QUIET_NAN ;


            mBladeEntry = BELFEM_QUIET_NAN ;
            mBladeEntryError = BELFEM_REAL_MAX ;

            mDeltaEtaVentilation = BELFEM_QUIET_NAN ;
            mDeltaEtaMixingAndExpansion = BELFEM_QUIET_NAN ;
            mEtaFullAdmission  = BELFEM_QUIET_NAN ;
            mEta  = BELFEM_QUIET_NAN ;

            mTt2s = BELFEM_QUIET_NAN ;
            mHt2s = BELFEM_QUIET_NAN ;

            real tT0 = mNozzleEntry.Tt() ;
            real tP0 = mNozzleEntry.pt() ;
            real tH0 = mNozzleEntry.ht() ;
            real tS  = mNozzleEntry.s() ;

            for( uint k=0; k<BELFEM_ENGINE_NUMSTATES; ++k )
            {
                mNozzleEntry.value( k ) = BELFEM_QUIET_NAN ;
                mTurbineEntryRotating.value( k ) = BELFEM_QUIET_NAN ;
                mTurbineDischarge.value( k ) = BELFEM_QUIET_NAN;
                mTurbineDischargeRotating.value( k ) = BELFEM_QUIET_NAN;
            }
            mNozzleEntry.value( BELFEM_ENGINE_STATE_TT ) = tT0 ;
            mNozzleEntry.value( BELFEM_ENGINE_STATE_PT ) = tP0 ;
            mNozzleEntry.value( BELFEM_ENGINE_STATE_HT ) = tH0 ;
            mNozzleEntry.value( BELFEM_ENGINE_STATE_S )  = tS ;
        }

//------------------------------------------------------------------------------

        void
        Turbine::check_sanity()
        {
            BELFEM_ERROR( mEntryFlag, "entry conditions have not been set" );

            BELFEM_ERROR( mNflag, "shaft speed n was not set" );
            BELFEM_ERROR( mDotMflag, "massflow was not set");

            BELFEM_ERROR( mYflag || mPflag,  "either power P or specific power Y must be set");

            BELFEM_ERROR( ! ( mYflag && mPflag ),
                         "power P and specific power Y must not be prescribet at the same time" );

            BELFEM_ERROR( mPsiFlag, "parameter psi was not set");

            BELFEM_ERROR( mPhi1Flag || mEpsilonFlag ,
                         "either phi1 or admission factor epsilon must be set" );

            BELFEM_ERROR( ! ( mPhi1Flag && mEpsilonFlag ),
                         "phi1 and blade admission factor epsilo must not be prescribed at the same time" );


            BELFEM_ERROR( ! ( mBflag && mBDflag ),
                         "Blade height b and blade diameter ratio bD must not be prescribed at the same time" );
            BELFEM_ERROR( mBflag || mBDflag ,
                         "either blade height b or blade diameter ratio must be set");
        }

//------------------------------------------------------------------------------

        void
        Turbine::compute_diameter()
        {
            if( mPflag )
            {
                mY = mP / mDotM ;
            }
            else if ( mYflag )
            {
                mP = mY * mDotM ;
            }
            else
            {
                BELFEM_ERROR( false, "Invalid call");
            }

            // radial velocity
            mUm = std::sqrt( mY / mPsi );

            // diameter
            mDm = mUm * 60.0 / ( constant::pi * mN );

            if ( mBflag )
            {
                mBD = mB1 / mDm ;
            }
            else if ( mBDflag )
            {
                mB1 = mBD * mDm ;
            }
            else
            {
                BELFEM_ERROR( false, "Invalid call");
            }
        }

//------------------------------------------------------------------------------

        bool
        Turbine::compute_turbine_entry( const real & aMu, const real & aPhi1, real & aEpsilon )
        {
            // remember phi
            mPhi1 = aPhi1 ;

            // get shortcuts
            const real & Tt0 = mNozzleEntry.value( BELFEM_ENGINE_STATE_TT );
            const real & pt0 = mNozzleEntry.value( BELFEM_ENGINE_STATE_PT );
            const real & ht0 = mNozzleEntry.value( BELFEM_ENGINE_STATE_HT );

            real & Tt1 = mTurbineEntry.value( BELFEM_ENGINE_STATE_TT );
            real & pt1 = mTurbineEntry.value( BELFEM_ENGINE_STATE_PT );
            real & ht1 = mTurbineEntry.value( BELFEM_ENGINE_STATE_HT );

            real & T1  = mTurbineEntry.value( BELFEM_ENGINE_STATE_T );
            real & p1  = mTurbineEntry.value( BELFEM_ENGINE_STATE_P );
            real & h1  = mTurbineEntry.value( BELFEM_ENGINE_STATE_H );
            real & s1  = mTurbineEntry.value( BELFEM_ENGINE_STATE_S );
            real & rho1  = mTurbineEntry.value( BELFEM_ENGINE_STATE_RHO );
            real & Ma1   = mTurbineEntry.value( BELFEM_ENGINE_STATE_MA );

            // assume no heat loss in nozzle
            Tt1 = Tt0 ;
            ht1 = ht0 ;

            mPhi2 = aMu * aPhi1 ;

            // axial velocity
            mCm1 = aPhi1 * mUm ;
            mCm2 = mPhi2 * mUm ;

            // trigonometry
            mC2 = mCm2 / std::sin( mAlpha2 );

            // David ( 28 )
            if ( std::abs( mAlpha2 - 0.5 * constant::pi ) > 1e-4 )
            {
                mCu2 = mUm * mPhi2 / std::tan( mAlpha2 );
            }
            else
            {
                mCu2 = 0.0 ;
            }

            // Euler / David ( 27 )
            mCu1 = mUm * mPsi + mCu2 ;

            // David ( 33 )
            mAlpha1 = std::atan( aPhi1 * mUm / mCu1 );

            // triginometry
            mC1 = mCm1 / std::sin( mAlpha1 );

            // David 1
            mDeltaHsNozzle = 0.5 * mC1 * mC1 / mEtaNozzle ;

            // compute help state at nozzle discharge
            real h1s = ht0 - mDeltaHsNozzle ;

            // check for validity
            if( this->check_temperature_range( h1s ) )
            {
                aEpsilon = BELFEM_QUIET_NAN ;
                return  true ;
            }

            // assuming ideal gas law
            real T1s = mGas.T_from_h( h1s, pt0 );

            // compute static pressure at nozzle discharge
            p1 = mGas.isen_p( Tt0, pt0, T1s );

            // static enthalpy at nozzle discharge / turbine entry
            h1  = ht1 - 0.5 * mC1 * mC1 ;


            // check for validity
            if( this->check_temperature_range( h1 ) )
            {
                aEpsilon = BELFEM_QUIET_NAN ;
                return  true ;

            }

            // static temperature at nozzle discharge / turbine entry
            T1 = mGas.T_from_h( h1, p1 );

            // check for plausibility
            if( T1 < T1s )
            {
                aEpsilon = BELFEM_QUIET_NAN ;
                return  true ;
            }

            // density
            rho1 = mGas.rho( T1, p1 );

            // total pressure at entry
            pt1 = mGas.isen_p( T1, p1, Tt1 );

            // entropy at entry
            s1 = mGas.s( Tt1, pt1 );

            // compute admission factor
            aEpsilon = mDotM /
                       ( mB1 * constant::pi * mDm * rho1 * aPhi1 * mUm );

            // compute cross section
            real tDo = mDm + mB1 ;
            real tDi = mDm - mB1 ;

            mA1 = aEpsilon * 0.25 * constant::pi * ( tDo * tDo - tDi * tDi );

            // Mach number at turbine entry
            Ma1 = mC1 / mGas.c( T1, p1 );

            return false ;
        }

//------------------------------------------------------------------------------

        void
        Turbine::compute_rotor_efficiency()
        {
            // David ( 3 )
            mWu1 = mCu1 - mUm ;
            mWu2 = mCu2 - mUm ;

            // Pythagoras
            mW1 = std::sqrt( mWu1 * mWu1 + mCm1 * mCm1 );
            mW2 = std::sqrt( mWu2 * mWu2 + mCm2 * mCm2 );

            // David ( 4 )
            mBeta1 = std::asin( mCm1 / mW1 );
            mBeta2 = std::acos( mWu2 / mW2 );


            mEtaRotorFullAdmission = polyval( mEtaPoly, std::abs( mBeta2 - mBeta1 ) );

            // initial guess
            mEtaRotor = mEtaRotorFullAdmission ;

            // initial guess: David ( 2 )
            mDeltaHsRotor = 0.5 * ( mW2 * mW2 / mEtaRotor - mW1 * mW1 );

            // initial guess
            mReaction = mDeltaHsRotor / ( mDeltaHsRotor + mDeltaHsNozzle );

        }

//------------------------------------------------------------------------------

        bool
        Turbine::predict_turbine_discharge()
        {
            // shortcuts
            real & Tt2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_TT );
            real & pt2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_PT );
            real & ht2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_HT );
            // real & s2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_S );

            real & T2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_T );
            real & p2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_P );
            real & h2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_H );
            real & rho2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_RHO );
            // real & Ma2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_MA );

            // reference state from enthalpy-entropy chart
            real h2s = mTurbineEntry.h() - mDeltaHsRotor ;

            // check temperature range
            if( this->check_temperature_range( h2s ) )
            {
                return true ;
            }

            real T2s = mGas.T_from_h( h2s, mTurbineEntry.p() );

            // static conditions at exit
            p2 = mGas.isen_p( mTurbineEntry.Tt(), mTurbineEntry.pt(), T2s );

            ht2 = mNozzleEntry.ht() - mY ;

            // check temperature range
            if( this->check_temperature_range( ht2 ) )
            {
                return true ;
            }

            h2  = ht2 - 0.5 * mC2 * mC2 ;

            // check temperature range
            if( this->check_temperature_range( h2 ) )
            {
                return true ;
            }

            T2  = mGas.T_from_h( h2, p2 );


            // total conditions at exit
            Tt2 = mGas.T_from_h( ht2, p2 );
            pt2 = mGas.isen_p( T2, p2, Tt2 );

            // reference state
            mTt2s = mGas.isen_T( mNozzleEntry.Tt(), mNozzleEntry.pt(), pt2 );
            mHt2s = mGas.h( mTt2s, pt2 );
            mYs = mNozzleEntry.ht() - mHt2s ;

            // efficiency for full admission
            mEtaFullAdmission = mY / mYs ;

            // density
            rho2 = mGas.rho( T2, p2 );

            // exit cross section
            mA2 = mDotM / ( mCm2 * rho2 );

            // exit blade height
            mB2 = mA2 / ( mEpsilon * constant::pi * mDm );

            return false ;
        }

//------------------------------------------------------------------------------

        bool
        Turbine::correct_turbine_discharge( real & aDeltEta, const real aOmega )
        {
            // shortcuts
            real & Tt2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_TT );
            real & pt2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_PT );

            real & T2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_T );
            real & p2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_P );

            // old rotor efficiency
            real tEtaRotor0 = mEtaRotor ;

            // compute additional losses
            this->compute_ventilation_losses();

            if( this->compute_mixing_and_expansion_losses() )
            {
                return false ;
            }

            // correct efficiency
            mEta = mEtaFullAdmission - mDeltaEtaMixingAndExpansion - mDeltaEtaVentilation ;

            // correct reference work
            mYs = mY / mEta ;

            // new isentropic reference state
            mHt2s = mNozzleEntry.ht() - mYs ;

            if( this->check_temperature_range( mHt2s ) )
            {
                return false ;
            }

            mTt2s = mGas.T_from_h( mHt2s, pt2 );

            // new total pressure (temperature Tt2 does not change)
            pt2 = mGas.isen_p( mNozzleEntry.Tt(), mNozzleEntry.pt(), mTt2s );

            // new static pressure (temperature T2 does not change)
            p2 = mGas.isen_p( Tt2, pt2, T2 );

            // isentropic state
            real T2s = mGas.isen_T( mTurbineEntry.Tt(), mTurbineEntry.pt(), p2 );
            real h2s = mGas.h( T2s, p2 );

            // new factor
            mDeltaHsRotor = mTurbineEntry.h() - h2s ;

            // correct reaction
            mReaction = mDeltaHsRotor / ( mDeltaHsRotor + mDeltaHsNozzle );

            // David ( 2 )
            real tEtaRotor1 = 0.5 * mW2 * mW2 / ( mDeltaHsRotor + 0.5 * mW1 * mW1 );

            mEtaRotor = ( 1.0 - aOmega ) * tEtaRotor0 + aOmega * tEtaRotor1 ;

            aDeltEta = tEtaRotor1 - tEtaRotor0 ;
            return false ;
        }

//------------------------------------------------------------------------------

        real
        Turbine::finalize_turbine_discharge()
        {
            const real & T2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_T );
            const real & p2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_P );

            real & rho2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_RHO );
            real & h2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_H );
            real & s2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_S ) ;
            real & Ma2 = mTurbineDischarge.value( BELFEM_ENGINE_STATE_MA );

            // new density
            rho2 = mGas.rho( T2, p2 );

            // new enthalpy
            h2 = mGas.h( T2, p2 );

            // new entropy
            s2 = mGas.s( T2, p2 );

            // new mach number
            Ma2 = mC2 / mGas.c( T2, p2 );

            // exit cross section
            mA2 = mDotM / ( mCm2 * rho2 );

            // exit blade height
            mB2 = mA2 / ( mEpsilon * constant::pi * mDm );

            // check value
            return  mB2 - mB1 * mB2B1 ;
        }

//------------------------------------------------------------------------------

        void
        Turbine::compute_pitch_and_chord()
        {
            // compute ratio
            if( ! mPitchChordRatioFlag )
            {
                mPitchChordRatio  = polyval( mPitchChordPoly, std::abs( mBeta2 - mBeta1 ) );
            }

            if( ! mZ2Flag )
            {
                // Aungier 10-17
                mZ2 = std::floor( 12.5 + 0.03 * std::pow( 33.0 - mAlpha1, 2 ));
            }

            // compute pitch, Aungier ( 4-1 )
            mPitch = constant::pi * mDm / mZ2 ;
            mChord = mPitch * mPitchChordRatio ;

            // see NASA SP 8110 ( 18 )
            mBladeOpening = std::sin( mBeta2 - 0.5 * constant::pi ) * mPitch ;

            // edge thickness
            real tt2 = mBladeOpening * ( 1.0 - mCe ) / mCe ;

            mBladeRadius1 = 0.5 * tt2 ;
            mBladeRadius2 = 0.5 * tt2 ;

            // mChordAngle = 0.5 * ( mBeta1 - mBeta2 + constant::pi );

            // Aungier ( 7-3 )
            /*mAxialChord = mChord * std::sin( mChordAngle )
                    + mBladeRadius1 * ( 1.0 - std::sin( mBeta1 ) )
                    + mBladeRadius2 * ( 1.0 - std::sin( mBeta2 ) ) ; */

            mAxialChord = mChord / ( 2.0 * (
                    1. + ( std::sin( mBeta1 ) * std::sin( mBeta2 )
                          + std::cos( mBeta1 ) * std::cos( mBeta2 ) ) ) );

            mChordAngle = std::acos( mAxialChord / mChord );
        }

//------------------------------------------------------------------------------

        void
        Turbine::compute_ventilation_losses()
        {
            // Traupel. p. 437, Eq. 8.4( 33 )

            // Kranz Frei
            // real tCbl = ( 0.045 + 0.58 * mBD ) * std::sin( mBeta2 );

            // kranz eingehüllt
            real tCbl = 0.0095 - 0.55 * std::pow( 0.125 - mBD, 2 ) ;

            mDeltaEtaVentilation = tCbl * ( 1.0 - mEpsilon ) / mEpsilon * mUm * mUm /
                    ( 2.0 * mPhi1 * mYs );

        }

//------------------------------------------------------------------------------

        bool
        Turbine::compute_mixing_and_expansion_losses()
        {
            real zeta2 = 2.0 * ( mTurbineDischarge.ht() - mHt2s ) /
                         ( mW1 * mW1 );

            const real & b2 = mAxialChord;
            const real & s2 = mPitch;

            // arc of admission
            real a = mEpsilon * constant::pi * mDm;

            // help factor
            real Kt1 = 2.0 * mW1 * a / ( mUm * b2 );

            // catch error
            if ( mW1 > 5000 )
            {
                return true ;
            }

            real m ;
            real n ;
            if( Kt1 > 700 )
            {
                real expKt1 = std::exp( Kt1 );
                n = 1. + 4. / (Kt1 * (1. + expKt1)) - 2. / Kt1;
                m = -1. - 1. / Kt1 + 2. / Kt1 *std::log(0.5 * (1. + expKt1) )
                    + 4. * expKt1 / (Kt1 * (1. + expKt1));
            }
            else
            {
                n =  1. - 2. / Kt1 ;
                //m =  1. - ( 3. - std::log( 2 ) ) / Kt1
                m = 1.0 - 2.306852819440055 / Kt1 ;
            }

            real A = mUm * mW1
                    * ( std::cos(mBeta1) - std::cos(mBeta2) )
                    * (1. - 0.5 * s2 / a) * (1. - n);

            real B = 0.5 * zeta2 * mW1 * mW1 * ( 1. - m * (1. - 0.5 * s2 / a ) );

            real C = 0.25 * s2 / a * mUm * mW1 * std::cos(mBeta2) ;

            mDeltaEtaMixingAndExpansion = ( A - B - C ) / mYs ;

            return false ;
        }

//-----------------------------------------------------------------------------

        bool
        Turbine::check_validity()
        {
            // check admission
            if ( mEpsilon > 1.1 || mEpsilon < 0.05 )
            {
                return false ;
            }

            // check haller
            if( mW2 / mW1 < 0.7 )
            {
                return false ;
            }

            // check reaction
            if( mReaction < -0.1 )
            {
                return false ;
            }

            // check compressibility
            if( mPhi2 / mPhi1 > 1.8 )
            {
                return false ;
            }

            return true ;
        }

//------------------------------------------------------------------------------

        bool
        Turbine::compute_nozzle_entry()
        {

            // assume constant cross section
            mA0 = mA1 ;

            // shortcuts
            const real & Tt = mNozzleEntry.value( BELFEM_ENGINE_STATE_TT );
            const real & pt = mNozzleEntry.value( BELFEM_ENGINE_STATE_PT );
            const real & ht = mNozzleEntry.value( BELFEM_ENGINE_STATE_HT );

            real & T = mNozzleEntry.value( BELFEM_ENGINE_STATE_T );
            real & p = mNozzleEntry.value( BELFEM_ENGINE_STATE_P );
            real & h = mNozzleEntry.value( BELFEM_ENGINE_STATE_H );
            real & rho = mNozzleEntry.value( BELFEM_ENGINE_STATE_RHO );
            real & Ma = mNozzleEntry.value( BELFEM_ENGINE_STATE_MA );

            real tPhi0 = 0.0 ;
            mPhi0 = mPhi1 ;

            uint tCount = 0 ;

            T = Tt ;
            p = pt ;
            rho = mGas.rho( T, p );

            while( abs( tPhi0 - mPhi0 ) > 1e-7 )
            {
                tPhi0 = mPhi0 ;
                mCm0 = mDotM / ( rho * mA0 );
                mPhi0 = mCm0 / mUm ;
                h = ht - 0.5 * mCm0 * mCm0 ;
                if( this->check_temperature_range( h ) )
                {
                    return true ;
                }
                T = mGas.T_from_h( h, p );
                p = mGas.isen_p( Tt, pt, T );
                rho = mGas.rho( T, p );

                BELFEM_ERROR( tCount++ < 100, "Too many iterations");
            }
            Ma = mCm0 / mGas.c( T, p );
            return false ;
        }

//------------------------------------------------------------------------------
        void
        Turbine::print()
        {

            std::fprintf( stdout, "power              P         : %8.3f MW\n",
                          ( double ) mP*1e-6 );

            std::fprintf( stdout, "massflow           dotm      : %8.3f kg/s\n",
                          ( double ) mDotM );

            std::fprintf( stdout, "admission          epsilon   : %8.3f\n",
                          ( double ) mEpsilon );

            std::fprintf( stdout, "head rise          psi        : %8.5f \n",
                          ( double ) mPsi );

            std::fprintf( stdout, "flow coefficient   phi0       : %8.5f \n",
                         ( double ) mPhi0 );


            std::fprintf( stdout, "flow coefficient   phi1       : %8.5f \n",
                         ( double ) mPhi1 );


            std::fprintf( stdout, "flow coefficient   phi2       : %8.5f \n",
                         ( double ) mPhi2 );

            std::fprintf( stdout, "\nNozzle Entry\n" );

            std::fprintf( stdout, "total temperature  Tt0        : %8.3f K \n",
                          ( double ) mNozzleEntry.Tt() );
            std::fprintf( stdout, "total pressure     pt0        : %8.3f bar \n",
                          ( double ) mNozzleEntry.pt() * 1e-5 );
            std::fprintf( stdout, "static temperature T0        : %8.3f K \n",
                          ( double ) mNozzleEntry.T() );
            std::fprintf( stdout, "static pressure    p0        : %8.3f bar \n",
                          ( double ) mNozzleEntry.p() * 1e-5 );
            std::fprintf( stdout, "Mach Number        Ma0        : %8.3f\n",
                          ( double ) mNozzleEntry.Ma() );
            std::fprintf( stdout, "Entropy            s0        : %8.3f J/(kg*K)\n",
                          ( double ) mNozzleEntry.s() );

            std::fprintf( stdout, "\nNozzle Discharge / Turbine Entry\n" );

            std::fprintf( stdout, "total temperature  Tt1        : %8.3f K \n",
                          ( double ) mTurbineEntry.Tt() );
            std::fprintf( stdout, "total pressure     pt1        : %8.3f bar \n",
                          ( double ) mTurbineEntry.pt() * 1e-5 );
            std::fprintf( stdout, "static temperature T1        : %8.3f K \n",
                          ( double ) mTurbineEntry.T() );
            std::fprintf( stdout, "static pressure    p1        : %8.3f bar \n",
                          ( double ) mTurbineEntry.p() * 1e-5 );
            std::fprintf( stdout, "Mach Number        Ma1        : %8.3f\n",
                          ( double ) mTurbineEntry.Ma() );
            std::fprintf( stdout, "Entropy            s1        : %8.3f J/(kg*K)\n",
                          ( double ) mTurbineEntry.s() );
            std::fprintf( stdout, "\nTurbine Discharge\n" );
            std::fprintf( stdout, "total temperature  Tt2        : %8.3f K \n",
                          ( double ) mTurbineDischarge.Tt() );
            std::fprintf( stdout, "total pressure     pt2        : %8.3f bar \n",
                          ( double ) mTurbineDischarge.pt() * 1e-5 );
            std::fprintf( stdout, "static temperature T2         : %8.3f K \n",
                          ( double ) mTurbineDischarge.T() );
            std::fprintf( stdout, "static pressure    p2         : %8.3f bar \n",
                          ( double ) mTurbineDischarge.p() * 1e-5 );
            std::fprintf( stdout, "Mach Number        Ma2        : %8.3f\n",
                          ( double ) mTurbineDischarge.Ma() );
            std::fprintf( stdout, "Entropy            s2         : %8.3f J/(kg*K)\n",
                          ( double ) mTurbineDischarge.s() );

            std::fprintf( stdout, "\nEfficiencies\n" );
            std::fprintf( stdout, "isentropic efficnency eta      : %8.3f\n",
                          ( double ) mEta );
            std::fprintf( stdout, "nozzle efficnency eta          : %8.3f\n",
                          ( double ) mEtaNozzle );
            std::fprintf( stdout, "rotor efficnency eta           : %8.3f\n",
                          ( double ) mEtaRotor );
            std::fprintf( stdout, "reaction                       : %8.3f\n",
                          ( double ) mReaction );

            std::fprintf( stdout, "Haller           w2/w1         : %8.3f\n",
                          ( double ) mW2 / mW1 );

            std::fprintf( stdout, "\nGeometry\n" );

            std::fprintf( stdout, "mean diameter    Dm             : %8.3f mm\n",
                          ( double ) mDm * 1000 );

            std::fprintf( stdout, "blade height     b              : %8.3f mm\n",
                          ( double ) mB1 * 1000 );

            std::fprintf( stdout, "number of blades Z2             : %u\n",
                          ( uint ) mZ2 );

            std::fprintf( stdout, "pitch            p              : %8.3f mm\n",
                          ( double ) mPitch * 1000 );

            std::fprintf( stdout, "chord            c              : %8.3f mm\n",
                          ( double ) mChord * 1000 );

            std::fprintf( stdout, "axial chord      a              : %8.3f mm\n",
                          ( double ) mAxialChord * 1000 );

            std::fprintf( stdout, "blade opening    o              : %8.3f mm\n",
                          ( double )  mBladeOpening * 1000 ) ;

            std::fprintf( stdout, "entry opening                   : %8.3f mm\n",
                          ( double ) mBladeEntry * 1000 );

            std::fprintf( stdout, "entry error                   : %8.3f mm\n",
                          ( double ) mBladeEntryError * 1000 );

            std::fprintf( stdout, "\nAngles\n" );

            std::fprintf( stdout, "alpha1                       : %8.3f °\n",
                          ( double ) mAlpha1 / constant::deg );
            std::fprintf( stdout, "alpha2                       : %8.3f °\n",
                          ( double ) mAlpha2 / constant::deg );
            std::fprintf( stdout, "beta1                        : %8.3f °\n",
                          ( double ) mBeta1 / constant::deg );
            std::fprintf( stdout, "beta2                        : %8.3f °\n",
                          ( double ) (  mBeta2 - 0.5 * constant::pi ) / constant::deg );
        }



//------------------------------------------------------------------------------
    }
}
