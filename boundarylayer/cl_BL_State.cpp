//
// Created by Christian Messe on 29.11.19.
//

#include "cl_BL_State.hpp"

namespace belfem
{
    namespace boundarylayer
    {
//------------------------------------------------------------------------------

        State::State( Gas & aGas ) :
            mGas( aGas )
        {

        }

//------------------------------------------------------------------------------

        void
        State::set_wall_temperature( const real & aTw )
        {
            // set the wall temperature
            mTw = aTw;

            // reset anything else wall related
            mRhow     = BELFEM_QUIET_NAN;
            mHw       = BELFEM_QUIET_NAN;
            mCpw      = BELFEM_QUIET_NAN;
            mMuw      = BELFEM_QUIET_NAN;
            mLambdaw  = BELFEM_QUIET_NAN;
            mPrandtlw = BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        void
        State::set_pressure_gradient( const real & adPdX )
        {
            mdPdX = adPdX;
        }

//------------------------------------------------------------------------------

        void
        State::set_loads(
                const real & aTauw,
                const real & aDotQ,
                const real & aHr )
        {

            mTauw      = aTauw;
            mDotQ      = aDotQ;
            mHr        = aHr;
        }

//------------------------------------------------------------------------------
        void
        State::compute( const real & aT, const real & aP, const real & aU )
        {
            mT  = aT;
            mP  = aP;
            mU  = aU;

            mMa = aU / mGas.c( aT, aP );

            // density
            mRho = mGas.rho( aT, aP );

            // enthalpy
            mH      = mGas.h( aT, aP );

            // entropy
            mS      = mGas.s( aT, aP );

            // total enthalpy
            mHt     = mH + 0.5 * mU * mU;

            // dynamic viscosity
            mMu     = mGas.mu( aT, aP );

            // thermal conductivity
            mLambda = mGas.lambda( aT, aP );

            // prandtl number
            mPrandtl = mMu * mGas.cp( aT, aP ) / mLambda;
        }

//------------------------------------------------------------------------------

        void
        State::compute_Cp( const State & aFreestream )
        {
            mCp = 2.0 * ( mP - aFreestream.p() ) /
                ( aFreestream.rho() * aFreestream.u() * aFreestream.u() );
        }

//------------------------------------------------------------------------------

        void
        State::compute_wall_state()
        {
            // density near wall
            mRhow = mGas.rho( mTw, mP );

            // wall enthalpy
            mHw  = mGas.h( mTw, mP );

            // specific heat capacity at the wall
            mCpw = mGas.cp( mTw, mP );

            // viscosity near the wall
            mMuw = mGas.mu( mTw, mP );

            // thermal conductivity near the wall
            mLambdaw = mGas.lambda( mTw, mP );

            // Prandtl number near the wall
            mPrandtlw = mMuw * mCpw / mLambdaw;
        }

//------------------------------------------------------------------------------

    }
}