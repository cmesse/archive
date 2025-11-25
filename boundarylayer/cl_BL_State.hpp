//
// Created by Christian Messe on 29.11.19.
//

#ifndef BELFEM_CL_BL_STATE_HPP
#define BELFEM_CL_BL_STATE_HPP

#include "typedefs.hpp"
#include "cl_Gas.hpp"

namespace belfem
{
    namespace boundarylayer
    {
        class State
        {
            Gas & mGas;

            real mT         = BELFEM_QUIET_NAN;
            real mP         = BELFEM_QUIET_NAN;
            real mU         = BELFEM_QUIET_NAN;
            real mMa        = BELFEM_QUIET_NAN;

            real mRho       = BELFEM_QUIET_NAN;

            real mH         = BELFEM_QUIET_NAN;
            real mS         = BELFEM_QUIET_NAN;

            real mHt        = BELFEM_QUIET_NAN;

            real mGamma     = BELFEM_QUIET_NAN;
            real mMu        = BELFEM_QUIET_NAN;
            real mLambda    = BELFEM_QUIET_NAN;

            real mPrandtl   = BELFEM_QUIET_NAN;

            // Pressure coefficient,
            // not specific heat capacity!
            real mCp        = BELFEM_QUIET_NAN;

            real mTw        = BELFEM_QUIET_NAN;
            real mRhow      = BELFEM_QUIET_NAN;
            real mHw        = BELFEM_QUIET_NAN;
            real mCpw       = BELFEM_QUIET_NAN; // specific heat capacity at the wall
            real mMuw       = BELFEM_QUIET_NAN;
            real mLambdaw   = BELFEM_QUIET_NAN;
            real mPrandtlw  = BELFEM_QUIET_NAN;

            real mDotQ      = BELFEM_QUIET_NAN;
            real mTauw      = BELFEM_QUIET_NAN;
            real mHr        = BELFEM_QUIET_NAN;

            real mdPdX      = 0.0 ;
            real mdUdX      = 0.0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            State( Gas & aGas );

//------------------------------------------------------------------------------

            ~State() = default;

//------------------------------------------------------------------------------

            /**
             * expose the gas pointer of the state
             */
            Gas &
            gas();

//------------------------------------------------------------------------------

            void
            compute( const real & aT, const real & aP, const real & aU );

//------------------------------------------------------------------------------

            /**
             * computes the pressure coefficient
             */
             void
             compute_Cp( const State & aFreestream );

//------------------------------------------------------------------------------

            /**
             * sets the wall temperature and resets all wall related
             * thermal properties. Must be called before compute_wall_stata
             */
            void
            set_wall_temperature( const real & aTw );

//------------------------------------------------------------------------------

            void
            set_pressure_gradient( const real & adPdX );

//------------------------------------------------------------------------------

            void
            set_loads(
                    const real & aTauw,
                    const real & aDotQ,
                    const real & aHr );

//------------------------------------------------------------------------------

            void
            compute_wall_state();

//------------------------------------------------------------------------------

            const real &
            T() const;

//------------------------------------------------------------------------------

            const real &
            p() const;

//------------------------------------------------------------------------------

            const real &
            rho() const;

//------------------------------------------------------------------------------

            const real &
            u() const;

//------------------------------------------------------------------------------

            const real &
            Ma() const;

//------------------------------------------------------------------------------

            const real &
            h() const;

//------------------------------------------------------------------------------

            const real &
            s() const;

//------------------------------------------------------------------------------

            const real &
            ht() const;

//------------------------------------------------------------------------------

            const real &
            gamma() const;

//------------------------------------------------------------------------------

            const real &
            mu() const;

//------------------------------------------------------------------------------

            const real &
            lambda() const;

//------------------------------------------------------------------------------

            const real &
            Pr() const;

//------------------------------------------------------------------------------

            const real &
            Pr_w() const;

//------------------------------------------------------------------------------

            const real &
            Tw() const;

//------------------------------------------------------------------------------

            const real &
            rho_w() const;

//------------------------------------------------------------------------------

            /**
              * enthalpy at the wall
              * @return
              */
            const real &
            hw() const;

//------------------------------------------------------------------------------

            /**
             * specific heat capacity at the wall
             * @return
             */
            const real &
            cpw() const;

//------------------------------------------------------------------------------

            const real &
            mu_w() const;

//------------------------------------------------------------------------------

            const real &
            lambda_w() const;

//------------------------------------------------------------------------------

            const real &
            dot_q() const;

//------------------------------------------------------------------------------

            const real &
            tau_w() const;

//------------------------------------------------------------------------------

            const real &
            hr() const;

//------------------------------------------------------------------------------

            // the pressure gradient
            const real &
            dpdx() const;

//------------------------------------------------------------------------------

            // the velocity gradient
            const real &
            dudx() const;

//------------------------------------------------------------------------------

            /**
             * the pressure coefficient
             * @return
             */
            const real &
            Cp() const;

//------------------------------------------------------------------------------

            void
            set_Cp( const real & aCp );

//------------------------------------------------------------------------------

            void
            set_dpdx( const real & adpdx );

//------------------------------------------------------------------------------

            void
            set_dudx( const real & adudx );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline Gas &
        State::gas()
        {
            return mGas;
        }

//------------------------------------------------------------------------------
        inline const real &
        State::T() const
        {
            return mT;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::p() const
        {
            return mP;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::rho() const
        {
            return mRho;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::u() const
        {
            return mU;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::Ma() const
        {
            return mMa;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::h() const
        {
            return mH;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::s() const
        {
            return mS;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::ht() const
        {
            return mHt;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::gamma() const
        {
            return mGamma;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::mu() const
        {
            return mMu;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::lambda() const
        {
            return mLambda;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::Pr() const
        {
            return mPrandtl;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::Pr_w() const
        {
            return mPrandtlw;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::dpdx() const
        {
            return mdPdX;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::dudx() const
        {
            return mdUdX;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::Cp() const
        {
            return mCp;
        }

//------------------------------------------------------------------------------

        inline void
        State::set_Cp( const real & aCp )
        {
            mCp = aCp;
        }

//------------------------------------------------------------------------------

        inline  void
        State::set_dpdx( const real & adpdx )
        {
            mdPdX = adpdx;
        }

//------------------------------------------------------------------------------

        inline  void
        State::set_dudx( const real & adudx )
        {
            mdUdX = adudx;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::Tw() const
        {
            return mTw;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::rho_w() const
        {
            return mRhow;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::hw() const
        {
            return mHw;
        }

// --------------------------------------------------------------------

        inline const real &
        State::cpw() const
        {
            return mCpw;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::mu_w() const
        {
            return mMuw;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::lambda_w() const
        {
            return mLambdaw;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::dot_q() const
        {
            return mDotQ;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::tau_w() const
        {
            return mTauw;
        }

//------------------------------------------------------------------------------

        inline const real &
        State::hr() const
        {
            return mHr;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_BL_STATE_HPP
