//
// Created by Christian Messe on 01.04.21.
//

#ifndef BELFEM_CL_EN_TURBINE_HPP
#define BELFEM_CL_EN_TURBINE_HPP

#include "typedefs.hpp"
#include "cl_EN_State.hpp"
#include "cl_Vector.hpp"
#include "cl_Gas.hpp"

namespace belfem
{
    namespace engine
    {
        class Turbine
        {
            Gas & mGas ;

            State mNozzleEntry ;
            State mTurbineEntry ;
            State mTurbineDischarge ;
            State mTurbineEntryRotating ;
            State mTurbineDischargeRotating ;

            // minimal enthalpy in lookup table
            real mHmin ;

            // maximal enthalpy lookup table
            real mHmax ;

            // shaft speed in RPM
            real mN ;

            // flag telling if N was set by user
            bool mNflag = false ;

            // specific power
            real mY ;

            // specific power at full efficiency
            real mYs ;

            // flag telling if Y was set by user
            bool mYflag = false ;

            // power
            real mP ;
            bool mPflag = false ;

            // mass flow
            real mDotM ;
            real mDotMflag = false ;

            // pressure coefficient
            real mPsi ;
            bool mPsiFlag = false ;

            // ratio of blade to diameters
            real mBD ;

            // flag telling of BD was set by user
            bool mBDflag = false ;


            // flag telling if B was set by user
            bool mBflag = false ;

            // blade height at nozzle entry
            real mB0 ;

            // blade height at turbine entry and nozzle discharge
            real mB1 ;

            // blade height at turbine discharge
            real mB2 ;

            // mead diameter in m
            real mDm ;

            // mean rotational speed in m/s
            real mUm ;

            real mAlpha1 ;

            // discharge angle
            real mAlpha2 = 90.0 * constant::deg ;

            real mBeta1 ;
            real mBeta2 ;

            // flag telling if entry has been set
            bool mEntryFlag = false ;


            // admission factor
            real mEpsilon ;
            bool mEpsilonFlag = false ;

            // flow coefficients
            real mPhi0 ;
            real mPhi1 ;
            real mPhi2 ;

            // cross sections
            real mA0 ;
            real mA1 ;
            real mA2 ;

            real mCm0 ;
            real mCm1 ;
            real mCm2 ;
            real mC1 ;
            real mC2 ;
            real mCu2 ;
            real mCu1 ;
            real mWu1 ;
            real mWu2 ;
            real mW1 ;
            real mW2 ;


            real mDeltaHsNozzle ;
            real mDeltaHsRotor ;
            real mReaction ;

            // flag telling if phi1 has been set
            bool mPhi1Flag ;

            // efficiency of nozzle, as suggested by David
            const real mEtaNozzle = 0.94 ;

            // trailing edge coefficient, see NASA SP 8110, Eq. ( 16 )
            const real mCe = 0.95 ;

            // efficiency for full admission
            real mEtaFullAdmission ;

            // actual efficiency
            real mEta ;

            // efficiency of rotor, computed later
            real mEtaRotorFullAdmission ;
            real mEtaRotor  ;


            // polynomial for rotor efficiency
            Vector< real > mEtaPoly ;

            // polynomial for pitch-to-chord ratio
            Vector< real > mPitchChordPoly ;

            // ratio for b2b1
            real mB2B1 = 1.0 ;

            // number of blades on rotor
            uint mZ2 ;
            bool mZ2Flag = false ;

            // pitch length
            real mPitch ;

            // chord length
            real mChord ;

            // axial chord
            real mAxialChord ;

            // blade opening, see NASA SP 8110 ( 18 )
            real mBladeOpening ;

            real mBladeEntry ;
            real mBladeEntryError ;

            // blade raduu
            real mBladeRadius1 ;
            real mBladeRadius2 ;

            // chord angle
            real mChordAngle ;

            real mDeltaEtaVentilation ;
            real mDeltaEtaMixingAndExpansion ;

            // help Temperature for ideal case
            real mTt2s ;

            // help enthalpy for ideal case
            real mHt2s ;

            int mErrorCode ;

            // user setting, optional
            real mPitchChordRatio ;
            bool mPitchChordRatioFlag = false ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Turbine( Gas & aGas );

//------------------------------------------------------------------------------

            ~Turbine() = default ;

//------------------------------------------------------------------------------

            /**
             * sets the rotational speed in rpm
             */
             void
             set_n( const real & aN );

//------------------------------------------------------------------------------

            /**
             * sets the mass flow
             */
             void
             set_massflow( const real & aDotM );

//------------------------------------------------------------------------------

            /**
             * sets the discharge angle in rad
             */
            void
            set_alpha2( const real & aAlpha2 );

//------------------------------------------------------------------------------

            /**
             * sets the power in W
             */
            void
            set_power( const real & aP );

//------------------------------------------------------------------------------

            /**
             * sets the turbine head in m^2 / s^2
             * @param aY
             */
            void
            set_Y( const real & aY );

//------------------------------------------------------------------------------

            /**
             * sets the number of blades
             * ( uses formula if not set )
             * @param aZ2
             */
            void
            set_Z2( const uint & aZ2 );

//------------------------------------------------------------------------------

            /**
             * sets the pitch-chord ratio
             * ( uses chart if not set )
             * @param aPitchChordRatio
             */
            void
            set_pitch_chord_ratio( const real & aRatio );

//------------------------------------------------------------------------------

            /**
             * sets the turbine head coefficient, US notation
             * @param aPsi
             */
            void
            set_psi( const real & aPsi );

//------------------------------------------------------------------------------

            /**
             * sets the flow coefficient at turbine entry
             * @param aPsi
             */
            void
            set_phi( const real & aPhi1 );

//------------------------------------------------------------------------------

            /**
             * sets the admission factor
             * @param aEpsilon
             */
            void
            set_epsilon( const real & aEpsilon );

//------------------------------------------------------------------------------

            /**
             * sets the blade-to-diameter-ratio
             * @param aBD
             */
            void
            set_bd( const real & aBD );

//------------------------------------------------------------------------------

            /**
             * sets the blade height
             * @param aB
             */
            void
            set_b( const real & aB );

//------------------------------------------------------------------------------

            /**
             * perform the computation
             * returns true if a result was found
             */
             bool
             compute();

//------------------------------------------------------------------------------

             /**
              * reset all member variables
              */
             void
             reset( const int aErrorCode=0 );

//------------------------------------------------------------------------------

            void
            set_entry( const real & aTt, const real & aPt );

            void
            print();


//------------------------------------------------------------------------------

            // return the pressure coefficient
            const real &
            psi() const ;

//------------------------------------------------------------------------------

            // return the flow factor at nozzle entry
            const real & phi0() const ;

//------------------------------------------------------------------------------

            // return the flow factor at nozzle discharge / turbine entry
            const real &
            phi1() const ;

//------------------------------------------------------------------------------

            // return the flow factor at turbine discharge
            const real &
            phi2() const ;

//------------------------------------------------------------------------------

            // return the reaction degree of the turbine
            const real &
            reaction() const ;

//------------------------------------------------------------------------------

            // return the efficiency
            const real &
            eta() const ;

//------------------------------------------------------------------------------

            // return the blade angle
            const real &
            alpha1() const ;

//------------------------------------------------------------------------------

            // return the blade angle
            const real &
            alpha2() const ;

//------------------------------------------------------------------------------

            // return the blade angle
            const real &
            beta1() const ;

//------------------------------------------------------------------------------

            // return the blade angle
            const real &
            beta2() const ;

//------------------------------------------------------------------------------

            // return the number of blades
            const uint &
            Z2() const ;

//------------------------------------------------------------------------------

            // return the pitch_to_chord_ratio
            const real &
            pitch_chord_ratio() const ;

//------------------------------------------------------------------------------

            // return the blade height
            const real &
            b() const ;

//------------------------------------------------------------------------------

            // return the mean diameter
            const real &
            Dm() const ;

//------------------------------------------------------------------------------

            // return the pitch
            const real &
            pitch() const ;

//------------------------------------------------------------------------------

            // return the chord length
            const real &
            chord() const ;

//------------------------------------------------------------------------------

            // return the axial chord length
            const real &
            axialchord() const ;

//------------------------------------------------------------------------------

            // return the Haller criterion
            const real
            haller() const ;

//------------------------------------------------------------------------------

            // return the admission factor
            const real &
            epsilon() const ;

//------------------------------------------------------------------------------

            // error at blade entry between geometry and mass conservation
            const real &
            blade_entry_error() const ;

//------------------------------------------------------------------------------

            const State *
            nozzle_entry() const;

//------------------------------------------------------------------------------

            const State *
            turbine_entry() const;

//------------------------------------------------------------------------------

            const State *
            turbine_discharge() const;

//------------------------------------------------------------------------------

            const State *
            turbine_entry_rotating() const ;


//------------------------------------------------------------------------------

            const State *
            turbine_discharge_rotating() const ;

//------------------------------------------------------------------------------

            const int &
            error_code() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            // make sure that user parameters make sense
            void
            check_sanity();

//------------------------------------------------------------------------------

            /**
             * compute basic geometry parameters
             */
             void
             compute_diameter();

//------------------------------------------------------------------------------

            /**
             * computes the state for the nozzle exit / turbine entry
             * writes the admission factor into aEpsilon
             * returns true if an error has been detected
             */
            bool
            compute_turbine_entry( const real & aMu, const real & aPhi1, real & aEpsilon );

//------------------------------------------------------------------------------

            void
            compute_rotor_efficiency();

//------------------------------------------------------------------------------

            /**
             * returns true if an error has been found
             * @return
             */
            bool
            predict_turbine_discharge();

//------------------------------------------------------------------------------

            bool
            correct_turbine_discharge( real & aDeltaEta, const real aOmega = 0.6 );

//------------------------------------------------------------------------------

            /**
             * returns the geometry error for the blade entry
             */
            real
            finalize_turbine_discharge();

//------------------------------------------------------------------------------

            /**
             *
             * @param aMu  ratio of phi2/phi1
             * @param aDeltaB2  difference between B2 and B1
             * @return     returns true if an error has been found
             */
            bool
            compute_blade_height( const real & aMu, real & aDeltaB2 );

//------------------------------------------------------------------------------

            void
            compute_pitch_and_chord();

//------------------------------------------------------------------------------

            void
            compute_ventilation_losses();

//------------------------------------------------------------------------------

            bool
            compute_mixing_and_expansion_losses();

//------------------------------------------------------------------------------

            /**
             * returns true if an error has been found
             */
            bool
            compute_nozzle_entry();

//------------------------------------------------------------------------------

            bool
            compute_rotating_states();

//------------------------------------------------------------------------------

            void
            compute_blade_entry();

//------------------------------------------------------------------------------

            bool
            check_validity();

//------------------------------------------------------------------------------

            /**
             * returns true if enthalpy is out of range
             */
            bool
            check_temperature_range( const real & aEnthalpy ) const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline void
        Turbine::set_n( const real & aN  )
        {
            mN = aN ;
            mNflag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_alpha2( const real & aAlpha2  )
        {
            mAlpha2 = aAlpha2 ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_massflow( const real & aDotM  )
        {
            mDotM = aDotM ;
            mDotMflag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_Y( const real & aY  )
        {
            mY = aY ;
            mYflag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_power( const real & aP  )
        {
            mP = aP ;
            mPflag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_Z2( const uint & aZ2  )
        {
            mZ2 = aZ2 ;
            mZ2Flag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_pitch_chord_ratio( const real & aRatio )
        {
            mPitchChordRatio = aRatio ;
            mPitchChordRatioFlag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_psi( const real & aPsi  )
        {
            mPsi = aPsi ;
            mPsiFlag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_phi( const real & aPhi1  )
        {
            mPhi1 = aPhi1 ;
            mPhi1Flag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_bd( const real & aBD  )
        {
            mBD = aBD ;
            mBDflag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_b( const real & aB  )
        {
            mB1 = aB ;
            mBflag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        Turbine::set_epsilon(  const real & aEpsilon  )
        {
            mEpsilon = aEpsilon;
            mEpsilonFlag = true ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::psi() const
        {
            return mPsi ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::phi0() const
        {
            return mPhi0 ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::phi1() const
        {
            return mPhi1 ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::phi2() const
        {
            return mPhi2;
        }

//------------------------------------------------------------------------------

        // return the reaction degree of the turbine
        inline const real &
        Turbine::reaction() const
        {
            return mReaction ;
        }

//------------------------------------------------------------------------------

        // return the efficiency
        inline const real &
        Turbine::eta() const
        {
            return mEta ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::alpha1() const
        {
            return mAlpha1;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::alpha2() const
        {
            return mAlpha2;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::beta1() const
        {
            return mBeta1;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::beta2() const
        {
            return mBeta2;
        }

//------------------------------------------------------------------------------

        inline const uint &
        Turbine::Z2() const
        {
            return mZ2 ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::pitch_chord_ratio() const
        {
            return mPitchChordRatio ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::b() const
        {
            return mB1 ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::Dm() const
        {
            return mDm ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::pitch() const
        {
            return mPitch ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::chord() const
        {
            return mChord ;
        }

//------------------------------------------------------------------------------

        // return the axial chord length
        inline const real &
        Turbine::axialchord() const
        {
            return mAxialChord ;
        }

//------------------------------------------------------------------------------

        inline const real
        Turbine::haller() const
        {
            return mW2 / mW1 ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::epsilon() const
        {
            return mEpsilon ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Turbine::blade_entry_error() const
        {
            return mBladeEntryError ;
        }

//------------------------------------------------------------------------------

        inline const State *
        Turbine::nozzle_entry() const
        {
            return & mNozzleEntry ;
        }

//------------------------------------------------------------------------------

        inline const State *
        Turbine::turbine_entry() const
        {
            return & mTurbineEntry ;
        }


//------------------------------------------------------------------------------

        inline const State *
        Turbine::turbine_discharge() const
        {
            return & mTurbineDischarge ;
        }

//------------------------------------------------------------------------------

        inline const State *
        Turbine::turbine_entry_rotating() const
        {
            return & mTurbineEntryRotating ;
        }


//------------------------------------------------------------------------------

        inline const State *
        Turbine::turbine_discharge_rotating() const
        {
            return & mTurbineDischargeRotating ;
        }

//------------------------------------------------------------------------------

        inline bool
        Turbine::check_temperature_range( const real & aEnthalpy ) const
        {
            return aEnthalpy < mHmin || aEnthalpy > mHmax ;
        }

//------------------------------------------------------------------------------

        inline const int &
        Turbine::error_code() const
        {
            return mErrorCode ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_EN_TURBINE_HPP
