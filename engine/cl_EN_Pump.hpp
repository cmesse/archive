//
// Created by Christian Messe on 19.01.21.
//

#ifndef BELFEM_CL_EN_PUMP_HPP
#define BELFEM_CL_EN_PUMP_HPP

#include "typedefs.hpp"
#include "constants.hpp"
#include "units.hpp"
#include "cl_EN_State.hpp"
#include "cl_Vector.hpp"
#include "cl_EN_StepanoffChart.hpp"
namespace belfem
{
    namespace engine
    {
        class Pump
        {
            Gas & mGas ;

            const bool mIsHydrogen ;

            StepanoffChart mStepanoff ;

            // state at entry
            State mEntry ;

            // state at exit if isotropic
            State mExitIsotropic ;

            // state at exit
            State mExit ;

            // exit diameter in m
            real mD2a = BELFEM_QUIET_NAN ;

            // flag telling if D2a was set by user
            bool mD2aFlag = false ;

            // mass flux in kg/s
            real mDotM = BELFEM_QUIET_NAN ;

            // rotational speed in RPM
            real mNRPM = BELFEM_QUIET_NAN ;

            // rotational speed in 1/s
            real mN = BELFEM_QUIET_NAN ;

            // pressure difference in Pa
            real mDeltaP = BELFEM_QUIET_NAN ;


            // flag telling if set_deltap was called
            bool mDeltaPFlag = false;

            // flag telling if set_pt2was called
            bool mPt2Flag = false;

            // total discharge pressure in Pa, only used if value was set by user
            real mPt2 = BELFEM_QUIET_NAN ;

            // volume flux in m^3/s
            real mDotV =  BELFEM_QUIET_NAN ;

            // volume flux in gal / min
            real mQ =  BELFEM_QUIET_NAN ;

            // isentropic enthalpy rise
            real mY = BELFEM_QUIET_NAN ;

            // isentropic enthalpy rise
            real mYth = BELFEM_QUIET_NAN ;

            // pump head in m
            real mH =  BELFEM_QUIET_NAN ;

            // pump head in ft
            real mHft =  BELFEM_QUIET_NAN ;

            // rotional specific speeds
            real mNs =  BELFEM_QUIET_NAN ;
            real mNq =  BELFEM_QUIET_NAN ;

            // specific diameter in us units
            real mDs = BELFEM_QUIET_NAN ;

            // exit diameter in m
            real mD2i = BELFEM_QUIET_NAN ;

            // flag telling is mD2 is equal to mD2a ;
            bool mD2Flag = false ;

            // mean exit diameter in m
            real mD2m = BELFEM_QUIET_NAN ;

            // mean entry diameter in m
            real mD1m = BELFEM_QUIET_NAN ;

            // pressure coefficient in US notation
            real mPsi = BELFEM_QUIET_NAN ;

            // flag telling of mPsi was set by user
            bool mPsiFlag = false ;

            // hub diameter
            real mDN =  0.012 ;

            // flag telling if DN was set by user
            bool mDNFlag = false ;

            // shaft diameter
            real mDw =  BELFEM_QUIET_NAN ;

            // inner diameter
            real mD1i = BELFEM_QUIET_NAN ;

            // tip diameterw
            real mD1a =  BELFEM_QUIET_NAN ;

            // flag telling if D1a was set by user
            bool mD1aFlag = false ;


            // radial entry velocity in m/s
            real mU1 =  BELFEM_QUIET_NAN ;

            // radial entry velocity in m/s
            real mW1 =  BELFEM_QUIET_NAN ;

            // angle in deg
            real mBeta1 = BELFEM_QUIET_NAN ;

            // radial discharge velocity in m/s
            real mU2 =  BELFEM_QUIET_NAN ;

            // suction cross section
            real mAs =  BELFEM_QUIET_NAN ;

            // suction velocity
            real mCs =  BELFEM_QUIET_NAN ;

            // Axial Velocity at entrance
            real mC1 =  BELFEM_QUIET_NAN ;

            // NPSH constant, depends on fluid
            real mKnpsh = BELFEM_QUIET_NAN ;

            // suppression head factor
            real mBetas = BELFEM_QUIET_NAN ;

            // required npsh
            real mNPSHr = BELFEM_QUIET_NAN ;

            // available npsh
            real mNPSHa = BELFEM_QUIET_NAN ;

            // safety margin for NPSH
            real mNPSHs = 6.0 ;

            // Wesche, Table 3.5 p 97 for R 9 ( Rz40 )
            real mAhydro = 0.99 ;

            // number of blades in inducer
            real mZ1 = 4.0 ;
            bool mZ1Flag = false ;

            // number of blades in impeller
            real mZ2 = BELFEM_QUIET_NAN ;
            bool mZ2Flag = false ;

            // tip flow coefficient
            real mPhi1a = BELFEM_QUIET_NAN ;

            real mPhi2 = BELFEM_QUIET_NAN ;

            // inlet factor
            real mEpsilon = BELFEM_QUIET_NAN ;

            // volumetric efficiency
            real mEtaV  =  BELFEM_QUIET_NAN ;

            // hydraulic efficiency
            real mEtaH = BELFEM_QUIET_NAN ;

            // mechanical efficiency
            real mEtaM = BELFEM_QUIET_NAN ;

            // isentropic efficiency
            real mEtaI = BELFEM_QUIET_NAN ;

            // total efficiency
            real mEta  = BELFEM_QUIET_NAN ;

            // friction power
            real mPr = BELFEM_QUIET_NAN ;

            // hydraulic power
            real mPh = BELFEM_QUIET_NAN ;

            // total power
            real mP = BELFEM_QUIET_NAN ;

            // mechanical power
            real mPm = BELFEM_QUIET_NAN ;

            // inner power
            real mPs = BELFEM_QUIET_NAN ;

            // discharge angle in deg
            real mBeta2 = BELFEM_QUIET_NAN ;

            // flag telling of beta2 was set
            bool mBeta2Flag = false ;

            // minderleistung
            real mMu = BELFEM_QUIET_NAN ;

            real mYthinf = BELFEM_QUIET_NAN ;

            real mC2u =  BELFEM_QUIET_NAN ;

            real mC2m =  BELFEM_QUIET_NAN ;

            real mW2 =  BELFEM_QUIET_NAN ;

            // blade thickness in m
            real mS2 = BELFEM_QUIET_NAN ;


            real mB2 =  BELFEM_QUIET_NAN ;


            // length of inducer in m
            real mInducerLength =  BELFEM_QUIET_NAN ;

            // length of impeller in m
            real mImpellerLength =  BELFEM_QUIET_NAN ;

            // value of haller criterion
            real mHaller = BELFEM_QUIET_NAN ;
            real mHallerFlag = false ;

            // suggestions for Z1
            Vector< real > mZ1Table = { 0, 1, 2, 3, 4, 5, 3, 4,
                                        4, 3, 5, 4, 4, 4, 5, 5,
                                        4, 6, 6, 5, 5, 3, 4, 4,
                                        4, 5, 3, 3, 4, 5, 5, 4 };

            // suggestions for Z2
            Vector< real > mZ2Table = { 0,  1,  2,  3,  4,  5,  6,  8,
                                        8,   9, 10, 12, 12, 12, 15, 15,
                                        16, 18, 18, 20, 20, 21, 24, 24,
                                        24, 25, 27, 27, 28, 30, 30, 32 };

// ------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Pump( Gas & aGas ) ;

//------------------------------------------------------------------------------

            ~Pump() = default ;

//------------------------------------------------------------------------------

            /**
             * set the mass flux in kg/s
             */
             void
             set_mass_flux( const real & aMassFlux ) ;

//------------------------------------------------------------------------------

           /**
            *  set the total conditions at the entry
            * @param aTt  : Total Temperature in K
            * @param aPt  : Total Pressure in Pa
            */
            void
            set_entry( const real & aTt, const real & aPt );

//------------------------------------------------------------------------------

            /**
             * set the rotational speed in RPM
             */
            void
            set_nrpm( const real & aNRPM ) ;

//------------------------------------------------------------------------------

            /**
             * set the pressure difference in Pa
             */
            void
            set_deltap( const real & aDeltaP ) ;

//------------------------------------------------------------------------------

            /**
             * set the exit pressure in Pa
             */
            void
            set_pt2( const real & aPt2 ) ;

//------------------------------------------------------------------------------

            /**
             * prescribe the discharge diameter
             * @param aD2a     diameter in m
             * @param aFixD2   flag telling if we enforce D2=D2m=D2a
             */
             void
             set_D2a( const real & aD2a, const bool aFixD2 = true );

//------------------------------------------------------------------------------

            /**
             * prescribe the pressure coefficient
             *
             * @param aPsi     pressure coefficient in US notation
             *
             * automatically sets D2=D2m=D2a
             */
            void
            set_psi( const real & aPsi );

//------------------------------------------------------------------------------

            /**
             * prescribe the hub diameter
             */
            void
            set_DN( const real & aDN );

//------------------------------------------------------------------------------

            /**
             * prescribe the suction diameter
             */
            void
            set_D1a( const real & aD1a );

//------------------------------------------------------------------------------

            /**
             * prescribe the number of blades for the inducer
             */
            void
            set_z1( const uint & aZ1 );

//------------------------------------------------------------------------------

            /**
             * prescribe the number of blades for the impeller
             */
            void
            set_z2( const uint & aZ2 );

//------------------------------------------------------------------------------

            /**
             * prescribe the blade thickness
             */
            void
            set_s2( const real & aS2 );

//------------------------------------------------------------------------------

            /**
             * prescribe the discharge angle
             * @param aBeta2  : discharge angle in deg
             */
            void
            set_beta2( const uint & aBeta2);

//------------------------------------------------------------------------------

            void
            compute();

//------------------------------------------------------------------------------

            void
            print();

//------------------------------------------------------------------------------

            /**
             * find a value for beta2 to fullfill the haller criterion
             */
             void
             set_haller( const real aHaller=0.7 );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            compute_isotropic_exit() ;

//------------------------------------------------------------------------------

            void
            compute_rotation_specific_speeds();

//------------------------------------------------------------------------------

            void
            compute_exit_diameters();

//------------------------------------------------------------------------------

            void
            compute_entry_diameters();

//------------------------------------------------------------------------------

            void
            compute_npsh();

//------------------------------------------------------------------------------

            real
            compute_npshr( const real & aD1a );

//------------------------------------------------------------------------------

            void
            compute_efficiencies();

//------------------------------------------------------------------------------

            void
            compute_triangles();

//------------------------------------------------------------------------------

            void
            compute_static( State & aState, const real & aA ) ;

//------------------------------------------------------------------------------

            void
            compute_axial_lengths();

//-------------------------------------------------------------------------

            /**
             * @param aBeta2 to be passed in rad
             * @return w2/w1
             */
            real
            compute_haller( const real & aBeta2 );

//-------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_EN_PUMP_HPP
