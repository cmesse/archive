//
// Created by Christian Messe on 25.11.20.
//
#ifndef BELFEM_CL_CH_BOUNDARYLAYER_HPP
#define BELFEM_CL_CH_BOUNDARYLAYER_HPP
#include "typedefs.hpp"
#include "cl_Spline.hpp"
#include "cl_Gas.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "CH_defines.hpp"
#include "CH_Enums.hpp"

namespace belfem
{
    namespace channel
    {
        enum class SigmaRecoveryMode
        {
            vanDriest,
            Petrukov
        };

        class Boundarylayer
        {
            Gas & mGas ;
            BoundaryLayerMethod mMethod ;

            const index_t mNumberOfCells ;
            const index_t mNumberOfNodes ;
            const index_t mCenter ;
            const real    mMeshRatio ;

            //! Kármán Constant
            const real mKarman = 0.41 ;

            //! constant for log-law
            real mBplus  = 5.0 ;

            const real mBplus0     = 5.0 ;   // value if there is no roughness
            const real mRoughConst = 3.4 ;   // for roughness correlation

            // flag if this is axisymmetric, otherwise 2D channel
            bool mIsAxisymmetric = true ;

            // value needed for the spalding function
            real mExpKB ;

            //! coles wake parameter
            real mPi = 0.0 ;

            //! technical roughness = 4.2 * Ra
            real mKtech  = 0 ;

            //! cross section
            real mA  = BELFEM_QUIET_NAN ;

            //! hydraulic diameter
            real mDh = BELFEM_QUIET_NAN ;

            //! averaged temperature
            real mTm ;

            //! averaged velocity
            real mUm ;

            //! averaged density
            real mRhom ;

            //! averaged enthalpy
            real mHm ;

            //! averaged entropy
            real mSm ;

            //! massflow
            real mDotM ;

            //! momentum flow
            real mDotI ;

            //! energy flow
            real mDotH ;

            //! shear velocity
            real mUtau     = BELFEM_QUIET_NAN;

            // non-dimensional grid
            Vector< real > mEta ;

            Matrix< real > mData ;

            // pressure ( must be constant )
            real mP = BELFEM_P_REF ;

            // reynolds number
            real mReDh  = BELFEM_QUIET_NAN ;
            real mReDhw = BELFEM_QUIET_NAN ;

            // mean prandtl number
            real mPrm    = BELFEM_QUIET_NAN ;

            // recovery number
            real mRecovery = BELFEM_QUIET_NAN ;

            // reynolds colburn analogy
            real mSigma    = BELFEM_QUIET_NAN ;

            // for Crocco-Busemann Equation
            cplx mPhi;
            cplx mPsi;
            cplx mChi;
            cplx mAlpha;
            cplx mBeta;
            real mCplus;  // dyplusdy

            uint mNumSplineSteps ;
            real mTmin ;
            real mTmax ;
            real mRhoMax ;

            SpMatrix mHelpMatrix ;
            Vector< real > mWorkTemperature ;
            Vector< real > mWorkVolume ;
            Vector< real > mWorkHeat ;
            Vector< real > mWorkMu ;
            Vector< real > mWorkLambda ;

            Spline * mVolumeSpline = nullptr ;
            Spline * mLambdaSpline = nullptr ;
            Spline * mHeatSpline   = nullptr ;
            Spline * mMuSpline     = nullptr ;

            real mTw1 ;
            real mTw2 ;
            //real mRadius1 ;
            //real mRadius2 ;
            //real mRhat    ;

            // if this switch is set, inflow and geometry are taken from
            // parameters object
            bool mUseParametersAsInput = false ;

            // container for balance constants. Only used my Messe Model
            Vector< real > mBalance ;

            void
            ( Boundarylayer::* mSigmaRecoveryFunction )();

            void
            ( Boundarylayer::* mFrictionFunction )( real & aDotQ, real & aTauw, real & aTr, real & aHr );

            // constant for Bartz equation: ( Dt / rc ) ^0.1
            real mBartzConst = BELFEM_QUIET_NAN ;

            Gas * mAlternateGas = nullptr ;

            /*
            // reference x-position for van Driest ( Nozzle only )
            real mXref = BELFEM_QUIET_NAN ;

            // offset x-position for van Driest ( Nozzle only )
            real mXoff = 0.0 ;

            // pressure gradient for van Driest ( Nozzle only )
            real mdPdX = 0.0 ;

            // turbulent Prandtl number vor van Direst ( Nozzle only )
            const real mPrTinf = 0.85 ;

            // Coles' Wake parameter, but for van Driest ( Nozzle only )
            real mPi_driest = 0.55 ; */

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Boundarylayer( Gas & aGas,
                           const BoundaryLayerMethod aBoundaryLayerMethod,
                           const SigmaRecoveryMode  aSigmaRecoveryMode,
                           const index_t aNumberOfCells=50,
                           const real    aMeshRatio = 1.2 ) ;

//------------------------------------------------------------------------------

            ~Boundarylayer() ;

//------------------------------------------------------------------------------

            void
            set_flow_conditions( const real & aT,
                                 const real & aP,
                                 const real & aU,
                                 const bool   aUpdateLookupTables=true );

//------------------------------------------------------------------------------

            void
            set_flow_conditions( const real & aT,
                                 const real & aP,
                                 const real & aU,
                                 const Matrix< real > & aHeatSplineData,
                                 const Matrix< real > & aViscositySplineData,
                                 const Matrix< real > & aConductivitySplineData );

//------------------------------------------------------------------------------

            void
            set_center_conditions(
                    const real & aT,
                    const real & aU );

//------------------------------------------------------------------------------

            void
            set_wall_temperature( const real & aTwall );

//------------------------------------------------------------------------------

            void
            set_hydraulic_diameter( const real & aDh );

//------------------------------------------------------------------------------

            void
            set_surface_roughness( const real & aRa );

//------------------------------------------------------------------------------

            /**
             * set the parameters needed for the Bartz equation
             */
             void
             set_bartz_geometry_params( const real & aDt, const real & aRc );

//------------------------------------------------------------------------------

            /**
             * if this is set, we don't need the setters for
             *  - set_flow_conditions
             *  - set_wall_temperature
             *  - set_hydraulic_diameter
             * @param aSwitch
             */
            void
            use_input_from_parameters( const bool aSwitch );

//------------------------------------------------------------------------------

            /**
             * returns the computed result as follows
             * @param aParameters
             *
             *  0 : x      - placeholder for reference position ( not written by this function )
             *  1 : A      - placeholder for cross section      ( not written by this function )
             *  2 : D_h    - hydraulic diameter in m
             *
             *  3 : T_m    - mean flow temperature  in K
             *  4 : P_m    - mean flow pressure     in Pa
             *  5 : u_m    - mean flow velocity     in m/2
             *  6 : Ma_m   - mean Mach number       -
             *
             *  7 : h_m    - mean specific enthalpy J/kg
             *  8 : s_m    - mean specific entropy  J/(kg K)
             *  9 : Pr_m   - mean Prandtl number
             *
             * 10 : Re_Dh  - Reynolds number with respect to mean values
             *
             * 11 : T_hat -  temperature at center of flow    in K
             * 12 : u_hat  - velocity at center of flow       in m/s
             *
             * 13 : err_M  - error in mass balance equation (non-dimensional)
             * 14 : err_I  - error in momentum balance equation (non-dimensional)
             * 15 : err_T  - temperature Error (non-dimensional)
             *
             * 16 : T_w    - wall temperature       in K      ( wall 1 )
             *
             * 17 : tau_w  - shear stress           in Pa     ( wall 1 )
             * 18 : dot_q  - heat flux              in W/m^2  ( wall 1 )
             * 19 : h_w    - wall enthalpy            in J/kg ( wall 1 )
             * 20 : Y+     - Y+ of first wall element         ( wall 1 )
             *
             * 21 : T_rec  - recovery temperature in K        ( wall 1 )
             * 22 : h_rec  - recovery enthalpy in J/kg        ( wall 1 )
             * 23 : alpha  - linear heat transfer factor      ( wall 1 )
             *
             * 24 : T_w    - wall temperature       in K      ( wall 2 )
             * 25 : tau_w  - shear stress           in Pa     ( wall 2 )
             * 26 : dot_q  - heat flux              in W/m^2  ( wall 2 )
             * 27 : h_w    - wall enthalpy            in J/kg ( wall 2 )
             * 28 : Y+     - Y+ of first wall element         ( wall 2 )
             *
             * 29 : T_rec  - recovery temperature in K        ( wall 2 )
             * 30 : h_rec  - recovery enthalpy in J/kg        ( wall 2 )
             * 31 : alpha  - linear heat transfer factor      ( wall 2 )
             */
            void
            compute( Vector< real > & aParameters, const bool aUpdateLookupTables=true );

//------------------------------------------------------------------------------

            void
            print();

//------------------------------------------------------------------------------

            /**
             * @return temperature at wall
             */
            const real &
            Tw() const ;

//------------------------------------------------------------------------------

            /**
             * @return density at wall
             */
            const real &
            rho_w() const ;

//------------------------------------------------------------------------------

            /**
             * @return dynamic viscosity at wall
             */
            const real &
            mu_w() const ;

//------------------------------------------------------------------------------

            /**
             * @return the shear stress at the wall
             */
            const real &
            tau_w() const ;

//------------------------------------------------------------------------------

            /**
             * set initial values for tau_w, recovery and reynolds-colburn
             */
            void
            compute_initial_guesses();

//------------------------------------------------------------------------------

            /**
             * expose the volume spline
             */
             Spline *
             volume_spline();

//------------------------------------------------------------------------------

            /**
             * expose the heat spline
             */
            Spline *
            heat_spline();

//------------------------------------------------------------------------------

            /**
             * expose the viscosity spline
             */
            Spline *
            viscosity_spline();

//------------------------------------------------------------------------------

            /**
             * expose the conductivity spline
             */
            Spline *
            conductivity_spline();

//------------------------------------------------------------------------------

            void
            update_lookup_tables();

//------------------------------------------------------------------------------

            /**
             * determines if the linearized heat transfer coefficient
             * refers to Tm or Tr.
             */
            void
            set_sigma_recovery_mode( const SigmaRecoveryMode aMode );

//------------------------------------------------------------------------------

            /**
             * set the friction model of this channel
             */
             void
             set_friction_method( const channel::BoundaryLayerMethod & aMethod );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            real
            compute_outer_step( const real & aPi, Vector< real > & aBalance );

 //------------------------------------------------------------------------------

            void
            compute_inner_step(
                          const real & aThat,
                          const real & aUhat,
                          const real & aPi,
                          Vector< real > & aBalance );

//------------------------------------------------------------------------------

            // create the nondimensional cell grid
            // and allocate the vectors
            void
            make_grid_and_allocate();

//------------------------------------------------------------------------------

            /**
             * adapt parameters to account for roughness
             * also computes B+ constant
             */
            void
            compute_wall_state();

//------------------------------------------------------------------------------

            /**
             * computes u+, u and dudy
             */
             void
             compute_velocity_profile();

//------------------------------------------------------------------------------

            /**
             * computes rho, T, cp, h, mu and lambda profiles
             */
             void
             compute_temperature_profile();

//------------------------------------------------------------------------------

            /**
             * computes mu_T, lambda_T, Pr, Pr_T, Pr_m, tau
             */
             void
             compute_turbulence();

//------------------------------------------------------------------------------

            void
            compute_sigma_recovery() ;

//------------------------------------------------------------------------------

             void
             compute_sigma_recovery_vandriest() ;

//------------------------------------------------------------------------------

             void
             compute_sigma_recovery_petrukov() ;

//------------------------------------------------------------------------------

            real
            compute_heatflux();

//------------------------------------------------------------------------------

            void
            check_balance( Vector< real > & aBalance );

//------------------------------------------------------------------------------

            /**
             * compute the derivative of a field to y
             */
             void
             derive( const uint aSource, const uint aTarget );

//------------------------------------------------------------------------------

            /**
             * integrate a field along y
             */
            void
            integrate( const uint aSource, const uint aTarget, const bool aIsAxisymmetric );

//------------------------------------------------------------------------------

            void
            compute_shear_stress();

//------------------------------------------------------------------------------

            void
            init_lookup_tables();

//------------------------------------------------------------------------------
// fricition functions
//------------------------------------------------------------------------------

            /**
             * calls the boundary layer model
             *
             * @param aDotQ
             * @param aTauw
             * @param aHr
             * @param aTr
             */
            void
            friction_messe( real & aDotQ,
                            real & aTauw,
                            real & aHr,
                            real & aTr );

//------------------------------------------------------------------------------

            /**
             * see Huzel & Huang: Modern Engineering for Design of
             *                     Liquid Rocket Engines.
             *     American Institute of Aeronautics and Astronautics, 1992
             *
             * @param aDotQ
             * @param aTauw
             * @param aHr
             * @param aTr
             */
            void
            friction_bartz( real & aDotQ,
                            real & aTauw,
                            real & aHr,
                            real & aTr );

//------------------------------------------------------------------------------

            /**
             * combines reference temperature method with moody chart
             *
             * @param aDotQ
             * @param aTauw
             * @param aHr
             * @param aTr
             */
            void
            friction_eckert( real & aDotQ,
                             real & aTauw,
                             real & aHr,
                             real & aTr );

//------------------------------------------------------------------------------

            /**
             * model for methane, see DOI: 10.1080/10407782.2015.1080575
             * @param aDotQ
             * @param aTauw
             * @param aHr
             * @param aTr
             */
            void
            friction_pizzarelli( real & aDotQ,
                                 real & aTauw,
                                 real & aHr,
                                 real & aTr );

//------------------------------------------------------------------------------

            /**
             * model for methane, same as the one used in RPA
             * @param aDotQ
             * @param aTauw
             * @param aHr
             * @param aTr
             */
            void
            friction_lebedinsky_kalmykov( real & aDotQ,
                                          real & aTauw,
                                          real & aHr,
                                          real & aTr );

//------------------------------------------------------------------------------

            /**
             * model for hotgas, nozzle only
             * @param aDotQ
             * @param aTauw
             * @param aHr
             * @param aTr
             */
             /*void
             friction_vandriest( real & aDotQ,
                                 real & aTauw,
                                 real & aHr,
                                 real & aTr ); */

//------------------------------------------------------------------------------

            real
            check_u_hat_from_t_hat( const real & aThat );

        };

//------------------------------------------------------------------------------

        inline const real &
        Boundarylayer::Tw() const
        {
            return mData( 0 , BELFEM_CHANNEL_T );
        }

//------------------------------------------------------------------------------

        inline const real &
        Boundarylayer::rho_w() const
        {
            return mData( 0 , BELFEM_CHANNEL_RHO );
        }

//------------------------------------------------------------------------------

        inline const real &
        Boundarylayer::mu_w() const
        {
            return mData( 0 , BELFEM_CHANNEL_MU );
        }

//------------------------------------------------------------------------------

        inline const real &
        Boundarylayer::tau_w() const
        {
            return mData( 0 , BELFEM_CHANNEL_TAU );
        }

//------------------------------------------------------------------------------

        inline Spline *
        Boundarylayer::volume_spline()
        {
            return mVolumeSpline ;
        }

//------------------------------------------------------------------------------

        inline Spline *
        Boundarylayer::heat_spline()
        {
            return mHeatSpline ;
        }

//------------------------------------------------------------------------------

        inline Spline *
        Boundarylayer::viscosity_spline()
        {
            return mMuSpline ;
        }

//------------------------------------------------------------------------------

        inline Spline *
        Boundarylayer::conductivity_spline()
        {
            return mLambdaSpline ;
        }

//------------------------------------------------------------------------------

        inline void
        Boundarylayer::compute_sigma_recovery()
        {
            ( this->*mSigmaRecoveryFunction)();
        }

//------------------------------------------------------------------------------

        /*inline real
        Boundarylayer::transform_u( const real & aUplus )
        {
                  cplx   A     = -mPhi ;
            const cplx & B     = mPsi ;
                  cplx   C     = std::sqrt( A );
            const real & u_hat = mData( mCenter, BELFEM_CHANNEL_U );


            //G = std::log( B + 2. * ( A * X + std::sqrt( A * ( 1. + B * X + A * X^2 ) ) )

            cplx G0 = std::log( B + 2. * ( A + std::sqrt( A ) ) );


            // aUplus = ( G - G0 ) * u_hat / mUtau ;

            cplx G = aUplus * mUtau / u_hat * C ; // + G0 ;

            // expression on log:
            cplx H = std::exp( G ) ;
            cplx K = H - B ;

            real X = std::real( (  0.25 * K - A) / ( A * ( B - K ) ) );

            return X * u_hat ;
        }*/

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_CL_CH_BOUNDARYLAYER_HPP
