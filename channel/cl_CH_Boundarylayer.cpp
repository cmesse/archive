//
// Created by Christian Messe on 25.11.20.
//
#include "cl_Logger.hpp"
#include "cl_CH_Boundarylayer.hpp"

#include "fn_gesv.hpp"
#include "fn_trans.hpp"

#include "fn_Mesh_ratio.hpp"
#include "cl_GM_Helmholtz.hpp"
#include "fn_BL_Spalding.hpp"
#include "fn_BL_g_plus.hpp"
#include "fn_BL_Kays_Crawford.hpp"
#include "fn_linspace.hpp"
#include "cl_GT_RefGas.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_polyval.hpp"
#include "en_GM_GasModel.hpp"
#include "fn_BL_ReferenceTempertaure.hpp"
#include "fn_BL_Moody.hpp"
#include "fn_BL_cf_flatplate_inc_turbulent.hpp"

namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        Boundarylayer::Boundarylayer(
                Gas & aGas,
                const BoundaryLayerMethod aBoundaryLayerMethod,
                const SigmaRecoveryMode aSigmaRecoveryMode,
                const index_t aNumberOfCells,
                const real    aMeshRatio) :
            mGas( aGas ),
            mMethod( aBoundaryLayerMethod ),
            mNumberOfCells( aNumberOfCells ),
            mNumberOfNodes( 2* aNumberOfCells + 1 ),
            mCenter( 2 * aNumberOfCells ),
            mMeshRatio( aMeshRatio )
        {
            // check if we need a help model
            if( aGas.gas_model() == GasModel::HELMHOLTZ )
            {
                mAlternateGas = new Gas(  aGas.component( 0 )->label(), GasModel::SRK );
            }

            // allocate vector for balance
            mBalance.set_size( 3, 0.0 );

            // link functions
            this->set_friction_method( aBoundaryLayerMethod );
            this->set_sigma_recovery_mode( aSigmaRecoveryMode );

            this->make_grid_and_allocate();
            this->init_lookup_tables();

            // set a default value for the wall temperature
            mTw1 = 300.0 ;
            mTw2 = 300.0 ;
            this->set_wall_temperature( 300.0 );
        }

//------------------------------------------------------------------------------

        Boundarylayer::~Boundarylayer()
        {
            if( mAlternateGas != nullptr )
            {
                delete mAlternateGas ;
            }

            delete mVolumeSpline ;
            delete mHeatSpline ;
            delete mMuSpline ;
            delete mLambdaSpline ;
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::make_grid_and_allocate()
        {
            ratio_adx2( mMeshRatio, 1.0, mNumberOfCells, mEta );
            mData.set_size( mNumberOfNodes, BELFEM_CHANNEL_N, BELFEM_QUIET_NAN  );
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::set_flow_conditions(
                const real & aT,
                const real & aP,
                const real & aU,
                const bool   aUpdateLookupTables )
        {
            mTm = std::max( std::min( aT , mTmax ), mTmin );
            mP  = aP ;
            mUm = aU ;

            mRhom = mGas.rho( aT, aP );
            mHm   = mGas.h( aT , aP );
            mSm   = mGas.s( aT , aP );

            // compute Prandtl number
            mPrm = mGas.Pr( aT, aP );

            if( aUpdateLookupTables )
            {
                this->update_lookup_tables();
            }
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::set_flow_conditions( const real & aT,
                             const real & aP,
                             const real & aU,
                             const Matrix< real > & aHeatSplineData,
                             const Matrix< real > & aViscositySplineData,
                             const Matrix< real > & aConductivitySplineData )
        {
            mTm = aT ;
            mP  = aP ;
            mUm = aU ;

            mRhom = mGas.rho( aT, aP );
            mHm   = mGas.h( aT, aP );
            mSm   = mGas.s( aT, aP );

            // write data into splines
            mHeatSpline->matrix_data() = aHeatSplineData ;
            mMuSpline->matrix_data() = aViscositySplineData ;
            mLambdaSpline->matrix_data() = aConductivitySplineData ;

            // compute Prandtl number
            mPrm = mHeatSpline->deval( mTm ) * mMuSpline->eval( mTm ) /
                    mLambdaSpline->eval( mTm );

        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::set_center_conditions(
                const real & aT,
                const real & aU )
        {
            mData( mCenter, BELFEM_CHANNEL_T )   = std::max( aT , mTmin );
            mData( mCenter, BELFEM_CHANNEL_U )   = aU ;

            if( mGas.is_idgas() )
            {
                mData( mCenter, BELFEM_CHANNEL_RHO ) = mGas.rho( aT, mP );
            }
            else
            {
                mData( mCenter, BELFEM_CHANNEL_RHO ) = 1.0 / mVolumeSpline->eval( mData( mCenter, BELFEM_CHANNEL_T ) ) ;
            }

            mData( mCenter, BELFEM_CHANNEL_H )   = mHeatSpline->eval( mData( mCenter, BELFEM_CHANNEL_T ) );
        }


//------------------------------------------------------------------------------

        void
        Boundarylayer::set_wall_temperature( const real & aTwall )
        {
            real tT = std::max( std::min( aTwall , mTmax  ), mTmin );
            mTw1 = tT ;
            mTw2 = tT;
            mData( 0 , BELFEM_CHANNEL_T )     = tT ;
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::set_hydraulic_diameter( const real & aDh )
        {
            mDh = aDh ;
            mA  = 0.25 * constant::pi * aDh * aDh ;

            real tRh = 0.5 * mDh ;

            // radius for integration
            for( uint k=0; k<mNumberOfNodes; ++k )
            {
                mData( k, BELFEM_CHANNEL_R ) = ( 1.0 - mEta( k ) ) * tRh ;
            }

            // wall distance
            for( uint k=0; k<mNumberOfNodes; ++k )
            {
                mData( k, BELFEM_CHANNEL_Y ) = mEta( k ) * tRh ;
            }
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::set_surface_roughness( const real & aRa )
        {
            mKtech = 4.2 * aRa ;
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::set_bartz_geometry_params( const real & aDt, const real & aRc )
        {
            mBartzConst = std::pow( aDt / aRc, 0.1 );
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::use_input_from_parameters( const bool aSwitch )
        {
            mUseParametersAsInput = aSwitch ;
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::compute( Vector< real > & aParameters, const bool aUpdateLookupTables )
        {
            BELFEM_ASSERT( aParameters.length() >= 24,
                          "Parameter Vector must be allocated with at least 24 entries" );

            if( mUseParametersAsInput )
            {
                this->set_hydraulic_diameter( aParameters( BELFEM_CHANNEL_DH ) );

                this->set_flow_conditions(
                        aParameters( BELFEM_CHANNEL_TM ),
                        aParameters( BELFEM_CHANNEL_PM ),
                        aParameters( BELFEM_CHANNEL_UM ),
                        aUpdateLookupTables );

                this->set_wall_temperature( aParameters( BELFEM_CHANNEL_TW1 ) );
            }

            // Averaged Reynolds Number
            mReDh = mRhom * mUm * mDh / mMuSpline->eval( mTm );

            // Wall Reynolds Number
            mReDhw = mGas.rho( this->Tw(), mP ) * mUm * mDh / mMuSpline->eval( this->Tw() );

            real h_r ;
            real dot_q ;
            real T_r ;
            real Ma_m = mUm / mGas.c( mTm, mP );

            const real & T_hat = mData( mCenter, BELFEM_CHANNEL_T ) ;
            const real & u_hat = mData( mCenter, BELFEM_CHANNEL_U ) ;
            //const real & h_hat = mData( mCenter, BELFEM_CHANNEL_H );
            const real & T_w   = mData( 0, BELFEM_CHANNEL_T );
            const real & h_w   = mData( 0, BELFEM_CHANNEL_H );
                  real & tau_w = mData( 0, BELFEM_CHANNEL_TAU );

            // this->compute_initial_guesses();

            // only needed for nozzle
            // mXref = aParameters( BELFEM_CHANNEL_S );

            // todo: this does not work. Remove!
            // mdPdX = aParameters( BELFEM_CHANNEL_DPDS );

            ( this->*mFrictionFunction )( dot_q, tau_w, h_r, T_r );

            aParameters(  2 ) = mDh ;

            aParameters(  3 ) = mTm ;
            aParameters(  4 ) = mP  ;
            aParameters(  5 ) = mUm ;
            aParameters(  6 ) = Ma_m ;

            aParameters(  7 ) = mHm  ;
            aParameters(  8 ) = mSm ;

            aParameters(  9 ) = mGas.Pr( mTm, mP );

            aParameters(  10 ) = mReDh ;

            aParameters( 11 ) = T_hat ;
            aParameters( 12 ) = u_hat ;

            if( mMethod == BoundaryLayerMethod::Messe )
            {
                aParameters( 13 ) = mBalance( 0 );  // massflow error
                aParameters( 14 ) = mBalance( 1 );  // momentum error
                aParameters( 15 ) = mBalance( 2 );  // energy error
            }
            /*else
            {
                // 13: reference position
                // 14: pressure gradient
                // 15: Wake Parameter
                aParameters( 15 ) = mPi_driest;
            }*/

            aParameters( 16 ) = T_w ;
            aParameters( 17 ) = tau_w ;
            aParameters( 18 ) = dot_q ;
            aParameters( 19 ) = h_w ;

            aParameters( 20 ) = mData( 1, BELFEM_CHANNEL_YPLUS );

            aParameters( 21 ) = T_r ;
            aParameters( 22 ) = h_r ;

            // linearized factor in reference to either Tm or Tw
            aParameters( 23 ) = dot_q / ( T_r - T_w );

        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::print()
        {
            for( uint k=0; k<mNumberOfNodes; ++k )
            {
                std::cout << mData( k, BELFEM_CHANNEL_YPLUS ) << " "
                          << mData( k, BELFEM_CHANNEL_UPLUS ) << " "
                          << mData( k, BELFEM_CHANNEL_Y ) << " "
                          << mData( k, BELFEM_CHANNEL_U ) << " "
                          << mData( k, BELFEM_CHANNEL_T ) << std::endl ;
            }

            std::cout << std::endl ;

            std::cout << "phi " << mPhi << std::endl ;
            std::cout << "psi " << mPsi << std::endl ;
            std::cout << "chi " << mChi << std::endl ;

            std::cout << "Y+1 " << mData( 1, BELFEM_CHANNEL_YPLUS ) << std::endl ;

            std::cout << "T_hat " << mData( mCenter, BELFEM_CHANNEL_T ) << std::endl ;
            std::cout << "u_hat " << mData( mCenter, BELFEM_CHANNEL_U ) << std::endl ;

            std::cout << "tau_w " << this->tau_w() << std::endl ;

            std::cout << "sigma " << mSigma<< std::endl ;
            std::cout << "r " << mRecovery << std::endl ;
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::set_sigma_recovery_mode( const SigmaRecoveryMode aMode )
        {
            switch( aMode )
            {
                case( SigmaRecoveryMode::vanDriest ) :
                {
                    mSigmaRecoveryFunction
                        = & Boundarylayer::compute_sigma_recovery_vandriest ;
                    break ;
                }
                case( SigmaRecoveryMode::Petrukov ) :
                {
                    mSigmaRecoveryFunction
                            = & Boundarylayer::compute_sigma_recovery_petrukov ;
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown Mode");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::set_friction_method( const channel::BoundaryLayerMethod & aMethod )
        {
            // remember method
            mMethod = aMethod ;

            switch( aMethod )
            {
                case( channel::BoundaryLayerMethod::Messe ) :
                {
                    mFrictionFunction
                        = & Boundarylayer::friction_messe ;
                    break ;
                }
                case( channel::BoundaryLayerMethod::Bartz ) :
                {
                    mFrictionFunction
                        = & Boundarylayer::friction_bartz ;
                    break ;
                }
                case( channel::BoundaryLayerMethod::Eckert ) :
                {
                    mFrictionFunction
                        = & Boundarylayer::friction_eckert ;
                    break ;
                }
                case( channel::BoundaryLayerMethod::LebedinskyKalmykov ) :
                {
                    mFrictionFunction
                        = & Boundarylayer::friction_lebedinsky_kalmykov ;
                    break ;
                }
                case( channel::BoundaryLayerMethod::Pizzarelli ) :
                {
                    mFrictionFunction
                        = & Boundarylayer::friction_pizzarelli ;
                    break ;
                }
                /*case( channel::BoundaryLayerMethod::vanDriest ) :
                {
                    mFrictionFunction
                        = & Boundarylayer::friction_vandriest ;
                    break ;
                }*/
                default:
                {
                    BELFEM_ERROR( false, "Invalid boundary layer method");
                }
            }
        }


//------------------------------------------------------------------------------

        real
        Boundarylayer::compute_outer_step( const real & aPi, Vector< real > & aBalance  )
        {
            // Relaxation Factors
            real tOmega = 0.5 ;
            real tOmegaT ;
            real tOmegaU ;
            real tAlpha ;

            // Data Containers
            Matrix< real > tJ( 2, 2 );
            Vector< real > tA( 2 );
            Vector< real > tB( 2 );
            Vector< real > tC( 2 );
            Vector< real > tD( 2 );
            Vector< real > tX( 2 );
            Vector< real > tY( 3 );
            Vector< int > tP( 2 );

            real tThat = mData( mCenter, BELFEM_CHANNEL_T ) ;
            real tUhat = mData( mCenter, BELFEM_CHANNEL_U );

            real tDeltaT = mTm * 0.0001 ;
            real tDeltaU = mUm * 0.0001 ;


            real tErr = BELFEM_REAL_MAX ;
            uint tCount = 0 ;

            //real tTmax = this->Tw() > mTm ? this->Tw() : mTmax ;
            real tTmax = mTmax ;

            while ( tErr > 1e-9 )
            {

                this->compute_inner_step( tThat - tDeltaT, tUhat, aPi, tA );
                this->compute_inner_step( tThat + tDeltaT, tUhat, aPi, tB );
                this->compute_inner_step( tThat, tUhat-tDeltaU, aPi, tC );
                this->compute_inner_step( tThat, tUhat+tDeltaU, aPi, tD );

                this->compute_inner_step( tThat, tUhat, aPi, tX );
                this->compute_sigma_recovery();

                for( uint i=0; i<2; ++i )
                {
                    tJ( i, 0 ) = ( tB( i ) - tA( i )) / ( 2.0 * tDeltaT );
                    tJ( i, 1 ) = ( tD( i ) - tC( i )) / ( 2.0 * tDeltaU );
                }

                tErr = std::sqrt( tX( 0 ) * tX( 0 ) + tX( 1 ) * tX( 1 ) );

                // catch crash
                if ( norm( tJ.col( 0 ) ) > 1e-12 )
                {
                    gesv( tJ, tX, tP );

                    // limit T
                    if( tThat < tX( 0 ) )
                    {
                        tOmegaT = tThat / tX( 0 ) * tOmega ;
                    }
                    else if ( tThat - tX( 0 ) > tTmax )
                    {
                        tOmegaT = ( tThat - tTmax ) / tX( 0 ) * tOmega ;
                    }
                    else if( abs( tX( 0 ) ) / tThat > 0.25 )
                    {
                        tOmegaT = 0.25 * tThat / abs( tX( 0 ) ) * tOmega ;
                    }
                    else
                    {
                        tOmegaT = tOmega ;
                    }

                    // limit U
                    if ( tUhat < tX( 1 ) )
                    {
                        tOmegaU = tUhat / tX( 1 ) * tOmega ;
                    }
                    else if( abs( tX( 1 ) ) / tUhat > 0.25 )
                    {
                        tOmegaU = 0.25 * tUhat / abs( tX( 1 ) ) * tOmega ;
                    }
                    else
                    {
                        tOmegaU = tOmega ;
                    }

                    tAlpha = tOmegaU < tOmegaT ? tOmegaU : tOmegaT ;

                    tThat -= tAlpha * tX( 0 );
                    tUhat -= tAlpha * tX( 1 );
                }
                else
                {
                    //std::cout << "warning: jacobi crash" << std::endl ;
                    //std::cout << "fail" << std::endl ;
                    //this->check_balance( tY  );
                    //tY.print("Balance");
                    break ;
                }

                // std::cout << " iterate " << tCount << " " << tThat << " " << tUhat << " " << tOmega << " " << tErr << std::endl ;

                BELFEM_ERROR( tCount++ < 500,"Too many iterations Tm=%12.3f K, um=%12.3f m/s, ,Tw=%12.3f K, T_hat=%12.3f K, u_hat=%12.3f m/s, Err=%12.3f",
                             ( double ) mTm , ( double  ) mUm , ( double  ) this->Tw() ,
                             ( double ) tThat , ( double  ) tUhat , ( double  ) tErr );

                if( tCount == 100 )
                {
                    tOmega = 0.1 ;
                    this->compute_initial_guesses() ;

                    tThat =  mData( mCenter, BELFEM_CHANNEL_T ) ;
                    tUhat =  mData( mCenter, BELFEM_CHANNEL_U ) ;
                }
            }

            this->check_balance( aBalance );

            return aBalance( 2 );
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::compute_inner_step(
                const   real   & aThat,
                const   real   & aUhat,
                const   real   & aPi,
                Vector< real > & aBalance  )
        {
            mPi = aPi ;
            this->set_center_conditions( aThat, aUhat );
            this->compute_wall_state() ;
            this->compute_velocity_profile();
            this->compute_temperature_profile() ;
            this->compute_turbulence();
            this->check_balance(aBalance );
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::compute_initial_guesses()
        {
            this->compute_wall_state() ;

            // Dittus-Boelther solution
            //mData( 0, BELFEM_CHANNEL_TAU )
            //    = 0.023 * std::pow( mReDhw, - 0.2 ) * this->rho_w() * mUm * mUm ;

            real tDotQ ;
            real tTauW ;
            real tHr ;
            real tTr ;
            this->friction_eckert( tDotQ, tTauW, tHr, tTr );

            mData( 0, BELFEM_CHANNEL_TAU ) = tTauW ;


            //mData( mCenter, BELFEM_CHANNEL_U ) = 1.1 * mUm ;

            this->compute_sigma_recovery_petrukov() ;

            /*real tThat = mTm ;
            while( this->check_u_hat_from_t_hat( tThat ) < 0.0 )
            {
                tThat += 0.05 * mTm ;
                std::cout << "try T_hat " << tThat << std::endl ;
            }*/

            //mData( mCenter, BELFEM_CHANNEL_T ) = tThat;

        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::compute_wall_state()
        {
            const real & tTw = this->Tw();

            // d'uh
            mData( 0, BELFEM_CHANNEL_Y ) = 0;
            mData( 0, BELFEM_CHANNEL_YPLUS ) = 0;
            mData( 0, BELFEM_CHANNEL_UPLUS ) = 0;
            mData( 0, BELFEM_CHANNEL_U ) = 0;

            // density at wall
            if ( mGas.is_idgas() )
            {
                mData( 0, BELFEM_CHANNEL_RHO ) = mGas.rho( tTw, mP );
            }
            else
            {
                mData( 0, BELFEM_CHANNEL_RHO ) = 1.0 / mVolumeSpline->eval( tTw );
            }

            // specific heat capacity at wall
            mData( 0, BELFEM_CHANNEL_CP )     = mHeatSpline->deval( tTw );

            // enthalpy at wall
            mData( 0, BELFEM_CHANNEL_H )      = mHeatSpline->eval( tTw );

            // viscosity at wall
            mData( 0, BELFEM_CHANNEL_MU )     = mMuSpline->eval( tTw );

            // thermal conductivity at wall
            mData( 0, BELFEM_CHANNEL_LAMBDA ) = mLambdaSpline->eval( tTw );
        }

//------------------------------------------------------------------------------


        void
        Boundarylayer::compute_velocity_profile()
        {
            // wall temperature
            real & T_w = mData( 0, BELFEM_CHANNEL_T );

            // wall density
            real & rho_w = mData( 0, BELFEM_CHANNEL_RHO );

            // shear stress at wall
            real & mu_w = mData( 0, BELFEM_CHANNEL_MU );

            // conductivity at wall
            real & lambda_w = mData( 0, BELFEM_CHANNEL_LAMBDA );

            // enthalpy at wall
            real & h_w = mData( 0, BELFEM_CHANNEL_H );

            // core temperature
            real & T_hat    = mData( mCenter, BELFEM_CHANNEL_T);

            // core density
            real & rho_hat = mData( mCenter, BELFEM_CHANNEL_RHO );

            // core velocity
            real & u_hat = mData( mCenter, BELFEM_CHANNEL_U );

            // core height
            real & y_hat = mData( mCenter, BELFEM_CHANNEL_Y );

            // core enthalpy
            real & h_hat = mData( mCenter, BELFEM_CHANNEL_H );

            // U+ in core
            //real & U_plus_hat = mData( mCenter, BELFEM_CHANNEL_UPLUS );

            // constants for Crocco-Buseman profile
            // ( 21 )
            mPsi = ( mu_w / ( mSigma * lambda_w ) )
                   * ( h_hat + mRecovery * 0.5 * u_hat * u_hat - h_w )
                   * mGas.alpha( T_w, mP );

            // ( 20 )
            mPhi = 1.0 + mPsi - rho_w / rho_hat;


            const cplx i( 0.0, 1.0 );
            const cplx tSqPhi = std::sqrt( mPhi ) ;
            cplx  tTheta ;

            //BELFEM_ERROR( std::abs( tSqPhi ) > 1e-12,
            //             "Correlation Fail for T=%10.3f K, p=%10.3f bar, u=%10.3f m/s, Tw=%10.3f",
            //             ( double ) T_hat, ( double ) mP , ( double ) u_hat, ( double ) T_w );

            mChi = std::sqrt(  mPsi * mPsi + 4.0 * mPhi );

            if( ! ( std::abs( mChi ) > 1e-12 || std::abs( tSqPhi ) > 1e-12 ) )
            {
                std::cout << "sigma " << mSigma << std::endl ;
                std::cout << "recovery " << mRecovery << std::endl ;
                std::cout << "alpha " << mGas.alpha( T_w, mP ) << std::endl ;
                std::cout << "phi " << mPhi << std::endl ;
                std::cout << "psi " << mPsi << std::endl ;
                std::cout << "h " << ( h_hat + mRecovery * 0.5 * u_hat * u_hat - h_w )<< std::endl ;
                std::cout << "lambda " << lambda_w << std::endl ;
            }

            // Check that this is not equal to zero
            BELFEM_ERROR( std::abs( mChi ) > 1e-12 || std::abs( tSqPhi ) > 1e-12 ,
                         "Correlation Fail for Dh=%10.3f mm, T=%10.3f K, p=%10.3f bar, u=%10.3f m/s, Tw=%10.3f",
                         ( double ) mDh, ( double ) T_hat, ( double ) mP , ( double ) u_hat, ( double ) T_w );

            const cplx tF0 = std::log( 2.0 + i / tSqPhi ) * i / tSqPhi ;
            cplx tChi ;
            real tUpsilon ;

            // initual guess for F+
            real tFplus = 0;
            real tdUplusdU;
            real tdFdYplus;
            real tdGdYplus;

            // backup center velocity
            real tUhat = u_hat ;

            mBeta = std::asin( mPsi / mChi );

            this->compute_shear_stress();

            mAlpha = std::sqrt( mPhi ) * mUtau / u_hat;

            for( index_t k=0; k<mNumberOfNodes; ++k )
            {
                // use last value as initial guess and compute the wall function
                tFplus = boundarylayer::spalding(
                        mBplus, mKarman, mExpKB, mData( k,BELFEM_CHANNEL_YPLUS ), tFplus );

                // compute velocity function
                mData( k, BELFEM_CHANNEL_UPLUS ) = tFplus +  boundarylayer::g_plus( mKarman, mPi, mEta( k ) );

                // catch special case
                //if( std::abs( mChi ) < 1e-12 )
                if( false )
                {
                    tTheta = tSqPhi * ( mData( k, BELFEM_CHANNEL_UPLUS ) * mUtau + tF0 ) * i;
                    mData( k, BELFEM_CHANNEL_U ) = std::real(
                            -(( std::exp( -2.0 * tTheta ) * std::pow( std::exp( tTheta ) + i * tSqPhi, 2 )) /
                              ( 4.0 * mPhi ) + 1.0 )
                            / ( mPsi - std::exp( -tTheta ) * ( std::exp( tTheta ) + i * tSqPhi ))) * u_hat;

                    tUpsilon = mData( k, BELFEM_CHANNEL_U ) / tUhat;

                    tChi = std::sqrt( 1.0 + tUpsilon * ( mPsi - mPhi * tUpsilon ));

                    tdUplusdU = std::real( i * (( mPsi - 2.0 * mPhi * tUpsilon ) / tChi - 2.0 * i * tSqPhi ) /
                                           ( tSqPhi * mUtau *
                                             ( 2.0 * tChi - i * ( 2.0 * mPhi * tUpsilon - 1.0 ) / tSqPhi )));

                }
                else
                {
                    // transform to real value
                    mData( k, BELFEM_CHANNEL_U ) =std::real(
                            ( std::sin( mData( k, BELFEM_CHANNEL_UPLUS ) * mAlpha - mBeta ) * mChi + mPsi )
                            / ( 2.0 * mPhi )) * u_hat;

                    // BELFEM_ASSERT( mData( k, BELFEM_CHANNEL_U ) > -1e-6, "Transformation Fail: k=%u, u=%12.3f", ( unsigned int ) k , ( double   )mData( k, BELFEM_CHANNEL_U )  );

                    // derivative
                    tdUplusdU = std::real( 2.0 * mPhi / ( mAlpha * mChi * u_hat *
                                                          std::cos( mBeta - mAlpha * mData( k, BELFEM_CHANNEL_UPLUS ))));
                }
                tdFdYplus = 1.0 / boundarylayer::spalding_dydf( mBplus, mKarman, mExpKB, tFplus );

                // dgdeta * detady * dydy
                tdGdYplus =  boundarylayer::dg_plus_deta( mKarman, mPi, mEta( k ) ) / ( y_hat * mCplus );

                // derivative of u
                mData( k, BELFEM_CHANNEL_DUDY ) = ( tdFdYplus + tdGdYplus ) / tdUplusdU * mCplus ;
            }

            // write center velocity back ( to avoid numerical dust )
            mData( 0, BELFEM_CHANNEL_U ) = 0.0 ;
            mData( mCenter, BELFEM_CHANNEL_U ) = tUhat ;
            mData( mCenter, BELFEM_CHANNEL_DUDY ) = 0.0 ;
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::compute_temperature_profile()
        {
            // compute densities
            real tUpsilon ;
            real tUhat = mData( mCenter, BELFEM_CHANNEL_U ) ;
            real tRhow = mData( 0, BELFEM_CHANNEL_RHO ) ;

            // save center density
            real tRho_hat = mData( mCenter, BELFEM_CHANNEL_RHO );

            real tF ;
            real tdF ;
            uint tCount ;

            // density profile
            for( uint k=1; k<mNumberOfNodes; ++k )
            {
                real & tT = mData( k, BELFEM_CHANNEL_T );
                real & tRho = mData( k, BELFEM_CHANNEL_RHO );

                tT = mData( k-1, BELFEM_CHANNEL_T );

                tUpsilon = mData( k, BELFEM_CHANNEL_U ) / tUhat ;

                tRho = std::min( tRhow / std::real( 1.0 + tUpsilon * ( mPsi - tUpsilon * mPhi ) ), mRhoMax );

                tF = BELFEM_REAL_MAX ;

                tCount = 0 ;
                if( mGas.is_idgas() )
                {
                    tT = mP / ( mGas.R( tT, mP ) * tRho );
                }
                else
                {
                    while ( abs( tF ) > 1e-9 )
                    {
                        tF = mVolumeSpline->eval( tT ) - 1.0 / tRho;
                        tdF = mVolumeSpline->deval( tT );

                        tT -= 0.9 * tF / tdF;
                        BELFEM_ERROR( tCount++ < 100, "Too many iterations" );
                    }
                }

                mData( k, BELFEM_CHANNEL_CP )     = mHeatSpline->deval( tT );
                mData( k, BELFEM_CHANNEL_H )      = mHeatSpline->eval( tT );
                mData( k, BELFEM_CHANNEL_MU )     = mMuSpline->eval( tT );
                mData( k, BELFEM_CHANNEL_LAMBDA ) = mLambdaSpline->eval( tT );
            }
            // backup to avoid numerical noise
            mData( mCenter, BELFEM_CHANNEL_RHO ) = tRho_hat ;
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::compute_turbulence()
        {

            for( uint k=0; k<mNumberOfNodes; ++k )
            {
                // laminar viscosity
                real & cp = mData( k, BELFEM_CHANNEL_CP );

                // laminar viscosity
                real & mu = mData( k, BELFEM_CHANNEL_MU );

                // laminar thermal conductivity
                real & lambda = mData( k, BELFEM_CHANNEL_LAMBDA );

                // prandtl number
                real & Pr = mData( k, BELFEM_CHANNEL_PR );

                // turbulent viscosity
                real & mu_T = mData( k, BELFEM_CHANNEL_MUT );

                // turbulent thermal conductivity
                real & lambda_T = mData( k, BELFEM_CHANNEL_LAMBDAT );

                // turbulent Prandtl number
                real & Pr_T = mData( k, BELFEM_CHANNEL_PRT );

                // mixed Prandtl number
                real & Pr_M = mData( k, BELFEM_CHANNEL_PRM );

                // compute the prandtl number
                Pr = mu * cp / lambda ;

                // compute the turbulent viscosity
                mu_T = mu * mKarman * mData( k, BELFEM_CHANNEL_YPLUS );

                // compute the turbulent Prandtl number
                Pr_T = boundarylayer::kays_crawford( Pr, mu, mu_T ) ;

                // turbulent thermal conductivity
                lambda_T = mu_T * cp / Pr_T ;

                // mixed prandtl number
                Pr_M = ( mu + mu_T ) * cp / ( lambda + lambda_T );

                // shear stress
                mData( k, BELFEM_CHANNEL_TAU ) = ( mu + mu_T ) * std::abs( mData( k, BELFEM_CHANNEL_DUDY ) );
            }

            this->derive( BELFEM_CHANNEL_TAU, BELFEM_CHANNEL_DTAUDY );
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::compute_sigma_recovery_vandriest()
        {

            real tRdh = 1.0 ; //mIsAxisymmetric ? 0.5 * mDh : 1.0 ;

            for( index_t k=0; k<mCenter; ++k )
            {
                mData( k, BELFEM_CHANNEL_WORK0 ) =
                        ( 1.0 - mData( k, BELFEM_CHANNEL_PRM ) ) * mData( k, BELFEM_CHANNEL_DTAUDY ) / ( tRdh * mData( k, BELFEM_CHANNEL_TAU ) );
            }

            // fix last value, it would be infinity since tau is zero in the center
            mData( mCenter, BELFEM_CHANNEL_WORK0 )
                =  mData( mCenter-2, BELFEM_CHANNEL_WORK0 )
                 + ( mData( mCenter, BELFEM_CHANNEL_Y )-mData( mCenter-2, BELFEM_CHANNEL_Y ) ) /
                   ( mData( mCenter-1, BELFEM_CHANNEL_Y )-mData( mCenter-2, BELFEM_CHANNEL_Y ) )
                 * ( mData( mCenter-1, BELFEM_CHANNEL_WORK0 ) - mData( mCenter-2, BELFEM_CHANNEL_WORK0 ) );

            // compute sigma
            this->integrate( BELFEM_CHANNEL_WORK0 , BELFEM_CHANNEL_WORK1, false );


            for( index_t k=0; k<mNumberOfNodes; ++k )
            {
                mData( k, BELFEM_CHANNEL_WORK0 )
                    =   mData( k, BELFEM_CHANNEL_PRM )
                      * std::exp( -mData( k, BELFEM_CHANNEL_WORK1 ) )
                      * mData( k, BELFEM_CHANNEL_DUDY ) / tRdh ;
            }
            this->integrate( BELFEM_CHANNEL_WORK0 , BELFEM_CHANNEL_WORK2, false );

            mSigma = mData( mCenter, BELFEM_CHANNEL_WORK2 ) / mData( mCenter, BELFEM_CHANNEL_U );

            // help function for recovery factor
            for( index_t k=0; k<mNumberOfNodes; ++k )
            {
                mData( k, BELFEM_CHANNEL_WORK2 )
                        =  std::exp( mData( k, BELFEM_CHANNEL_WORK1 ) ) * mData( k, BELFEM_CHANNEL_DUDY ) / tRdh ;
            }
            this->integrate( BELFEM_CHANNEL_WORK2 , BELFEM_CHANNEL_WORK1, false );

            // integral for recovery factor
            for( index_t k=0; k<mNumberOfNodes; ++k )
            {
                mData( k, BELFEM_CHANNEL_WORK2 ) =
                        mData( k, BELFEM_CHANNEL_WORK0 ) * mData( k, BELFEM_CHANNEL_WORK1 ) ;
            }

            this->integrate( BELFEM_CHANNEL_WORK2 , BELFEM_CHANNEL_WORK0, false );

            mRecovery = 2.0 * mData( mCenter, BELFEM_CHANNEL_WORK0 ) /
                   (  mData( mCenter, BELFEM_CHANNEL_U ) * mData( mCenter, BELFEM_CHANNEL_U ) );

            if( ( std::abs( mSigma ) > 10 ) ||  ( std::abs( mRecovery ) > 10 ) )
            {
                this->print();
            }

            /*BELFEM_ERROR( ( std::abs( mSigma ) < 10 ) &&  ( std::abs( mRecovery ) < 10 ),
                         "Invalid recovery values for averaged values of Dh=%10.6f mm, T=%10.3f K, p=%10.3f bar, u=%10.3f m/s, Tw=%10.3f, Y+=%12.4f",
                         ( double ) mDh * 1000,
                         ( double ) mTm,
                         ( double ) mP *1e-5 ,
                         ( double ) mUm,
                         ( double ) this->Tw(),
                         ( double ) mData( 1, BELFEM_CHANNEL_YPLUS )
                         ); */

            if( ( std::abs( mSigma ) > 20 ) ||  ( std::abs( mRecovery ) > 20 ) )
            {
                std::fprintf( stdout , "    Warning: Invalid recovery values for averaged values of Dh=%10.6f mm, \n    T=%10.3f K, p=%10.3f bar, u=%10.3f m/s, Tw=%10.3f, Y+=%12.4f, \n    sigma = %12.4f, rec = %12.4f",
                              ( double ) mDh * 1000,
                              ( double ) mTm,
                              ( double ) mP *1e-5 ,
                              ( double ) mUm,
                              ( double ) this->Tw(),
                              ( double ) mData( 1, BELFEM_CHANNEL_YPLUS ),
                              ( double ) mSigma,
                              ( double ) mRecovery );

                mRecovery = std::pow( mData( 0, BELFEM_CHANNEL_PR ), 1.0/3.0 );
                mSigma = mRecovery * mRecovery ;

            }

            }
//------------------------------------------------------------------------------

        void
        Boundarylayer::compute_sigma_recovery_petrukov()
        {
            // friction factor with respect to center conditions
            real tCf =  2.0 * this->tau_w() / ( mRhom * mUm * mUm );

            mRecovery = std::pow( mPrm, 1.0/3.0 );

            // Petrukov equation, see 10.1016/S0065-2717(08)70153-9
            // or my thesis Eq. ( 2. 97 ), using reference conditions
            mSigma = ( 1.0 + 13.6 * tCf ) + ( 11.7  + 1.8 / mRecovery )
                 * ( mRecovery * mRecovery - 1.0 ) * std::sqrt( 0.5 * tCf );
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::check_balance( Vector< real > & aBalance )
        {
            mIsAxisymmetric = true ;

            real tScale = 2.0 * constant::pi ;

            // integrate balance equations for mass and momentum
            for( index_t k=0; k<mNumberOfNodes; ++k )
            {
                // Mass equation
                mData( k, BELFEM_CHANNEL_WORK0 )
                    = mData( k, BELFEM_CHANNEL_RHO ) * mData( k, BELFEM_CHANNEL_U ) * tScale ;

                // momentum equation
                mData( k, BELFEM_CHANNEL_WORK1 )
                        = mData( k, BELFEM_CHANNEL_WORK0 ) * mData( k, BELFEM_CHANNEL_U ) ;
            }

            this->integrate( BELFEM_CHANNEL_WORK0, BELFEM_CHANNEL_WORK2, mIsAxisymmetric );
            aBalance( 0 ) = ( mData( mCenter, BELFEM_CHANNEL_WORK2 ) - mDotM ) / mDotM ;

            this->integrate( BELFEM_CHANNEL_WORK1, BELFEM_CHANNEL_WORK2, mIsAxisymmetric );

            aBalance( 1 ) = ( mData( mCenter, BELFEM_CHANNEL_WORK2 ) - mDotI ) / mDotI ;

            if( aBalance.length() > 2 )
            {
                // integrate balance equations for energy
                for ( index_t k = 0; k < mNumberOfNodes; ++k )
                {
                    mData( k, BELFEM_CHANNEL_WORK1 )
                            = ( mData( k, BELFEM_CHANNEL_H ) +
                                0.5 * mData( k, BELFEM_CHANNEL_U ) * mData( k, BELFEM_CHANNEL_U ) )
                              * mData( k, BELFEM_CHANNEL_WORK0 )  ;

                }
                this->integrate( BELFEM_CHANNEL_WORK1, BELFEM_CHANNEL_WORK2, mIsAxisymmetric );
                real tH = mData( mCenter, BELFEM_CHANNEL_WORK2 ) / mDotM
                        - 0.5 * mData( mCenter, BELFEM_CHANNEL_U ) * mData( mCenter, BELFEM_CHANNEL_U ) ;

                real tT = mGas.T_from_h( tH, mP ) ;

                aBalance( 2 ) = ( tT - mTm ) / mTm;
            }
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::derive( const uint aSource, const uint aTarget )
        {
            index_t tOff = 0;

            mData( 0, aTarget ) = 0.0 ;

            real tLength ;

            // loop over all cells
            for( index_t e=0; e<mNumberOfCells; ++e )
            {
                // get cell length
                tLength = mData( tOff + 2 , BELFEM_CHANNEL_Y ) - mData( tOff, BELFEM_CHANNEL_Y );

                // get values
                real & tVal0 = mData( tOff, aSource );
                real & tVal1 = mData( tOff+2, aSource );
                real & tVal2 = mData( tOff+1, aSource );

                // compute derivatives
                mData( tOff,   aTarget ) += ( 2.0 * tVal2 - 1.5 * tVal0 - 0.5 * tVal1 ) / tLength ;
                mData( tOff+1, aTarget ) = ( tVal1 - tVal0 ) / tLength ;
                mData( tOff+2, aTarget ) = ( 0.5 * tVal0 + 1.5 * tVal1 - 2.0 * tVal2 ) / tLength ;

                tOff += 2 ;
            }

            mData( 0, aTarget ) *= 2 ;
            mData( mCenter, aTarget ) *= 2 ;
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::compute_shear_stress()
        {
            // relaxation factor
            real tOmega = 0.9 ;
            real tAlpha ;
            real dx ;

            // tauw
            real & tau_w = mData( 0, BELFEM_CHANNEL_TAU );

            // density at wall
            const real & rho_w = this->rho_w();

            // shear stress at wall
            const real & mu_w = this->mu_w() ;

            // velocity at center of flow
            const real & u_hat = mData( mCenter, BELFEM_CHANNEL_U );

            const real y = 0.5 * mDh ;

            /**
             * The original function crashes since the asine function is not bijective!
             * apparently, some information gets lost when inverting the sin function.
             *
             * The solution is to use the log functions instead. Blimey!
             */
            //real A = std::real( ( std::asin(
            //        (2.0 * mPhi - mPsi ) / mChi ) + mBeta ) / std::sqrt( mPhi ) );

            const cplx i( 0.0, 1.0 );
            cplx a = std::sqrt( mPhi );

            cplx cA0 = i * std::log( ( 2.0 + i * mPsi / a ) ) ;
            cplx cA1 = i * std::log( 2.0 * std::sqrt( 1.0 + mPsi - mPhi ) + i * ( mPsi - 2.0 * mPhi ) / a );
            real A = std::real( ( cA1 - cA0 ) / a );

            real & u_tau = mUtau ;

            real U ;
            real Y ;
            real K ;
            real L ;
            real LB ;
            real & B = mBplus ;
            real G = boundarylayer::g_plus( mKarman, mPi, 1.0 );
            real F = BELFEM_REAL_MAX ;

            real dY;
            real dU ;
            real dL ;
            real dK ;
            real dB ;
            real dF ;

            index_t tCount = 0 ;

            // compute shear stress
            u_tau = std::sqrt( tau_w / rho_w );

            real x = 1.0 / u_tau ;

            dU = A * u_hat ;
            L = 0.0 ;

            while ( abs( F ) >1e-12 )
            {
                // compute the shear stress
                u_tau = 1.0 / x ;

                // U+
                U = dU * x ;

                // Y+
                Y = rho_w * u_tau * y  / mu_w ;
                dY = -Y / x ;

                // k-plus term for roughness
                K = rho_w * mKtech * u_tau / mu_w;
                dK = -K / x ;

                // B-Term offset for loglaw
                B = mBplus0 - std::log( 1.0 + K / mRoughConst ) / mKarman;
                dB = -dK / ( mKarman * ( K + mRoughConst ));


                // Log-Term for loglaw
                // L = std::log( Y ) / mKarman ;
                mExpKB = std::exp( -mKarman * mBplus );
                LB = boundarylayer::spalding( B, mKarman, mExpKB, Y, L );
                dL = dY / ( Y * mKarman );

                // the full function
                F = LB + G - U;
                dF = dL + dB - dU;

                dx = F/dF ;

                if( x < dx )
                {
                    tAlpha = x / dx * tOmega ;
                }
                else
                {
                    tAlpha = tOmega ;
                }

                x -= tAlpha * dx ;


                if( tCount >= 100 )
                {
                    std::cout << "tau_w fail" << std::endl ;
                    this->print();
                }
                BELFEM_ERROR( tCount++ < 100, "Too many iterations , tau_w=%8.3f", u_tau * u_tau * rho_w );
            }

            tau_w =  u_tau * u_tau * rho_w ;

            mCplus = this->rho_w() * mUtau / this->mu_w() ;

            // compute wall distance
            for( index_t k=0; k<mNumberOfNodes; ++k )
            {
                mData( k, BELFEM_CHANNEL_YPLUS ) = mCplus * mData( k, BELFEM_CHANNEL_Y ) ;
            }

        }

//------------------------------------------------------------------------------
        void
        Boundarylayer::integrate( const uint aSource, const uint aTarget, const bool aIsAxisymmetric )
        {
            index_t tOff = 0;

            mData( 0, aTarget ) = 0.0;

            real tLength;

            // loop over all cells
            if ( aIsAxisymmetric )
            {
                for ( index_t e = 0; e < mNumberOfCells; ++e )
                {
                    // get cell length
                    tLength = mData( tOff + 2, BELFEM_CHANNEL_Y ) - mData( tOff, BELFEM_CHANNEL_Y );

                    // get nodal values
                    real tVal0 = mData( tOff, aSource ) * mData( tOff, BELFEM_CHANNEL_R );
                    real tVal1 = mData( tOff + 2, aSource ) * mData( tOff + 2, BELFEM_CHANNEL_R );
                    real tVal2 = mData( tOff + 1, aSource ) * mData( tOff + 1, BELFEM_CHANNEL_R );

                    // integrate up to the middle point
                    mData( tOff + 1, aTarget ) =
                            mData( tOff, aTarget ) + ( 5.0 * tVal0 - tVal1 + 8.0 * tVal2 ) * tLength / 24.0;

                    // integrate along the whole cell using the simpson rule
                    mData( tOff + 2, aTarget ) =
                            mData( tOff, aTarget ) + ( tVal0 + tVal1 + 4.0 * tVal2 ) * tLength / 6.0;
                    tOff += 2;
                }
            }
            else
            {
                for ( index_t e = 0; e < mNumberOfCells; ++e )
                {
                    // get cell length
                    tLength = mData( tOff + 2, BELFEM_CHANNEL_Y ) - mData( tOff, BELFEM_CHANNEL_Y );

                    // get nodal values
                    real & tVal0 = mData( tOff, aSource );
                    real & tVal1 = mData( tOff + 2, aSource );
                    real & tVal2 = mData( tOff + 1, aSource );

                    // integrate up to the middle point
                    mData( tOff + 1, aTarget ) =
                            mData( tOff, aTarget ) + ( 5.0 * tVal0 - tVal1 + 8.0 * tVal2 ) * tLength / 24.0;

                    // integrate along the whole cell using the simpson rule
                    mData( tOff + 2, aTarget ) =
                            mData( tOff, aTarget ) + ( tVal0 + tVal1 + 4.0 * tVal2 ) * tLength / 6.0;

                    tOff += 2;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::init_lookup_tables()
        {
            real tTmin ;
            real tTmax ;

            // get minumum temperature
            if ( mGas.gas_model() == GasModel::HELMHOLTZ )
            {
                gasmodels::Helmholtz * tEoS
                        = reinterpret_cast< gasmodels::Helmholtz * > ( mGas.eos() ) ;
                mTmin = tEoS->T_min() ;
                tTmin = 0.75 * tEoS->T_min() ;
                mTmax = 2000 ; //1.25 * tEoS->T_max() ;
                tTmax = 2200 ; //1.5 * tEoS->T_max() ;
                mNumSplineSteps = 101 ;
            }
            else
            {
                mTmin = 200.0  ;
                mTmax = 6000.0   ;
                tTmin = 150.0  ;
                tTmax = 6200.0 ;
                mNumSplineSteps = 201 ;
            }
            real tDeltaT = ( tTmax - tTmin ) / mNumSplineSteps ;

            spline::create_helpmatrix( mNumSplineSteps, tDeltaT, mHelpMatrix );

            linspace( tTmin, tTmax, mNumSplineSteps, mWorkTemperature );

            mWorkVolume.set_size( mNumSplineSteps );
            mWorkHeat.set_size( mNumSplineSteps );
            mWorkMu.set_size( mNumSplineSteps );
            mWorkLambda.set_size( mNumSplineSteps );

            mVolumeSpline = new Spline( mNumSplineSteps, tTmin, tTmax );
            mHeatSpline   = new Spline( mNumSplineSteps, tTmin, tTmax );
            mMuSpline     = new Spline( mNumSplineSteps, tTmin, tTmax );
            mLambdaSpline = new Spline( mNumSplineSteps, tTmin, tTmax );
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::update_lookup_tables()
        {
            if( mGas.gas_model() == GasModel::HELMHOLTZ )
            {
                gasmodels::Helmholtz * tEoS
                        = reinterpret_cast< gasmodels::Helmholtz * > ( mGas.eos() ) ;

                // maximum temperture where original model can be used
                real tT0     = tEoS->T_max() - 50.0 ;

                // temperature where alternate model begins
                real tT1    = 1000.0 ;

                // step for derivative
                real tDeltaT = 0.001 ;

                // create volume polynomial
                real tY0     = mGas.v( tT0, mP );
                real tdYdT0  = tY0 * mGas.alpha( tT0, mP );
                real tY1     = mAlternateGas->v( tT1, mP );
                real tdYdT1  = tY1 * mAlternateGas->alpha( tT1, mP );

                Vector< real > tVolumePoly;
                create_beam_poly( tT0, tY0, tdYdT0, tT1, tY1, tdYdT1, tVolumePoly );

                // get the reference gas
                gastables::RefGas * tRefgas = mAlternateGas->component( 0 ) ;

                // create enthalpy polynomial

                // enthalpy offset

                tY0     = mGas.h( tT0, mP );
                tdYdT0  = mGas.cp( tT0, mP );

                real tDeltaH = tY0 - mAlternateGas->h( tT0, mP );

                tY1     = mAlternateGas->h( tT1, mP ) + tDeltaH ;
                tdYdT1  = mAlternateGas->cp( tT1, mP );

                Vector< real > tHeatPoly;
                create_beam_poly( tT0, tY0, tdYdT0, tT1, tY1, tdYdT1, tHeatPoly );


                // create viscosity polynomial
                tY0    = mGas.mu( tT0, mP );
                real tDeltaMu = tY0 - tRefgas->mu( tT0 );
                tY1     = tRefgas->mu( tT1 ) + tDeltaMu ;

                // derivatives
                tdYdT0 = (
                          mGas.mu( tT0+tDeltaT, mP )
                        - mGas.mu( tT0 - tDeltaT, mP ) ) / ( 2.0 * tDeltaT );


                tdYdT1 = tRefgas->dmudT( tT1 );

                Vector< real > tViscosityPoly;
                create_beam_poly( tT0, tY0, tdYdT0, tT1, tY1, tdYdT1, tViscosityPoly );

                // create conductivity polynomial
                tY0    = mGas.lambda( tT0, mP );
                real tDeltaLambda = tY0 - tRefgas->lambda( tT0 );
                tY1     = tRefgas->lambda( tT1 ) + tDeltaLambda ;

                // derivatives
                tdYdT0 = (  mGas.lambda( tT0+tDeltaT, mP )
                          - mGas.lambda( tT0 - tDeltaT, mP ) ) / ( 2.0 * tDeltaT );

                tdYdT1 = tRefgas->dlambdadT( tT1 );

                Vector< real > tConductivityPoly;
                create_beam_poly( tT0, tY0, tdYdT0, tT1, tY1, tdYdT1, tConductivityPoly );

                uint k = 0 ;

                real tT = mWorkTemperature( 0 ) ;

                // cold part
                while( tT < tT0 )
                {
                    mWorkVolume( k ) = mGas.v( tT, mP );
                    mWorkHeat( k ) = mGas.h( tT, mP );
                    mWorkMu( k ) = mGas.mu( tT, mP );
                    mWorkLambda( k ) = mGas.lambda( tT, mP );

                    ++k ;
                    tT = mWorkTemperature( k );
                }

                // transition
                while( tT < tT1 )
                {
                    mWorkVolume( k ) = polyval( tVolumePoly, tT );
                    mWorkHeat( k )   = polyval( tHeatPoly, tT ) ;
                    mWorkMu( k )     = polyval( tViscosityPoly, tT );
                    mWorkLambda( k ) = polyval( tConductivityPoly, tT );

                    ++k ;
                    tT = mWorkTemperature( k );
                }

                // alternate model
                while( k < mNumSplineSteps )
                {
                    tT = mWorkTemperature( k );

                    mWorkVolume( k ) = mAlternateGas->v( tT, mP );
                    mWorkHeat( k )   = mAlternateGas->h( tT, mP ) + tDeltaH ;
                    mWorkMu( k )     = tRefgas->mu( tT ) + tDeltaMu ;
                    mWorkLambda( k ) = tRefgas->lambda( tT ) + tDeltaLambda ;

                    ++k ;
                }
            }
            else
            {
                // populate Properties
                for ( uint k = 0; k < mNumSplineSteps; ++k )
                {
                    const real & tT = mWorkTemperature( k );

                    mWorkVolume( k ) = mGas.v( tT, mP );
                    mWorkHeat( k ) = mGas.h( tT, mP );
                    mWorkMu( k ) = mGas.mu( tT, mP );
                    mWorkLambda( k ) = mGas.lambda( tT, mP );
                }
            }

            mVolumeSpline->update_data( mHelpMatrix, mWorkVolume );

            mHeatSpline->update_data( mHelpMatrix, mWorkHeat );

            mMuSpline->update_data( mHelpMatrix, mWorkMu );
            mLambdaSpline->update_data( mHelpMatrix, mWorkLambda );

            mRhoMax = mGas.rho( mTmin, mP );
        }

//------------------------------------------------------------------------------

        real
        Boundarylayer::compute_heatflux()
        {
            const real & rho_w = mData( 0, BELFEM_CHANNEL_RHO );
            const real & mu_w  = mData( 0, BELFEM_CHANNEL_MU );
            const real & lambda_w  = mData( 0, BELFEM_CHANNEL_LAMBDA );

            real dudy = this->tau_w() / mu_w ;

            // compute gradient of drho_dy
            real drhody = - std::real( mPsi ) * rho_w * dudy / mData( mCenter, BELFEM_CHANNEL_U ) ;

            real dvdT = mGas.alpha( this->Tw(), mP ) / rho_w ;
            real drhodT = - rho_w * rho_w * dvdT ;

            real dTdy = drhody / drhodT ;

            return lambda_w * dTdy ;
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::friction_messe( real & aDotQ, real & aTauw, real & aHr, real & aTr )
        {
            // Balances
            mDotM = mRhom * mUm * mA ;
            mDotI = ( mRhom * mUm * mUm ) * mA ;
            mDotH = ( mHm + 0.5 * mUm * mUm ) * mDotM ;

            this->compute_outer_step( mPi, mBalance );

            const real & u_hat = mData( mCenter, BELFEM_CHANNEL_U ) ;
            const real & h_hat = mData( mCenter, BELFEM_CHANNEL_H );

            // compute recovery enthalpy
            aHr = h_hat + 0.5 * mRecovery * u_hat * u_hat;

            aTauw = mData( 0, BELFEM_CHANNEL_TAU );

            aDotQ = aTauw * ( aHr - mData( 0, BELFEM_CHANNEL_H ) ) / ( mSigma * u_hat );

            aTr = mGas.T_from_h( aHr, mP );

        }

//------------------------------------------------------------------------------
        void
        Boundarylayer::friction_bartz( real & aDotQ, real & aTauw, real & aHr, real & aTr )
        {
            // make sure that this is a combustion gas
            BELFEM_ASSERT( mGas.number_of_components() > 1 && mGas.is_idgas(),
                          "The Bartz correlation can only be used with a combustion gas" );

            BELFEM_ASSERT( ! std::isnan( mBartzConst ),
                          "The Bartz specific geometry parameters have not been set" );

            // get the wall temperature
            const real & Tw = this->Tw() ;

            // specific heat
            real cp = mGas.cp( mTm, mP );

            // adiabatic coefficient
            real k = mGas.gamma( mTm, mP );

            // Mach number
            real Ma = mUm / mGas.c( mTm, mP );

            // viscosity
            real mu = mGas.mu( mTm, mP );

            // Prandtl Number
            real Pr = mGas.Pr( mTm, mP );

            // update the average Reynolds number
            mReDh = mRhom * mUm * mDh / mu ;

            // power law constant
            // that yields mu / mu_ref = ( T / T_ref )^omega
            real omega = std::log( mu / mGas.mu( Tw, mP ) ) / std::log( mTm / Tw ) ;

            // correction factor
            real sigma = std::pow( 0.5 * Tw / mTm * ( 1.0 + 0.5 * ( k-1. ) * Ma * Ma ) + 0.5, 0.2 * omega - 0.8 )
                    * std::pow( 1.0 + 0.5 * ( k - 1. ) * Ma * Ma, -0.2 * omega );

            // Reynolds-colburn analogy
            mSigma = std::pow( Pr, 0.6 ) ;

            // heat transfer coefficient this is really the Bartz equation,
            // we just don't refer to the throat conditions. This is OK!
            real alpha = 0.026 * cp / mSigma * std::pow( mReDh, -0.2 ) * mBartzConst
                    * mRhom * mUm * sigma ;

            // compute heat flux
            aDotQ = alpha * ( mTm - Tw );

            // recovery enthalpy
            mRecovery = 1.0 ;
            aHr = mHm + mRecovery * 0.5 * mUm * mUm ;

            // recovery temperature
            aTr = mGas.T_from_h( aHr, mP );

            // stanton number
            real St = aDotQ / ( mRhom * mUm * ( aHr - mGas.h( Tw, mP ) ) );

            // friction factor
            real cf = 2.0 * St * mSigma ;

            // shear stress
            aTauw = 0.5 * cf * mRhom * mUm * mUm ;
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::friction_eckert( real & aDotQ, real & aTauw, real & aHr, real & aTr )
        {
            // wall temperature
            const real & Tw = this->Tw() ;

            // density
            real rho = mGas.rho( mTm, mP ) ;

            // wall enthalpy
            real hw = mGas.h( Tw, mP );

            // reference temperature
            real T_ref = boundarylayer::reference_temperature( mGas, mTm, mP, mUm, Tw, true );

            // reference density
            real rho_ref = mGas.rho( T_ref, mP );

            // recovery factor
            mRecovery = std::pow( mGas.Pr( T_ref, mP ), 1./3. );

            // Reynolds-Colburn
            mSigma = mRecovery * mRecovery ;

            // Reynolds number
            real Re_Dh = rho_ref * mUm * mDh / mGas.mu( T_ref, mP );

            // read friction from moody chart
            real cf = boundarylayer::cf_moody( Re_Dh, mDh, mKtech ) * rho_ref / rho ;

            // shear stress
            aTauw = 0.5 * cf * rho * mUm * mUm ;

            // recovery enthalpy
            aHr = mHm + 0.5 * mRecovery * mUm * mUm ;

            // recovery temperature
            aTr = mGas.T_from_h( aHr, mP );

            // heat load
            aDotQ = aTauw / ( mSigma * mUm ) * ( aHr - hw );
        }

//------------------------------------------------------------------------------

        void
        Boundarylayer::friction_pizzarelli( real & aDotQ, real & aTauw, real & aHr, real & aTr )
        {
            BELFEM_ASSERT( mGas.component( 0 )->label() == "CH4" && mGas.number_of_components() == 1,
                "the Pizzarelli correlation can only be used with methane" );

            // properties at bulk temperature
            const real & Tb = mTm ;
            const real & p  = mP ;

            // critical pressure
            const real & pc = mGas.component( 0 )->data()->p_crit() ;

            // density
            real rho_b = mGas.rho( Tb, p );

            // specific heat
            real cp_b = mGas.cp( Tb, p ) ;

            // viscosity
            real mu_b = mGas.mu( Tb, p );

            // thermal conductivity
            real k_b = mGas.lambda( Tb, p );


            // prandtl number
            real Pr = mGas.Pr( Tb, p );

            // properties at wall
            const real & Tw = this->Tw() ;

            // density at wall
            real rho_w = mGas.rho( Tw, p );

            // enthalpy at wall
            real hw = mGas.h( Tw, p );

            // viscosity at wall
            real mu_w = mGas.mu( Tw, p );

            // thermal conductivit at wall
            real k_w = mGas.lambda( Tw, p );

            // mean specific heat
            real cp_m = ( hw - mHm ) / ( Tw - Tb );

            // Nusselt number, Eq. (29) in paper
            real Nu = 0.0272
                    * std::pow( mReDh,         0.8 )
                    * std::pow( Pr,            0.353 )
                    * std::pow( Tw / Tb,      -0.607 )
                    * std::pow( rho_w / rho_b, 0.357 )
                    * std::pow( mu_w / mu_b,  -0.662 )
                    * std::pow( k_w / k_b,     0.397 )
                    * std::pow( cp_m / cp_b,   0.351 )
                    * std::pow( p / pc ,       0.042 );

            // Stanton Number
            real St = Nu / ( mReDh * Pr );

            // Reynolds-Colburn ( 0.647 = 1 - 0.353 )
            mSigma = std::pow( Pr,            0.647 );

            // Recovery factor is neglected
            mRecovery = 0 ;

            // recovery enthalpy
            aHr = mHm ; //+ 0.5 * mRecovery * mUm * mUm ;

            // heat load
            aDotQ = St * rho_b * mUm * ( aHr - hw );

            // friction factor
            real cf = 2.0 * mSigma * St ;

            // shear stress
            aTauw = 0.5 * cf * rho_b * mUm * mUm ;
        }

 //------------------------------------------------------------------------------

        void
        Boundarylayer::friction_lebedinsky_kalmykov( real & aDotQ, real & aTauw, real & aHr, real & aTr )
        {
            BELFEM_ASSERT( mGas.component( 0 )->label() == "CH4" && mGas.number_of_components() == 1,
                          "the Lebedinsky-Kalmykov correlation can only be used with methane" );

            const real & Tw = this->Tw() ;

            real Pr = mGas.Pr( mTm, mP );

            // equation for Methane
            real Nu = 0.0185 * std::pow( mReDh, 0.8 ) * std::pow( Pr, 0.4 )
                      * std::pow( mTm / Tw, 0.1 );
            real St = Nu / ( mReDh * Pr );

            aDotQ = St * mRhom * mUm * ( mHm - mGas.h( Tw, mP ) );

            // help magnitude for cf
            real Xcf =  2.0 / ( mRhom * mUm * mUm );

            aHr = mHm ;
            aTr = mTm ;

            uint tCount = 0 ;

            aTauw = 2.0 * mSigma * St * Xcf;

            real tau_w_old = BELFEM_REAL_MAX ;

            while( std::abs( tau_w_old - aTauw ) / aTauw > 1e-6 )
            {
                // shift tau
                tau_w_old = aTauw;

                this->compute_sigma_recovery_petrukov();

                aTauw *= 0.1;
                aTauw += 0.9 * 2.0 * mSigma * St * Xcf;

                BELFEM_ERROR( tCount++ < 100, "Too many iterations" );
            }

            // recovery factor is ignored here
            mRecovery = 0.0 ;
        }

//------------------------------------------------------------------------------
/*
        void
        Boundarylayer::friction_vandriest(
                real & aDotQ,
                real & aTauw,
                real & aHr,
                real & aTr )
        {
            // note: for better readability, we don't use the usual
            //       notation for member objects and temporary objects
            //       in this context

            // relaxation factor
            real relax = 0.3 ;

            // create shortcuts

            // heat flux
            real & dotQ = aDotQ ;

            // shear stress at wall
            real & tau_w = aTauw ;

            // old stress
            real tau_w_old = 0.0 ;

            // recovery temperature
            real & h_r = aHr ;

            // recovery temperature
            real & T_r = aTr ;

            // recovery factor
            real & r = mRecovery ;

            // Reynolds-Colburn analogy
            real & sigma = mSigma ;

            // Crocco-Buseman parameters
            cplx & phi = mPhi ;
            cplx & psi = mPsi ;
            cplx & chi = mChi ;

            // offset for log-law
            real & B = mBplus ;

            // Kármán constant
            const real & k = mKarman ;

            const real & Pr_T = mPrTinf ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // we begin with computin the freestream condition
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // freestream temperature
            const real & T = mTm ;

            // freestream pressure
            const real & p = mP ;

            // freestream velocity
            const real & u = mUm ;

            // density
            real rho = mGas.rho( T, p );

            // enthalpy
            real h = mGas.h( T, p );

            // viscosity
            real mu = mGas.mu( T, p );

            // mach number
            real Ma = mUm / mGas.c( T, p );

            // adiabatic coefficient
            real gamma = mGas.gamma( T, p );


            // Prandtl Number
            real Pr = mGas.Pr( T, p );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Now the properties at the wall
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // the wall temperature
            const real & T_w = this->Tw();

            // the density at the wall
            real rho_w = mGas.rho( T_w, p );

            // specific heat at the wall
            //real cp_w = mGas.cp( T_w, p );

            // enthalpy at the wall
            real h_w = mGas.h( T_w, p );

            // viscosity at the wall
            real mu_w = mGas.mu( T_w, p );

            // thermal conductivity at the wall
            real lambda_w = mGas.lambda( T_w, p );


            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Now the Reynolds number
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // reference position for Reynolds
            real x = mXref + mXoff ;

            // Reynolds number at this point
            real Re_x = rho * u * x / mu ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Displacement Thickness
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // based on Hirschel: Selected Aerothermodynamic Design Problems of
            // Hypersonic flight vecicles

            real delta_1_ic = 0.0504 * x / std::pow( Re_x, 0.2 );

            // help magnitude
            real delta_help = delta_1_ic * (
                      0.129 + 0.871 * T_w/T + 0.648 * 0.5 * ( gamma - 1. ) * Ma * Ma ) ;

            // reference enthalpy according to eckert
            real h_ref ;

            // reference temperature according to Eckert
            real T_ref ;

            // displacement thickness
            real delta_1 ;

            real beta ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Help variables for direst coefficents
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            const real xi   = Pr / mPrTinf - 1.0 ;
            const real eta  = std::log( 1. + 5./6. * xi );
            const real zeta = 1. + 0.875 * xi ;
            const real epsilon
            = constant::pi * constant::pi / 6. + 1.5 * ( 1. - mPrTinf );
            const real ln6 = std::log( 6. );
            real alpha ;


            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Initial guesses for the other values
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // Using Ackermann's equation for first guess
            real r0 = std::pow( mGas.Pr( T_w, p ), 1./3. );
            r = r0 ;

            // first guess for Reynolds-Colburn Analogy
            sigma = r * r ;

            // Coles Wake Parameter
            real & Pi = mPi_driest ;

            // initial guess for friction parameter
            real cf = 0.4768 / std::pow( std::log10( Re_x ), 2.591 );
            tau_w = 0.5 * cf * rho * u * u ;

            // initial guess for log-law offset
            B = mBplus0 ;

            // conversion factor for Reynolds number
            real F_Re_x ;

            // conversion factor for friction parameter
            real F_cf ;

            // help factor
            real G ;

            // term for roughness
            real K ;

            // parameter for power law
            real omega = std::log( mu / mu_w ) / std::log( T / T_w ) ;

            // initialize counter
            uint tCount = 0 ;

            // shear velocity
            real u_tau ;

            real cf_inc ;


            // begin loop
            while( std::abs( tau_w - tau_w_old ) > 1e-7 )
            {
                // update recovery enthalpy
                h_r = h + 0.5 * r * u * u;

                // compute recovery temperature
                T_r = mGas.T_from_h( h_r, p );

                // Crocco-Busemann parameters
                psi = mu_w * ( h_r - h_w ) * rho_w * mGas.alpha( T_w, p ) /
                      ( lambda_w * sigma );

                phi = 1.0 + psi - rho_w / rho;

                chi = std::sqrt( psi * psi + 4. * phi );

                // van Driest transformation factors
                G = std::real( std::asin(( 2. * phi - psi ) / chi ) + std::asin( psi / chi ));
                F_cf = ( T_r / T - 1. ) / ( G * G );
                F_Re_x = mu / ( mu_w * F_cf );



                // reference temperature according to Eckert
                // Parameters from Meador and Smart
                h_ref = 0.34 * h + 0.16 * h_r + 0.5 * h_w;
                T_ref = mGas.T_from_h( h_ref, p );

                // displacement thickness
                delta_1 = delta_help * std::pow( T_ref / T, 0.2 * ( omega - 4. ));

                // Clauser parameter
                beta = delta_1 * mdPdX / tau_w;

                // Wake parameter, using Christian's correlation
                // Pi *= 0.1 ;
                // Pi += 0.9 * ( 0.55 + beta * ( 0.52049 - 0.017960 * beta ) );

                // compute incompressible friction factor
                cf_inc = boundarylayer::cf_flatplate_inc_turbulent( Re_x * F_Re_x, B, k, Pi );
                cf = cf_inc / F_cf ;


                // shift shear stress
                tau_w_old = tau_w;

                // compute shear stress
                tau_w *= ( relax - 1.0 );
                tau_w += relax * 0.5 * cf * rho * u * u;

                // help factor for sigma and recovery
                if( tCount < 100 )
                {
                    alpha = std::sqrt( cf * 0.5 );

                    // Driest 1954, Eq. ( 41 )
                    sigma = Pr_T * ( 1.0 + 5.0 * alpha * ( 0.2 / k * ( 1.0 - Pr_T ) * epsilon + xi + eta ));

                    // Driest 1954, Eq. ( 54 )
                    r = Pr_T * ( 1.0 + 2. / k * alpha * ( 1.0 - Pr_T ) * epsilon + 12.5 * cf
                                                                                   * (( xi + 2.0 * eta +
                                                                                        ln6 * std::log( zeta ))));
                }
                else
                {
                    r = r0 ;
                    sigma = r * r;
                }

                std::cout << tCount << " " << x << " " << cf << " " << r << " " << sigma << " " << beta << " " << Pi << std::endl ;

                // not sure if cf or cf_inc is correct here
                u_tau = u * std::sqrt( 0.5 * cf_inc );

                K = rho_w * mKtech * u_tau / mu_w;

                B = mBplus0 - std::log( 1.0 + K / mRoughConst ) / mKarman;

                BELFEM_ERROR( tCount++ < 200, "Too many iterations" );
            }

            // geat flux
            dotQ = tau_w * ( h_r - h_w ) / ( sigma * u );
        } */

//------------------------------------------------------------------------------

        real
        Boundarylayer::check_u_hat_from_t_hat( const real & aThat )
        {
            // write data into containers
            if( mGas.is_idgas() )
            {
                mData( mCenter, BELFEM_CHANNEL_RHO ) = mGas.rho( aThat, mP );
            }
            else
            {
                mData( mCenter, BELFEM_CHANNEL_RHO ) = 1.0 / mVolumeSpline->eval( aThat );
            }

            mData( mCenter, BELFEM_CHANNEL_H ) = mHeatSpline->eval( aThat );
            mData( mCenter, BELFEM_CHANNEL_T ) = aThat ;

            const real & u_hat    = mData( mCenter, BELFEM_CHANNEL_U ) ;
            const real & h_hat    = mData( mCenter, BELFEM_CHANNEL_H ) ;
            const real & rho_hat  = mData( mCenter, BELFEM_CHANNEL_RHO ) ;

            const real & T_w      = mData( 0, BELFEM_CHANNEL_T) ;

            const real & rho_w    = mData( 0, BELFEM_CHANNEL_RHO ) ;

            const real & h_w      = mData( 0, BELFEM_CHANNEL_H ) ;
            const real & mu_w     = mData( 0, BELFEM_CHANNEL_MU ) ;
            const real & lambda_w = mData( 0, BELFEM_CHANNEL_LAMBDA ) ;

            // ( 21 )
            mPsi = ( mu_w / ( mSigma * lambda_w ) )
                   * ( h_hat + mRecovery * 0.5 * u_hat * u_hat - h_w )
                   * mGas.alpha( T_w, mP );

            // ( 20 )
            mPhi = 1.0 + mPsi - rho_w / rho_hat;

            mChi = std::sqrt(  mPsi * mPsi + 4.0 * mPhi );

            mBeta = std::asin( mPsi / mChi );

            this->compute_shear_stress();

            mAlpha = std::sqrt( mPhi ) * mUtau / u_hat;

            real tYplus = rho_w * mUtau * mData( mCenter, BELFEM_CHANNEL_Y ) / mu_w ;

            real tFplus = 1.0 / mKarman * std::log( tYplus ) + mBplus ;

            tFplus = boundarylayer::spalding( mBplus, mKarman, mExpKB, tYplus, tFplus );

            // compute velocity function
            real tUplus = tFplus +  boundarylayer::g_plus( mKarman, mPi, 1.0 );

            return std::real( ( std::sin( tUplus * mAlpha - mBeta ) * mChi + mPsi )
                    / ( 2.0 * mPhi )) - 1.0 ;
        }

//------------------------------------------------------------------------------
    }
}