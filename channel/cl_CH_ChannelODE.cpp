//
// Created by Christian Messe on 25.02.20.
//

#include "fn_gesv.hpp"
#include "cl_CH_ChannelODE.hpp"
#include "cl_CH_Element.hpp"
#include "cl_CH_Boundarylayer.hpp"
#include "cl_GT_RefGas.hpp"
#include "fn_dot.hpp"
namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        ChannelODE::ChannelODE(
                Geometry & aGeometry,
                Gas & aGas,
                const ChannelMode aMode ) :
            ODE( 3 ),
            mMode( aMode ),
            mGeometry( &aGeometry ),
            mGas( aGas ),
            mElementMode( false )
        {
            if ( aMode == ChannelMode::Channel )
            {
                // allocate Jacobi matrix
                mJacobi.set_size( 3, 3 );

                // allocate Pivot Vector
                mPivot.set_size( 3 );

                // link compute matrix
                if( mGas.is_idgas() )
                {
                    mJacobiFunction
                            = & ChannelODE::compute_jacobi_idgas ;
                }
                else
                {
                    mJacobiFunction
                            = & ChannelODE::compute_jacobi_realgas ;
                }

                mComputeFunction = & ChannelODE::compute_channel_ode ;
            }
            else
            {
                BELFEM_ERROR( mGas.is_idgas(), "Gas must be ideal gas in combustor mode" );
                mComputeFunction = & ChannelODE::compute_combustor_ode ;
            }

            // link geometry function
            mGeometryFunction
                = & ChannelODE::compute_geometry_from_geometry_object ;

            mdYdx.set_size( mGas.number_of_components(), 0.0 );


            // link friction function
            mFrictionFunction
                = & ChannelODE::compute_friction_dittus_boelter ;

        }

//------------------------------------------------------------------------------

        // constructor
        ChannelODE::ChannelODE( Gas & aGas,
                const ChannelMode aMode  ):
                ODE( 3 ),
                mMode( aMode ),
                mGas( aGas ),
                mElementMode( true )
        {
            if ( aMode == ChannelMode::Channel )
            {
                // allocate Jacobi matrix
                mJacobi.set_size( 3, 3 );

                // allocate Pivot Vector
                mPivot.set_size( 3 );

                // link compute matrix
                if( mGas.is_idgas() )
                {
                    mJacobiFunction
                            = & ChannelODE::compute_jacobi_idgas ;
                }
                else
                {
                    mJacobiFunction
                            = & ChannelODE::compute_jacobi_realgas ;
                }

                mComputeFunction = & ChannelODE::compute_channel_ode ;
            }
            else
            {
                BELFEM_ERROR( mGas.is_idgas(), "Gas must be ideal gas in combustor mode" );
                mComputeFunction = & ChannelODE::compute_combustor_ode ;
            }

            mWorkN.set_size( 3 );
            mWorkV.set_size( 3 );

            // link geometry function
            mGeometryFunction
                    = & ChannelODE::compute_geometry_from_element ;

            // link friction function
            mFrictionFunction
                    = & ChannelODE::compute_friction_element ;
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::set_wall_temperature( const real & aTw )
        {
            BELFEM_ASSERT( ! mElementMode,
                          "can't set wall temperature when ODE is in Element-Mode");

            mTwall = aTw ;
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::set_combustion( const real & adRdxR, const real & adwdx )
        {
            mdRdxR = adRdxR ;
            mdwdx  = adwdx ;
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::set_composition_change( const real & adRdxR, const Vector< real >  & adYdx )
        {
            mCombust = 0.0 ;
            mdRdxR = adRdxR ;
            mdYdx = adYdx ;
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::link_element( Element * aElement )
        {
            // todo: to make this an assert, element mode must vanish in non-debug mode
            BELFEM_ERROR( mElementMode , "can't link when ODE is in Geometry-Mode" );
            mElement = aElement ;
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::link_geometry( Geometry * aGeometry, const bool aReverse )
        {
            mGeometry = aGeometry ;

            // link geometry function
            if( aReverse )
            {
                mReverse = true ;
                mGeometryFunction
                        = & ChannelODE::compute_geometry_from_geometry_object_reversex ;
            }
            else
            {
                mReverse = false ;
                mGeometryFunction
                        = & ChannelODE::compute_geometry_from_geometry_object;
            }
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute( const real & aT,
                             const Vector< real > & aY,
                             Vector< real > & adYdT )
        {
            ( this->*mComputeFunction )( aT, aY, adYdT );
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute_channel_ode(
                const real          & aT,
                const Vector <real> & aY,
                      Vector <real> & adYdT )
        {
            // get inverse density
            const real & tV = aY( 0 );

            // get velocity
            const real & tU = aY( 1 );

            // get temperature
            const real & tT = aY( 2 );

            // compute pressure
            const real tP = mGas.p( tT, tV );

            // compute geometry data
            ( this->*mGeometryFunction )( aT, mDh, mA, mdAdx );

            // mass flux
            real tDotM = mA * tU / tV ;

            // local shear stress
            real tTauW ;

            // heat flux
            real tDotQ ;

            // Compute Jacobian Matrix
            ( this->*mJacobiFunction )( tV, tU, tT, tP );

            // compute shear stress
            ( this->*mFrictionFunction ) (
                    aT,        // x
                    tV,   // v
                    tU,   // u
                    tT,   // T
                    tP,   // p
                    mTwall, // wall
                    tTauW,
                    tDotQ ) ;

            // assemble right hand side
            adYdT( 0 ) = mdAdx / mA ;
            adYdT( 1 ) = - 4.0 * tTauW / ( mDh * tP ) - mdRdxR ;
            adYdT( 2 ) = - 4.0 * mA * tDotQ / ( mDh * tDotM ) - mdwdx ;

            if( mCombust )
            {
                this->compute_combustion( tT, tP, adYdT );
            }
            if( mReverse )
            {
                adYdT *= -1.0 ;
            }

            // solve system
            gesv( mJacobi, adYdT, mPivot );

            // scale output
            adYdT( 0 ) *= aY( 0 );
            adYdT( 1 ) *= aY( 1 );
            adYdT( 2 ) *= aY( 2 );
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute_combustor_ode(
                               const real          & aT,
                               const Vector <real> & aY,
                               Vector <real>       & adYdT )
        {
            // get inverse density
            const real & v = aY( 0 );

            // get velocity
            const real & u = aY( 1 );

            // get temperature
            const real & T = aY( 2 );

            real & dvdx = adYdT( 0 );
            real & dudx = adYdT( 1 );
            real & dTdx = adYdT ( 2 );

            // pressure
            real p = mGas.p( T, v );

            // compute geometry data
            ( this->*mGeometryFunction )( aT, mDh, mA, mdAdx );

            // shear stress
            real tau_w ;

            // heatflux
            real dot_q ;

            // compute shear stress
            ( this->*mFrictionFunction ) (
                    aT,        // x
                    v,   // v
                    u,   // u
                    T,   // T
                    p,   // p
                    mTwall, // wall
                    tau_w,
                    dot_q ) ;

            // friction factor
            real cf = 2.0 * v * v * tau_w / ( u * u );

            // specific heat capacity
            real cp = mGas.cp( T, p );

            // ratio of specific heats
            real k = mGas.gamma( T, p );

            // speed of sound
            real c = mGas.c( T, p );

            // gas constant
            real R = mGas.R( T, p );

            // Mach number
            real Ma = u / c ;

            // squared mach number
            real Ma2 = Ma * Ma ;

            // combustion part
            if( mCombust )
            {
                // perimeter
                real U = 4.0 * mA / mDh ;

                mWorkV.fill( 0.0 );
                this->compute_combustion( T, p, mWorkV );
                mdwdx = mWorkV( 2 )-dot_q * U * v / ( mA * u );
            }

            // help variable
            real xi = ( 1. + k * Ma2 ) * mdMdxM
                     + k * ( 2.0 * cf / mDh * Ma2 - mdIdx )
                     + mdRdxR - mdAdx / mA ;

            // help variable
            real eta = mdwdx / ( cp * T );

            // temperature change
            dTdx = T / ( Ma2 - 1.0 ) *
                    ( ( k-1.) * xi * Ma2 + ( k * Ma2 - 1 ) * eta );

            // change in Mach Number
            real dMadx = Ma * ( ( xi + eta ) / ( 1. - Ma2 ) - 0.5 * ( dTdx/T + mdRdxR ) );

            // change in the ration of specific heats
            real dcpdT = mGas.dcpdT( T, p );
            real dkdT = ( dcpdT - cp * dcpdT / ( cp - R ) ) / ( cp - R );
            real dkdx = dkdT * dTdx ;

            // change in the speed of sound
            real dcdx = ( R * T * dkdx + c * c * mdRdxR + k * R * dTdx ) / ( 2. * c );

            // change in velocity
            dudx = dMadx * c + Ma * dcdx ;

            // change in inverse density
            dvdx = v * ( mdAdx / mA + dudx / u - mdMdxM );

            if ( mReverse )
            {
                adYdT *= -1 ;
            }
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute_geometry_from_geometry_object(
                const real & aX,
                      real & aDh,
                      real & aA,
                      real & adAdX )
        {
            aDh   = mGeometry->Dh( aX );
            aA    = mGeometry->A( aX );
            adAdX = mGeometry->dAdx( aX );
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute_geometry_from_geometry_object_reversex(
                const real & aX,
                real & aDh,
                real & aA,
                real & adAdX )
        {
            real tX = mGeometry->length();

            aDh   = mGeometry->Dh( tX-aX );
            aA    = mGeometry->A( tX-aX );
            adAdX = mGeometry->dAdx( tX-aX );

        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute_geometry_from_element(
                const real & aX,
                      real & aDh,
                      real & aA,
                      real & adAdX )
        {
            // compute the interpolation function
            mElement->compute_N( aX, mWorkN );

            // Hydraulic Diameter
            mElement->collect_data( BELFEM_CHANNEL_DH, mWorkV );
            aDh = dot( mWorkN, mWorkV );

            // Cross Section
            mElement->collect_data( BELFEM_CHANNEL_A, mWorkV );
            aA = dot( mWorkN, mWorkV );

            // Cross section derivative
            mElement->compute_B( aX, mWorkN );
            adAdX = dot( mWorkN, mWorkV );
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute_jacobi_idgas(
                const real & aV,
                const real & aU,
                const real & aT,
                const real & aP )
        {
            // see 10.2514/6.2017-4989, Eq. 20
            mJacobi( 0, 0 ) =  1.0;
            mJacobi( 1, 0 ) = -1.0;
            mJacobi( 2, 0 ) =  0.0;

            mJacobi( 0, 1 ) = -1.0;
            mJacobi( 1, 1 ) = aU * aU / ( aP * aV );
            mJacobi( 2, 1 ) = aU * aU ;

            mJacobi( 0, 2 ) = 0.0;
            mJacobi( 1, 2 ) = 1.0;
            mJacobi( 2, 2 ) = mGas.cp( aT, aP ) * aT ;
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute_jacobi_realgas(
                const real & aV,
                const real & aU,
                const real & aT,
                const real & aP )
        {
            const real tAlpha = mGas.alpha( aT, aP );
            const real tBeta  = mGas.beta( aT, aP );
            const real tKappa = mGas.kappa( aT, aP );

            // see 10.2514/6.2017-4989, Eq. 17

            mJacobi( 0, 0 ) =  1.0 ;
            mJacobi( 1, 0 ) = - tBeta / tAlpha ;
            mJacobi( 2, 0 ) = ( aT * tAlpha - 1.0 ) * aV / tKappa ;

            mJacobi( 0, 1 ) =  -1.0 ;
            mJacobi( 1, 1 ) = aU * aU / ( aP * aV );
            mJacobi( 2, 1 ) = aU * aU ;

            mJacobi( 0, 2 ) =  0.0 ;
            mJacobi( 1, 2 ) = tBeta * aT ;
            mJacobi( 2, 2 ) = ( mGas.cv( aT, aP ) + aP * aV * tBeta ) * aT ;
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute_friction_dittus_boelter(
                const real & aX,
                const real & aV,
                const real & aU,
                const real & aT,
                const real & aP,
                const real & aTwall,
                real & aTauW,
                real & aDotQ )
        {
                // Reynolds number
                real tReDh = aU * mDh /
                        ( aV * mGas.mu( aT, aP ) );

                // recovery factor
                real tRecovery = std::pow( mGas.Pr( aT, aP ), 1.0 / 3.0 );

                // Reynolds-Colburn analogy
                real tSigma = tRecovery * tRecovery ;

                // shear stress
                aTauW = 0.023 * aU * aU /
                        ( std::pow( tReDh, 0.2 ) * aV );

                // heat flux
                aDotQ = aTauW *
                        ( mGas.h( aT, aP ) + 0.5 * tRecovery * aU * aU
                          - mGas.h( aTwall, aP ) ) / ( tSigma * aU );
        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute_friction_element(
                const real & aX,
                const real & aV,
                const real & aU,
                const real & aT,
                const real & aP,
                const real & aTwall,
                real & aTauW,
                real & aDotQ )
        {
            // compute the interpolation function
            mElement->compute_N( aX, mWorkN );

            // shear stress
            mElement->collect_data( BELFEM_CHANNEL_TAUW, mWorkV );
            aTauW = dot( mWorkN, mWorkV );

            // Heat Flux
            mElement->collect_data( BELFEM_CHANNEL_DOTQ, mWorkV );
            aDotQ = dot( mWorkN, mWorkV );

            // change of gas constant
            //mElement->collect_data( BELFEM_CHANNEL_RM, mWorkV );

            //mElement->compute_B( aX, mWorkN );
            //mdRdxR = dot( mWorkN, mWorkV ) / mWorkV( 2 );

        }

//------------------------------------------------------------------------------

        void
        ChannelODE::compute_combustion( const real & aT, const real & aP, Vector<real> & aRHS )
        {
            uint tN = mGas.number_of_components() ;

            // molar mass, does not depend from T and p in this case
            real tM = mGas.M( aT, aP );
            for( uint k=0; k<tN; ++k )
            {
                gastables::RefGas * tGas = mGas.component( k );

                aRHS( 1 ) -= tM / tGas->M() * mdYdx( k );
                aRHS( 2 ) -= tGas->h( aT ) * mdYdx( k );
            }

        }

//------------------------------------------------------------------------------
    }
}