//
// Created by Christian Messe on 08.12.20.
//

#include "cl_CH_GeometryNozzle.hpp"
#include "cl_Matrix.hpp"
#include "fn_gesv.hpp"

namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        GeometryNozzle::GeometryNozzle(
                const real & aThroatDiameter,
                const real & aOpeningAngle,
                const real & aExhaustAngle,
                const real & aExpansionRatio,
                const real & aCircleRadius ) :
                Geometry( true ),
                mType( NozzleType::Rao ),
                mThroatDiameter( aThroatDiameter ),
                mOpeningAngle( aOpeningAngle ),
                mExhaustAngle( aExhaustAngle ),
                mExpansionRatio( aExpansionRatio ),
                mCircleRadius( aCircleRadius )
        {
            // compute the tiny circle
            this->compute_rao_circle();

            // compute the rao nozzle
            this->compute_rao_coefficients();

            // link the functions
            mRadiusFunction = & GeometryNozzle::radius_rao ;
            mSlopeFunction  = & GeometryNozzle::dradius_rao ;
        }

//------------------------------------------------------------------------------

        // Constructor for Bezier Nozzle
        GeometryNozzle::GeometryNozzle(
                const real & aThroatDiameter,
                const real & aOpeningAngle,
                const real & aExhaustAngle,
                const real & aExpansionRatio,
                const real & aCircleRadius,
                const real & aLength,
                const real & aXi,
                const real & aEta )  :
                Geometry( true ),
                mType( NozzleType::Bezier ),
                mThroatDiameter( aThroatDiameter ),
                mOpeningAngle( aOpeningAngle ),
                mExhaustAngle( aExhaustAngle ),
                mExpansionRatio( aExpansionRatio ),
                mCircleRadius( aCircleRadius ),
                mEx( aLength )
        {
            // compute the tiny circle
            this->compute_rao_circle();

            this->create_bezier_object( aXi, aEta ) ;


            // link the functions
            mRadiusFunction = & GeometryNozzle::radius_bezier ;
            mSlopeFunction  = & GeometryNozzle::dradius_bezier ;

        }

//------------------------------------------------------------------------------

        GeometryNozzle::~GeometryNozzle()
        {
            if ( mBezier != nullptr )
            {
                delete mBezier ;
            }
        }

//------------------------------------------------------------------------------
        void
        GeometryNozzle::compute_rao_circle()
        {
            // center of nozzle circle
            mMx = 0.0 ;
            mMr = 0.5 * mThroatDiameter + mCircleRadius ;

            // end of nozzle circle
            mNx = mMx + std::sin( mOpeningAngle ) * mCircleRadius ;
            mNr = mMr - std::cos( mOpeningAngle ) * mCircleRadius ;

            // cross section at throat
            real tAt = 0.25 * constant::pi * mThroatDiameter * mThroatDiameter ;

            // cross section at exit
            real tAe = tAt * mExpansionRatio ;

            // radius at exit
            mEr = std::sqrt( tAe / constant::pi );

        }

//------------------------------------------------------------------------------


        void
        GeometryNozzle::compute_rao_coefficients()
        {
            // the Vandermonde matrix
            Matrix< real > tV( 3, 3 );

            tV( 0, 0 ) = mNr * mNr ;
            tV( 1, 0 ) = 2.0 * mNr ;
            tV( 2, 0 ) = 2.0 * mEr ;
            tV( 0, 1 ) = mNr ;
            tV( 1, 1 ) = 1.0 ;
            tV( 2, 1 ) = 1.0 ;
            tV( 0, 2 ) = 1.0 ;
            tV( 1, 2 ) = 0.0 ;
            tV( 2, 2 ) = 0.0 ;

            // the right hand side
            Vector< real > tX( 3 );

            tX( 0 ) = mNx ;
            tX( 1 ) = 1.0 / std::tan( mOpeningAngle );
            tX( 2 ) = 1.0 / std::tan( mExhaustAngle );

            // pivot for gesv
            Vector< int > tP( 3 );

            // solve system
            gesv( tV, tX, tP );

            // write the coefficients for the nozzle
            mA = tX( 0 );
            mB = tX( 1 );
            mC = tX( 2 );

            // compute the length
            mEx = mEr * ( mA * mEr + mB ) + mC ;
        }
//------------------------------------------------------------------------------

        void
        GeometryNozzle::create_bezier_object( const real & aXi, const real & aEta )
        {
            // create the bezier object with arbitrary parameters
            mBezier = new Bezier( 0., 0., 1., 1., 1., 1. );

            // help points
            real tA =  aXi * ( mEx - mNx ) ;
            real tB = ( 1.0 - aEta ) * ( mEx - mNx ) ;

            // manually compute the reference points along the X-axis
            Vector< real > & tX = mBezier->basis_x() ;
            tX( 0 ) = mNx ;
            tX( 1 ) = mNx + tA ;
            tX( 2 ) = mNx + tB ;
            tX( 3 ) = mEx ;

            // compute the reference points along the radius
            Vector< real > & tR = mBezier->basis_y() ;
            tR( 0 ) = mNr ;
            tR( 1 ) = mNr + tA * std::tan( mOpeningAngle );
            tR( 2 ) = mEr - tB * std::tan( mExhaustAngle );
            tR( 3 ) = mEr ;
        }

//------------------------------------------------------------------------------

        real
        GeometryNozzle::radius_circle( const real & aX ) const
        {
            return mMr - std::sqrt(
                    ( mCircleRadius + aX - mMx )
                  * ( mCircleRadius - aX + mMx ) );
        }

        real
        GeometryNozzle::dradius_circle( const real & aX ) const
        {
            return ( aX - mMx ) / std::sqrt(
                       ( mCircleRadius + aX - mMx )
                    *  ( mCircleRadius - aX + mMx ) );
        }

//------------------------------------------------------------------------------

        real
        GeometryNozzle::radius_bezier( const real & aX ) const
        {
            return mBezier->y( aX );
        }

        real
        GeometryNozzle::dradius_bezier( const real & aX ) const
        {
            return mBezier->dydx( aX );
        }

//------------------------------------------------------------------------------

        real
        GeometryNozzle::radius_rao( const real & aX ) const
        {
            return ( - mB + std::sqrt( mB * mB - 4.0 * mA * ( mC - aX ) ) ) /
                    ( 2.0 * mA );
        }

//------------------------------------------------------------------------------

        real
        GeometryNozzle::dradius_rao( const real & aX ) const
        {
            return 1.0 / std::sqrt( mB * mB - 4.0 * mA * ( mC - aX ) ) ;
        }


//------------------------------------------------------------------------------

    }
}