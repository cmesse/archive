//
// Created by Christian Messe on 24.11.19.
//

#include "cl_CH_GeometryCombustor.hpp"

namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        GeometryCombustor::GeometryCombustor() :
            Geometry( false )
        {
            this->initialize();
        }

//------------------------------------------------------------------------------

        void
        GeometryCombustor::initialize()
        {
            // update length
            mLength = mLength0 + mLength1;

            // Major Points for geometry
            mPoints.set_size( 2, 9 );

            // beginning of combustor
            mPoints( 0, 0 ) = 0.0;
            mPoints( 1, 0 ) = 0.5 * mHeight0;

            // end of part of constant height
            mPoints( 0, 1 ) = mLength0;
            mPoints( 1, 1 ) = mPoints( 1, 0 );

            // supporting point for kink
            mPoints( 0, 4 ) = mPoints( 0, 1 );
            mPoints( 1, 4 ) = mPoints( 1, 1 ) + mKinkRadius;

            // endpoint for rounded kink
            mPoints( 0, 2 ) = mPoints( 0, 4 ) + mKinkRadius * std::sin( mAngle );
            mPoints( 1, 2 ) = mPoints( 1, 4 ) - mKinkRadius * std::cos( mAngle );

            mLines.set_size( 2, 4 );

            // first line
            mLines( 0, 0 ) = 0.0;
            mLines( 1, 0 ) = mPoints( 1, 0 );

            // second line
            mLines( 0, 1 ) = std::tan( mAngle );
            mLines( 1, 1 ) = mPoints( 1, 2 ) - mPoints( 0, 2 ) * mLines( 0, 1 );

            // endpoint for combustor
            mPoints( 0, 3 ) = mLength0 + mLength1;
            mPoints( 1, 3 ) = mLines( 0, 1 ) * mPoints( 0, 3 ) + mLines( 1, 1 );

            // point for beginning of injector
            mPoints( 0, 5 ) = mInjectorPosition - mInjectorLength;
            mPoints( 1, 5 ) = 0.0;

            // ramp of injector
            mLines( 0, 2 ) = mInjectorHeight / mInjectorLength;
            mLines( 1, 2 ) = mPoints( 1, 5 ) - mPoints( 0, 5 ) * mLines( 0, 2 );

            // top of injector
            mPoints( 0, 6 ) = mInjectorPosition - 0.5 * mInjectorLength;
            mPoints( 1, 6 ) = 0.5 * mInjectorHeight;

            // end of injector
            mPoints( 0, 7 ) = mInjectorPosition;
            mPoints( 1, 7 ) = 0.5 * mInjectorHeight;

            mPoints( 0, 8 ) = mInjectorPosition;
            mPoints( 1, 8 ) = 0;

            mLines( 0, 3 ) = 0.0;
            mLines( 1, 3 ) = mPoints( 1, 7 );
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::R( const real & aX ) const
        {
            if( aX < mPoints( 0, 1 ) )
            {
                return aX * mLines( 0, 0 ) + mLines( 1, 0 );
            }
            else if ( aX < mPoints( 0, 2 ) )
            {
                return mPoints( 1, 4 ) - std::sqrt(  ( mKinkRadius + aX - mPoints( 0, 4 ) )
                   * ( mKinkRadius - aX + mPoints( 0, 4 ) ) );
            }
            else
            {
                return aX * mLines( 0, 1 ) + mLines( 1, 1 );
            }
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::dRdx( const real & aX ) const
        {
            if( aX < mPoints( 0, 1 ) )
            {
                return mLines( 0, 0 );
            }
            else if ( aX < mPoints( 0, 2 ) )
            {
                return ( aX - mPoints( 0, 4 ) ) /
                    std::sqrt(  ( mKinkRadius + aX - mPoints( 0, 4 ) )
                              * ( mKinkRadius - aX + mPoints( 0, 4 ) ) );
            }
            else
            {
                return mLines( 0, 1 );
            }
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::r( const real & aX ) const
        {
            return -this->R( aX );
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::drdx( const real & aX ) const
        {
            return  -this->dRdx( aX );
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::r_injector( const real & aX ) const
        {
            if( aX <= mPoints( 0, 5 ) )
            {
                return 0.0;
            }
            else if ( aX < mPoints( 0, 6 ) )
            {
                return mLines( 0, 2 ) * aX + mLines( 1, 2 );
            }
            else if( aX <= mInjectorPosition )
            {
                return mLines( 0, 3 ) * aX + mLines( 1, 3 );
            }
            else
            {
                return 0.0 ;
            }
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::drdx_injector( const real & aX ) const
        {
            if( aX <= mPoints( 0, 5 ) )
            {
                return 0.0;
            }
            else if ( aX < mPoints( 0, 6 ) )
            {
                return mLines( 0, 2 );
            }
            else if( aX <= mInjectorPosition )
            {
                return mLines( 0, 3 );
            }
            else
            {
                return 0.0 ;
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        GeometryCombustor::points() const
        {
            return mPoints;
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::injector_position() const
        {
            return mInjectorPosition ;
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::A( const real & aX ) const
        {
            if( aX < mInjectorPosition-mInjectorLength || aX > mInjectorPosition )
            {
                return mWidth * 2.0 * this->R( aX ) + ( constant::pi - 4.0 ) * mCornerRadius * mCornerRadius ;
            }
            else
            {
                return mWidth * 2.0 * ( this->R( aX ) - this->r_injector( aX ) )
                       + ( constant::pi - 4.0 ) * mCornerRadius * mCornerRadius ;
            }
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::dAdx( const real & aX ) const
        {
            if( aX < mInjectorPosition-mInjectorLength  || aX > mInjectorPosition )
            {
                return mWidth * 2.0 * this->dRdx( aX );
            }
            else
            {
                return mWidth * 2.0 * ( this->dRdx( aX ) - this->drdx_injector( aX ) );
            }
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::P( const real & aX ) const
        {
            if( aX < mInjectorPosition-mInjectorLength  || aX > mInjectorPosition )
            {
                return 2.0 * mWidth + 4.0 * this->R( aX ) + ( 2.0 * constant::pi - 8.0 ) * mCornerRadius;
            }
            else
            {
                return 2.0 * mWidth + 4.0 * ( this->R( aX ) - this->r_injector( aX ) )
                    + ( 2.0 * constant::pi - 8.0 ) * mCornerRadius;
            }
        }

//------------------------------------------------------------------------------

        real
        GeometryCombustor::dPdx( const real & aX ) const
        {
            if( aX < mInjectorPosition-mInjectorLength  || aX > mInjectorPosition )
            {
                return 4.0 * this->dRdx( aX );
            }
            else
            {
                return 4.0 * ( this->dRdx( aX ) - this->drdx_injector( aX ) );
            }
        }

//------------------------------------------------------------------------------
    }
}
