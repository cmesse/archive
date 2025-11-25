#include "cl_CH_GeometryCylinderCombustor.hpp"


namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        GeometryCylinderCombustor::GeometryCylinderCombustor(
                const real & aThroatDiameter,
                const real & aChamberDiameter,
                const real & aCylinderLength,
                const real & aChamberLength,
                const real & aKinkRadius,
                const real & aCurvatureRadius ) :
            Geometry( true )
        {
            this->set_length( aChamberLength );

            mR  = 0.5 * aChamberDiameter ;
            mRk = aKinkRadius ;
            mRc = aCurvatureRadius ;

            // center point for kink circle, Point 3
            mKx = aCylinderLength ;
            mKr = 0.5 * aChamberDiameter - aKinkRadius ;

            // center point for circle at throat, Point 6
            mMx = aChamberLength ;
            mMr = 0.5 * aThroatDiameter + aCurvatureRadius ;

            // help lengths
            real tDX = aChamberLength - aCylinderLength ;
            real tDR = mMr - mKr ;
            real tSR = aKinkRadius + aCurvatureRadius ;

            // contraction angle
            real tAlpha = std::asin( tSR / tDX
                    - ( tDR * ( tDX * std::sqrt( tDR * tDR - tSR * tSR + tDX * tDX )
                       + tDR * tSR ) ) / ( tDX * ( tDR * tDR + tDX * tDX ) ) );


            real tS = std::sin( tAlpha ) ;
            real tC = std::cos( tAlpha ) ;

            // Point 4
            mPx = mKx + tS * aKinkRadius ;
            mPr = mKr + tC * aKinkRadius ;

            // Point 5
            mQx = mMx - tS * aCurvatureRadius ;
            mQr = mMr - tC * aCurvatureRadius ;

            // slope
            mA = ( mQr - mPr ) / ( mQx - mPx ) ;

            // offset
            mB = mPr - mA * mPx ;
        }

//------------------------------------------------------------------------------

        real
        GeometryCylinderCombustor::R( const real & aX ) const
        {
            if ( aX <= mKx )
            {
                return mR ;
            }
            else if ( aX < mPx )
            {
                return mKr + std::sqrt( ( mRk + mKx - aX ) * ( mRk - mKx + aX ) );
            }
            else if( aX < mQx )
            {
                return mA * aX + mB ;
            }
            else if ( aX < mMx )
            {
                return mMr - std::sqrt( ( mRc + mMx - aX ) * ( mRc - mMx + aX ) );
            }
            else
            {
                return mMr - mRc ;
            }
        }

//------------------------------------------------------------------------------

        real
        GeometryCylinderCombustor::dRdx( const real & aX ) const
        {
            if ( aX <= mKx )
            {
                return 0.0 ;
            }
            else if ( aX < mPx )
            {
                return ( mKx - aX ) / std::sqrt( ( mRk + mKx - aX ) * ( mRk - mKx + aX ) );
            }
            else if( aX < mQx )
            {
                return mA ;
            }
            else if ( aX < mMx )
            {
                return ( aX - mMx ) / std::sqrt( ( mRc + mMx - aX ) * ( mRc - mMx + aX ) );
            }
            else
            {
                return 0.0 ;
            }
        }

//------------------------------------------------------------------------------

        real
        GeometryCylinderCombustor::A( const real & aX ) const
        {
            return constant::pi * std::pow( this->R( aX ), 2 );
        }

//------------------------------------------------------------------------------

        real
        GeometryCylinderCombustor::dAdx( const real & aX ) const
        {
            return 2.0 * constant::pi * this->R( aX )  * this->dRdx( aX ) ;
        }

//------------------------------------------------------------------------------

        real
        GeometryCylinderCombustor::P( const real & aX ) const
        {
            return  2.0 * constant::pi * this->R( aX );
        }

//------------------------------------------------------------------------------

        real
        GeometryCylinderCombustor::dPdx( const real & aX ) const
        {
            return  2.0 * constant::pi * this->dRdx( aX );
        }

//------------------------------------------------------------------------------

        real
        GeometryCylinderCombustor::r( const real & aX ) const
        {
            return 0.0 ;
        }

//------------------------------------------------------------------------------

        real
        GeometryCylinderCombustor::drdx( const real & aX ) const
        {
            return 0.0 ;
        }

//------------------------------------------------------------------------------
    }
}