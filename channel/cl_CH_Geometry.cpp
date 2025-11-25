//
// Created by Christian Messe on 24.11.19.
//


#include "cl_CH_Geometry.hpp"
#include "assert.hpp"
#include "constants.hpp"

namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        Geometry::Geometry( const bool aIsAxisymmetric ) :
            mIsAxisymmetric( aIsAxisymmetric )
        {
            if ( mIsAxisymmetric )
            {
                mFP    = & Geometry::P_axisymmetric;
                mFdPdx = & Geometry::dPdx_axisymmetric;
                mFp    = & Geometry::p_axisymmetric;
                mFdpdx = & Geometry::dpdx_axisymmetric;

                mFA    = & Geometry::A_axisymmetric;
                mFdAdx = & Geometry::dAdx_axisymmetric;

                mHaveSecondWall = false;

            }
            else
            {
                mFP    = & Geometry::P_plane;
                mFdPdx = & Geometry::dPdx_plane;
                mFp    = & Geometry::p_plane;
                mFdpdx = & Geometry::dpdx_plane;

                mFA    = & Geometry::A_plane;
                mFdAdx = & Geometry::dAdx_plane;

                mHaveSecondWall = true;
            }
        }

//------------------------------------------------------------------------------

        bool
        Geometry::is_axisymmetric() const
        {
            return mIsAxisymmetric;
        }

//------------------------------------------------------------------------------

        const real &
        Geometry::width() const
        {
            BELFEM_ASSERT(
                    ! this->is_axisymmetric(),
                    "Geometry::width() can not be called for an axisymmetric geometry" );

            return mWidth;
        }

//------------------------------------------------------------------------------

        const real &
        Geometry::length() const
        {
            return mLength;
        }

//------------------------------------------------------------------------------

        real
        Geometry::R( const real & aX ) const
        {
            BELFEM_ERROR( false, "Geometry::R not implemented for this channel geometry" );
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        Geometry::r( const real & aX ) const
        {
            return -this->R( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::dRdx( const real & aX ) const
        {
            BELFEM_ERROR( false, "Geometry::dRdx not implemented for this channel geometry" );
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        Geometry::drdx( const real & aX ) const
        {
            return -this->dRdx( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::P( const real & aX ) const
        {
            return ( this->*mFP )( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::dPdx( const real & aX ) const
        {
            return ( this->*mFdPdx )( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::p( const real & aX ) const
        {
            return ( this->*mFp )( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::dpdx( const real & aX ) const
        {
            return ( this->*mFdpdx )( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::A( const real & aX ) const
        {
            return( this->*mFA ) ( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::dAdx( const real & aX ) const
        {
            return( this->*mFdAdx ) ( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::Dh( const real & aX ) const
        {
            return 4.0 * this->A( aX ) /
                    ( this->P( aX ) + this->p( aX ) );
        }

//------------------------------------------------------------------------------

        real
        Geometry::dDhdx( const real & aX ) const
        {
            real tP = this->P( aX ) + this->p( aX );

            return 4.0 * ( this->dAdx( aX ) - this->A( aX )
                   * ( this->dPdx( aX ) + this->dpdx( aX ) ) / tP ) / tP;
        }

//------------------------------------------------------------------------------

        real
        Geometry::P_axisymmetric( const real & aX ) const
        {
            return 2.0 * constant::pi * this->R( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::P_plane( const real & aX ) const
        {
            return mWidth;
        }

//------------------------------------------------------------------------------

        real
        Geometry::dPdx_axisymmetric( const real & aX ) const
        {
            return 2.0 * constant::pi * this->dRdx( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::dPdx_plane( const real & aX ) const
        {
            return 0.0;
        }

//------------------------------------------------------------------------------

        real
        Geometry::p_axisymmetric( const real & aX ) const
        {
            return 2.0 * constant::pi * this->r( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::p_plane( const real & aX ) const
        {
            return mWidth;
        }

//------------------------------------------------------------------------------

        real
        Geometry::dpdx_axisymmetric( const real & aX ) const
        {
            return 2.0 * constant::pi * this->drdx( aX );
        }

//------------------------------------------------------------------------------

        real
        Geometry::dpdx_plane( const real & aX ) const
        {
            return 0;
        }

//------------------------------------------------------------------------------

        real
        Geometry::A_axisymmetric( const real & aX ) const
        {
            return constant::pi * (
                      std::pow( this->R( aX ), 2 )
                    - std::pow( this->r( aX ), 2 ) );
        }

//------------------------------------------------------------------------------

        real
        Geometry::dAdx_axisymmetric( const real & aX ) const
        {
            return 2.0 * constant::pi * (
                     this->R( aX ) * this->dRdx( aX )
                   - this->r( aX ) * this->drdx( aX ) );
        }

//------------------------------------------------------------------------------

        real
        Geometry::A_plane( const real & aX ) const
        {
            std::cout << " A " << aX << " " << this->R( aX ) << " " << this->r( aX ) << std::endl ;
            return mWidth * ( this->R( aX ) - this->r( aX ) );
        }

//------------------------------------------------------------------------------

        real
        Geometry::dAdx_plane( const real & aX ) const
        {
            return mWidth * ( this->dRdx( aX ) - this->drdx( aX ) );
        }

//------------------------------------------------------------------------------

        bool
        Geometry::has_second_wall() const
        {
            return mHaveSecondWall;
        }

//------------------------------------------------------------------------------

        void
        Geometry::set_length( const real & aLength )
        {
            mLength = aLength ;
        }

//------------------------------------------------------------------------------
    }
}