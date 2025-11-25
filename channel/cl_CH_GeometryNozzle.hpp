//
// Created by Christian Messe on 08.12.20.
//

#ifndef BELFEM_CL_CH_GEOMETRYNOZZLE_HPP
#define BELFEM_CL_CH_GEOMETRYNOZZLE_HPP

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_CH_Geometry.hpp"
#include "cl_Vector.hpp"
#include "cl_Bezier.hpp"
namespace belfem
{
    namespace channel
    {

        enum class NozzleType
        {
            Rao,
            Bezier,
            UNDEFINED
        };

//------------------------------------------------------------------------------

        class GeometryNozzle : public Geometry
        {
            const NozzleType mType ;

            // Throat diameter in m
            const real mThroatDiameter ;

            // opening angle in rad
            const real mOpeningAngle ;

            // exhaust engla in rad
            const real mExhaustAngle ;

            // expansion ratio
            const real mExpansionRatio ;

            // radius of the small circle
            const real mCircleRadius ;

            // center for nozzle circle
            real mMx;
            real mMr;

            // end of nozzle circle
            real mNx;
            real mNr;

            // end point
            real mEx ;
            real mEr ;

            // offset for x-coordinate
            real mXoff = 0.0;

            // polynomial, if this is a rao nozzle
            real mA = BELFEM_QUIET_NAN ;
            real mB = BELFEM_QUIET_NAN ;
            real mC = BELFEM_QUIET_NAN ;

            // the bezier object
            Bezier * mBezier = nullptr ;

//------------------------------------------------------------------------------

            // Funciton pointer for radius object
            real
            ( GeometryNozzle:: * mRadiusFunction )
                    ( const real          & aX )  const;

//------------------------------------------------------------------------------

            // for derivative of radius
            real
            ( GeometryNozzle:: * mSlopeFunction )
                    ( const real          & aX ) const;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // Constructor for Rao Nozzle
            GeometryNozzle(
                    const real & aThroatDiameter,
                    const real & aOpeningAngle,
                    const real & aExhaustAngle,
                    const real & aExpansionRatio,
                    const real & aCircleRadius );

            // Constructor for Bezier Nozzle
            GeometryNozzle(
                    const real & aThroatDiameter,
                    const real & aOpeningAngle,
                    const real & aExhaustAngle,
                    const real & aExpansionRatio,
                    const real & aCircleRadius,
                    const real & aLength,
                    const real & aXi,
                    const real & aEta );

            ~GeometryNozzle();

//------------------------------------------------------------------------------

            /**
             * sets the offset of this nozzle ( = combustor length )
             */
             void
             set_offset( const real & aXoff );

//------------------------------------------------------------------------------

            /**
             * returns the type of this nozzle
             */
            const NozzleType &
            type() const ;

//------------------------------------------------------------------------------

            /**
             * returns the length of this nozzle
             */
            const real &
            length() const ;

//------------------------------------------------------------------------------

            /**
             * outer radius
             */
            real
            R( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * derivative of outer radius
             */
            real
            dRdx( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * cross section
             */
            real
            A( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * derivative of cross section
             */
            real
            dAdx( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            P( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            dPdx( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            r( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            drdx( const real & aX ) const;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            compute_rao_circle();

            void
            compute_rao_coefficients();

            void
            create_bezier_object( const real & aXi, const real & aEta );

//------------------------------------------------------------------------------

            // function for small circle
            real
            radius_circle( const real & aX ) const ;

            // function for small circle
            real
            dradius_circle( const real & aX ) const ;

//------------------------------------------------------------------------------

            // funcition for rao nozzle
            real
            radius_rao( const real & aX ) const ;

            // derivative function for rao nozzle

            real
            dradius_rao( const real & aX ) const ;

//------------------------------------------------------------------------------

            // funcition for rao nozzle
            real
            radius_bezier( const real & aX ) const ;

            // derivative function for rao nozzle

            real
            dradius_bezier( const real & aX ) const ;

//------------------------------------------------------------------------------


        };

//------------------------------------------------------------------------------

        inline void
        GeometryNozzle::set_offset( const real & aXoff )
        {
            mXoff = aXoff ;
        }

//------------------------------------------------------------------------------

        inline const NozzleType &
        GeometryNozzle::type() const
        {
            return mType ;
        }

//------------------------------------------------------------------------------

        inline const real &
        GeometryNozzle::length() const
        {
            return mEx ;
        }

//------------------------------------------------------------------------------
        inline real
        GeometryNozzle::R( const real & aX ) const
        {
            if ( aX < mNx + mXoff )
            {
                return this->radius_circle( aX -mXoff );
            }
            else
            {
                return ( this->*mRadiusFunction ) ( aX-mXoff );
            }
        }

//------------------------------------------------------------------------------

        inline real
        GeometryNozzle::dRdx( const real & aX ) const
        {
            if ( aX < mNx + mXoff )
            {
                return this->dradius_circle( aX -mXoff );
            }
            else
            {
                return ( this->*mSlopeFunction ) ( aX-mXoff );
            }
        }

//------------------------------------------------------------------------------

        inline real
        GeometryNozzle::A( const real & aX ) const
        {
            return constant::pi * std::pow( this->R( aX ) , 2 );
        }

//------------------------------------------------------------------------------

        inline real
        GeometryNozzle::dAdx( const real & aX ) const
        {
            return 2.0 * constant::pi *this->R( aX ) * this->dRdx( aX );
        }

//------------------------------------------------------------------------------

        inline real
        GeometryNozzle::P( const real & aX ) const
        {
            return 2.0 * constant::pi * this->R( aX );
        }

//------------------------------------------------------------------------------

        inline real
        GeometryNozzle::dPdx( const real & aX ) const
        {
            return 2.0 * constant::pi * this->dRdx( aX );
        }

//------------------------------------------------------------------------------

        inline real
        GeometryNozzle::r( const real & aX ) const
        {
            return 0.0 ;
        }

//------------------------------------------------------------------------------

        inline real
        GeometryNozzle::drdx( const real & aX ) const
        {
            return 0.0 ;
        }

//------------------------------------------------------------------------------

    }
}

#endif //BELFEM_CL_CH_GEOMETRYNOZZLE_HPP