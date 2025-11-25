//
// Created by Christian Messe on 24.11.19.
//

#ifndef BELFEM_CL_CH_GEOMETRY_HPP
#define BELFEM_CL_CH_GEOMETRY_HPP

#include "typedefs.hpp"


namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        /**
         * The Geometry Class provides functions that specify the geometry
         * of the channel. A Geometry function can either be assigned
         * globally, or, if there are multiple blocks in the channels,
         * blockwise.
         */
        class Geometry
        {

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

           const bool mIsAxisymmetric;

           // flag telling if this function creates an inner / lower wall
           bool mHaveSecondWall;

           // width of channel, only relevant in 2D plane case
           real mWidth = 1.0;

           // length of channel, optional
           real mLength = 1.0;

           // pointers to the default functions for inner / lower radius
           // perimeter and cross section

           real
           ( Geometry:: * mFP )( const real & aX ) const;

           real
           ( Geometry:: * mFdPdx )( const real & aX ) const;

           real
           ( Geometry:: * mFp )( const real & aX ) const;

           real
           ( Geometry:: * mFdpdx)( const real & aX ) const;

           real
           ( Geometry:: * mFA )( const real & aX ) const;

           real
           ( Geometry:: * mFdAdx)( const real & aX ) const;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Geometry( const bool aIsAxisymmetric );

//------------------------------------------------------------------------------

            virtual ~Geometry() = default;

//------------------------------------------------------------------------------

            /**
             * Flag telling if this geometry is axisymmetric.
             */
            bool
            is_axisymmetric() const;

//------------------------------------------------------------------------------

            /**
             * the width of this channel. Only for axisymmetric case
             * @return
             */
            const real &
            width() const;

//------------------------------------------------------------------------------

            /**
             * return the length of this channel
             */
            const real &
            length() const;

//------------------------------------------------------------------------------

            /**
             * axisymmetric: outer radius
             * 2D plane    : upper radius
             */
            virtual real
            R( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * axisymmetric: outer perimeter
             * 2D plane    : upper perimeter
             */
            virtual real
            P( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * axisymmetric: inner radius
             * 2D plane    : lower radius
             */
            virtual real
            r( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * axisymmetric: inner perimeter
             * 2D plane    : lower perimeter
             */
            virtual real
            p( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * cross section
             */
            virtual real
            A( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * hydraulic diameter
             */
            virtual real
            Dh( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * axisymmetric: outer radius, axial derivative
             * 2D plane    : upper radius, axial derivative
             */
            virtual real
            dRdx( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * axisymmetric: outer perimeter, axial derivative
             * 2D plane    : upper perimeter, axial derivative
             */
            virtual real
            dPdx( const real & aX ) const;


//------------------------------------------------------------------------------

            /**
             * axisymmetric: inner radius, axial derivative
             * 2D plane    : lower radius, axial derivative
             */
            virtual real
            drdx( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * axisymmetric: inner perimeter, axial derivative
             * 2D plane    : lower perimeter, axial derivative
             */
            virtual real
            dpdx( const real & aX ) const;


//------------------------------------------------------------------------------

            /**
             * cross section, axial derivative
             */
            virtual real
            dAdx( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
             * hydraulic diameter derivative
             */
            real
            dDhdx( const real & aX ) const;

//------------------------------------------------------------------------------

            /**
              * Tells if this geometry has an inner / lower wall
              * should be true if not axisymmetric
              */
            bool
            has_second_wall() const;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            real
            P_axisymmetric( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            P_plane( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            p_axisymmetric( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            p_plane( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            dPdx_axisymmetric( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            dPdx_plane( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            dpdx_axisymmetric( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            dpdx_plane( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            A_axisymmetric( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            dAdx_axisymmetric( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            A_plane( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            dAdx_plane( const real & aX ) const;

//------------------------------------------------------------------------------

            void
            set_length( const real & aLength );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_CH_GEOMETRY_HPP
