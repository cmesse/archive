//
// Created by Christian Messe on 24.11.19.
//

#ifndef BELFEM_CL_CH_GEOMETRYCOMBUSTOR_HPP
#define BELFEM_CL_CH_GEOMETRYCOMBUSTOR_HPP

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_CH_Geometry.hpp"

namespace belfem
{
    namespace channel
    {
        /**
         * This is am example class that creates the geometry of the ITLR combustor,
         * as seen in DOI: 10.18419/opus-9381, Fig. 6.1
         */
        class GeometryCombustor : public Geometry
        {
            // initial height
            real mHeight0 = 0.04;

            // width
            real mWidth = 0.065;

            // length until opening
            real mLength0 = 0.03;

            // flat length
            real mLength1 = 0.9;

            // corner radius in perpendicular plane
            real mCornerRadius = 0.01;

            // bending radius
            real mKinkRadius = 0.1;

            // opening angle in rad
            real mAngle = 1.0 * constant::deg;

            // position of injector
            real mInjectorPosition = 0.116;

            // Length of injector
            real mInjectorLength = 0.086;

            // height of injector
            real mInjectorHeight = 0.007;

            // Geometry points
            Matrix< real > mPoints;

            // Coefficients for Lines
            Matrix< real > mLines;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

            GeometryCombustor();

//------------------------------------------------------------------------------

            virtual ~GeometryCombustor() = default;

//------------------------------------------------------------------------------

            real
            R( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            dRdx( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            A( const real & aX ) const;

//------------------------------------------------------------------------------

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

            const Matrix< real > &
            points() const;

//------------------------------------------------------------------------------

            real
            injector_position() const;

//------------------------------------------------------------------------------
//        private:
//------------------------------------------------------------------------------

            real
            r_injector( const real & aX ) const;

//------------------------------------------------------------------------------

            real
            drdx_injector( const real & aX ) const;

//------------------------------------------------------------------------------


            /**
             * initialize relevant constants for the function
             */
            void
            initialize();

//------------------------------------------------------------------------------

        };
    }
}
#endif //BELFEM_CL_CH_GEOMETRYCOMBUSTOR_HPP
