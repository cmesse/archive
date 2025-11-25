//
// Created by Christian Messe on 02.12.20.
//

#ifndef BELFEM_CL_CH_GEOMETRYCYLINDERCOMBUSTOR_HPP
#define BELFEM_CL_CH_GEOMETRYCYLINDERCOMBUSTOR_HPP

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_CH_Geometry.hpp"

namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        class GeometryCylinderCombustor : public Geometry
        {
            // chamber radius
            real mR ;

            // kink radius
            real mRk ;

            // curvature radius
            real mRc ;

            // slope between P and Q
            real mA ;

            // offset for line between P and Q
            real mB ;

            // center for kink circle
            real mKx;
            real mKr ;

            // help point, end of kink circle
            real mPx ;
            real mPr ;

            // help point, beginning of nozzle circle
            real mQx ;
            real mQr ;

            // center for nozzle circle
            real mMx;
            real mMr ;


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            GeometryCylinderCombustor(
                    const real & aThroatDiameter,
                    const real & aChamberDiameter,
                    const real & aCylinderLength,
                    const real & aChamberLength,
                    const real & aKinkRadius,
                    const real & aCurvatureRadius );

            ~GeometryCylinderCombustor() = default ;

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
        };

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_CH_GEOMETRYCYLINDERCOMBUSTOR_HPP
