
//
// Created by Christian Messe on 22.01.21.
//

#ifndef BELFEM_CL_EN_STEPANOFFCHART_HPP
#define BELFEM_CL_EN_STEPANOFFCHART_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace engine
    {
        class StepanoffChart
        {
            // lower polynomial for Km1
            Vector< real > mKm1Poly0 ;

            // lower polynomial for Km2
            Vector< real > mKm2Poly0 ;

            // transitioning poly for Km1
            Vector< real > mKm1Poly1 ;

            // transitioning poly for Km2
            Vector< real > mKm2Poly1 ;

            // common poly for Km1 and Km2
            Vector< real > mKmPoly2 ;

            const real mSwitch0 = 8.2 ;
            const real mSwitch1 = 9.2 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            StepanoffChart() ;

//------------------------------------------------------------------------------

            ~StepanoffChart() = default ;

//------------------------------------------------------------------------------

            real
            Km1( const real & aNs ) const ;

//------------------------------------------------------------------------------

            real
            Km2( const real & aNs ) const ;

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_EN_STEPANOFFCHART_HPP
