//
// Created by Christian Messe on 22.01.21.
//
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_create_beam_poly.hpp"
#include "cl_EN_StepanoffChart.hpp"

namespace belfem
{
    namespace engine
    {
        StepanoffChart::StepanoffChart()
        {
            // lower polynomial for Km1
            mKm1Poly0 = { 2.348135800197E-02, + 1.126407158777E-01, -3.880095474057E+00 };

            // lower polynomial for Km2
            mKm2Poly0 = { 9.819000157735E-02, -9.356590933327E-01,  -5.309983165966E-01 };

            // upper polynomial for Km1 and Km2
            mKmPoly2 = { 2.332232691271E-01, - 3.512801133446E+00, 1.178348572838E+01};

            // create transitioning Polynomial for Km1
            create_beam_poly( mSwitch0,
                              polyval( mKm1Poly0, mSwitch0 ),
                              dpolyval( mKm1Poly0, mSwitch0),
                              mSwitch1,
                              polyval( mKmPoly2, mSwitch1 ),
                              polyval( mKmPoly2, mSwitch1 ),
                              mKm1Poly1 );

            // create transitioning Polynomial for Km2
            create_beam_poly( mSwitch0,
                              polyval( mKm2Poly0, mSwitch0 ),
                              dpolyval( mKm2Poly0, mSwitch0),
                              mSwitch1,
                              polyval( mKmPoly2, mSwitch1 ),
                              polyval( mKmPoly2, mSwitch1 ),
                              mKm2Poly1 );

        }
//------------------------------------------------------------------------------

        real
        StepanoffChart::Km1( const real & aNs ) const
        {
            real tX = std::log( aNs );

            if( tX < mSwitch0 )
            {
                return std::exp( polyval( mKm1Poly0, tX ) );
            }
            else if ( tX < mSwitch1 )
            {
                return std::exp( polyval( mKm1Poly1, tX ) );
            }
            else
            {
                return std::exp( polyval( mKmPoly2, tX ) );
            }
        }

//------------------------------------------------------------------------------

        real
        StepanoffChart::Km2( const real & aNs ) const
        {
            real tX = std::log( aNs );

            if( tX < mSwitch0 )
            {
                return std::exp( polyval( mKm2Poly0, tX ) );
            }
            else if ( tX < mSwitch1 )
            {
                return std::exp( polyval( mKm2Poly1, tX ) );
            }
            else
            {
                return std::exp( polyval( mKmPoly2, tX ) );
            }
        }

//------------------------------------------------------------------------------
    }
}