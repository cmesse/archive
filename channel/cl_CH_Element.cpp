//
// Created by Christian Messe on 27.11.20.
//

#include "cl_CH_Element.hpp"

namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        Element::Element( Segment & aSegment0,
                          Segment & aSegment1,
                          Segment & aSegment2 ) :
            mSegment0( aSegment0 ),
            mSegment1( aSegment1 ),
            mSegment2( aSegment2 ),
            mLength( aSegment1.x() - aSegment0.x() )
        {

        }

//------------------------------------------------------------------------------

        void
        Element::collect_data( const uint & aIndex, Vector< real > & aValues ) const
        {
            aValues( 0 ) = mSegment0.value( aIndex );
            aValues( 1 ) = mSegment1.value( aIndex );
            aValues( 2 ) = mSegment2.value( aIndex );
        }

//------------------------------------------------------------------------------

        void
        Element::compute_N( const real & aX , Vector< real > & aN ) const
        {
            real tXi = ( aX - this->x0() ) / mLength ;

            aN( 0 ) = tXi * ( 2.0 * tXi - 3.0 ) + 1.0 ;
            aN( 1 ) = tXi * ( 2.0 * tXi - 1.0 );
            aN( 2 ) = 4.0 * tXi * ( 1.0 - tXi );

            /*aN( 0 ) = 1.0 - tXi ;
            aN( 1 ) = tXi ;
            aN( 2 ) = 0.0 ; */
        }

//------------------------------------------------------------------------------

        void
        Element::compute_B( const real & aX , Vector< real > & aB ) const
        {
            real tXi = ( aX - this->x0() ) / mLength ;

            aB( 0 ) = 4.0 * tXi - 3.0 ;
            aB( 1 ) = 4.0 * tXi - 1.0 ;
            aB( 2 ) = 4.0 - 8.0 * tXi ;

            /*aB( 0 ) = -1.0 ;
            aB( 1 ) =  1.0 ;
            aB( 2 ) =  0.0 ; */

            aB /= mLength ;
        }

//------------------------------------------------------------------------------
    }
}