//
// Created by Christian Messe on 27.11.20.
//

#ifndef BELFEM_CL_CH_CELL_HPP
#define BELFEM_CL_CH_CELL_HPP
#include "cl_CH_Segment.hpp"

namespace belfem
{
    namespace channel
    {
        /**
         * an element contains of three segments
         */
        class Element
        {
            // the entry segment
            Segment & mSegment0 ;

            // the middle segment
            Segment & mSegment1 ;

            // the exit segment
            Segment & mSegment2 ;

            const real mLength ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Element( Segment & aSegment0,
                     Segment & aSegment1,
                     Segment & aSegment2 );

//------------------------------------------------------------------------------

            ~Element() = default ;

//------------------------------------------------------------------------------

            /**
             * return the length of this element
             */
            const real &
            length() const ;

//------------------------------------------------------------------------------

            /**
             * coordinate of entry point
             * @return
             */
            const real &
            x0() const ;

//------------------------------------------------------------------------------

            /**
              * coordinate of exit point
              * @return
              */
            const real &
            x1() const ;

//------------------------------------------------------------------------------

             /**
              * coordinate of middle point
              * @return
              */
            const real &
            x2() const ;

//------------------------------------------------------------------------------

            void
            compute_N( const real & aX , Vector< real > & aN ) const ;

//------------------------------------------------------------------------------

            void
            compute_B( const real & aX , Vector< real > & aN ) const;

//------------------------------------------------------------------------------

            void
            collect_data( const uint & aIndex, Vector< real > & aValues ) const;

//------------------------------------------------------------------------------

            // return the segment at the entry
            Segment *
            segment0();

//------------------------------------------------------------------------------

            // return the segment at the exit
            Segment *
            segment1();

//------------------------------------------------------------------------------

            // return the segment in the center
            Segment *
            segment2();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const real &
        Element::length() const
        {
            return mLength ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Element::x0() const
        {
            return mSegment0.value( 0 );
        }

//------------------------------------------------------------------------------

        inline const real &
        Element::x1() const
        {
            return mSegment1.value( 0 );
        }

//------------------------------------------------------------------------------

        inline const real &
        Element::x2() const
        {
            return mSegment2.value( 0 );
        }

//------------------------------------------------------------------------------

        // return the segment at the entry
        inline Segment *
        Element::segment0()
        {
            return & mSegment0 ;
        }

//------------------------------------------------------------------------------

        // return the segment at the exit
        inline Segment *
        Element::segment1()
        {
            return & mSegment1 ;
        }

//------------------------------------------------------------------------------

        // return the segment in the center
        inline Segment *
        Element::segment2()
        {
            return & mSegment2 ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_CH_CELL_HPP
