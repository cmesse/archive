//
// Created by Christian Messe on 21.11.20.
//

#ifndef BELFEM_CL_CH_SEGMENT_HPP
#define BELFEM_CL_CH_SEGMENT_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_Vector.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"
#include "cl_CH_Wall.hpp"

namespace belfem
{
    namespace channel
    {
        class Factory ;

        class Segment
        {
            // id for this segment
            const id_t mID ;

            // number of connected walls
            const uint mNumWalls ;

            Vector< real > mData ;

            Cell< Wall * > mWalls ;

            friend Factory ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Segment( const id_t & aID,
                     const real & aX,
                     const real & aA,
                     const real & aU,
                     const uint & aNumWalls ) ;

//------------------------------------------------------------------------------

            ~Segment();

//------------------------------------------------------------------------------

            void
            print();

//------------------------------------------------------------------------------

            /**
             * return the x-position of this segment
             */
            const real &
            x() const ;


//------------------------------------------------------------------------------

            /**
             * return the cross section of this segment
             */
            const real &
            cross_section() const ;

//------------------------------------------------------------------------------

            /**
             * return the perimeter of this segment
             */
            const real
            perimeter() const ;

//------------------------------------------------------------------------------

            /**
             * return the hydraulic diameter of this segment
             */
            const real &
            hydraulic_diameter() const ;

//------------------------------------------------------------------------------

            /**
             * expose the data container
             */
             Vector< real > &
             data() ;

//------------------------------------------------------------------------------

            /**
             * expose one single value
             */
            const real &
            value( const index_t & aIndex ) const ;

//------------------------------------------------------------------------------

            /**
             * overwrite one specific value
             */
             void
             set_value( const index_t & aIndex, const real & aValue );

//------------------------------------------------------------------------------

            /**
             * get wall temperatures from mesh
             */
             void
             pull_surface_temperatures();

//------------------------------------------------------------------------------

            /**
             * send the heatloads to the mesh
             */
            void
            push_heatloads();

//------------------------------------------------------------------------------

            /**
             * send the pressures to the mesh
             */
            void
            push_flowdata();

//------------------------------------------------------------------------------

            /**
             * return the number of walls
             */
             const uint &
             num_walls() const;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            void
            add_wall( const uint & aWallIndex, Wall * aWall );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline const real &
        Segment::x() const
        {
            return mData( 0 );
        }
//------------------------------------------------------------------------------
        inline const real &
        Segment::cross_section() const
        {
            return mData( 1 );
        }

//------------------------------------------------------------------------------

        inline const real
        Segment::perimeter() const
        {
            return 4.0 * mData( 1 ) / mData( 2 );
        }

//------------------------------------------------------------------------------

        inline const real &
        Segment::hydraulic_diameter() const
        {
            return mData( 2 );
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Segment::data()
        {
            return mData ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Segment::value( const index_t & aIndex ) const
        {
            return mData( aIndex );
        }

//------------------------------------------------------------------------------

        inline void
        Segment::set_value( const index_t & aIndex, const real & aValue )
        {
            mData( aIndex ) = aValue ;
        }

//------------------------------------------------------------------------------

        inline const uint &
        Segment::num_walls() const
        {
            return mNumWalls ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_CH_SEGMENT_HPP
