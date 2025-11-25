//
// Created by Christian Messe on 21.11.20.
//

#ifndef BELFEM_CL_CH_FACTORY_HPP
#define BELFEM_CL_CH_FACTORY_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_HDF5.hpp"

#include "CH_Enums.hpp"
#include "cl_CH_Wall.hpp"
#include "cl_CH_Segment.hpp"
#include "cl_CH_Geometry.hpp"

namespace belfem
{
    namespace channel
    {
        class Factory
        {
            HDF5 & mDatabase ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Factory( HDF5 & aDatabase );

//------------------------------------------------------------------------------

            ~Factory() = default ;

//------------------------------------------------------------------------------

            void
            create_channels(
                    const string      & aGroup,
                    Mesh              & aMesh,
                    Cell< Segment * > & aSegments ) ;

//------------------------------------------------------------------------------

            void
            create_cylinder_segments(
                    Geometry          * aGeometry,
                    Mesh              & aMesh,
                    Cell< Segment * > & aSegments,
                    bool                aReverse );

//------------------------------------------------------------------------------

            void
            create_nozzle_segments(
                    Geometry          * aGeometry,
                    Mesh              & aMesh,
                    Cell< Segment * > & aSegments );

//------------------------------------------------------------------------------

            channel::Geometry *
            create_cylinder_geometry();

//------------------------------------------------------------------------------

            channel::Geometry *
            create_nozzle_geometry();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /**
             * this subroutine takes the X and R data and computes a local
             * reference coordinate
             */
             void
             compute_reference_coordinate(
                     const string & aMatrixLabel,
                     Vector< real > & aS );

//------------------------------------------------------------------------------

        };
    }
}
#endif //BELFEM_CL_CH_FACTORY_HPP
