//
// Created by Christian Messe on 23.11.20.
//

#ifndef BELFEM_CL_CH_WALL_HPP
#define BELFEM_CL_CH_WALL_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_Vector.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"

namespace belfem
{
    namespace channel
    {
        /**
         * The channel wall has two functions
         *
         * 1. Compute the average surface temperature
         * 2. Impose Boundary Condition of Mesh
         */
        class Wall
        {
            Mesh & mMesh ;

            // number of nodes on this wall
            const index_t mNumNodes ;

            // number of elements on this wall
            const index_t mNumElements ;

            // Cell of nodes on this mesh, owned and destroyed by mesh
            Cell< mesh::Node * >  mNodes ;

            // Cell of elements, owned and destroyed by wall
            Cell< mesh::Element * > mElements ;

            // container with length of each element
            Vector< real > mElementLengths ;

            // length of segment along wall ( should be half perimeter )
            real mSegmentLength ;

            // Shape Function for integration
            Vector< real > mIntegrationWeights ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Wall( Mesh & aMesh, Vector< id_t > & aNodeIDs );

//------------------------------------------------------------------------------

            ~Wall();

//------------------------------------------------------------------------------

            /**
             * grab the surface temperature from the mesh and average
             * along segment
             * @return
             */
            real
            average_surface_temperature();

//------------------------------------------------------------------------------

            /**
             * grab the surface temperature from the mesh and average
             * along segment
             * @return
             */
            real
            average_heatload( const real & aAlpha, const real & aTinf );

//------------------------------------------------------------------------------

            /**
             * impose pressures the field on the mesh
             * careful: fields must exist
             */
            void
            set_flowdata( const real & aT, const real & aP, const real & aMa ) ;

//------------------------------------------------------------------------------

            /**
             * return the lenth of the segment
             */
             const real &
             segment_length() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            collect_nodes_from_mesh( const Vector< id_t > & aNodeIDs );

//------------------------------------------------------------------------------

            void
            create_integration_elements();

//------------------------------------------------------------------------------

            void
            initialize_integration_weights();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline const real &
        Wall::segment_length() const
        {
            return mSegmentLength ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_CH_WALL_HPP
