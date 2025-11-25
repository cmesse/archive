//
// Created by Christian Messe on 13.06.20.
//

#ifndef BELFEM_CL_BL_STREAMLINE_HPP
#define BELFEM_CL_BL_STREAMLINE_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Element.hpp"
#include "cl_Mesh.hpp"
#include "cl_Gas.hpp"
#include "cl_BL_Panel.hpp"
#include "cl_BL_StagnationPoint.hpp"

namespace belfem
{
    namespace boundarylayer
    {
//----------------------------------------------------------------------------

        class StreamLine
        {

            StagnationPoint & mStagnationPoint ;
            Gas   & mGas ;
            State & mFreestream ;
            State & mStagnation ;

            Mesh *     mMesh ;
            const id_t mID ;
            const ElementType mElementType ;

            // all panels
            Cell< Panel * > mPanels ;

            // lower and upper panels
            Cell< Panel * > mLowerPanels ;
            Cell< Panel * > mUpperPanels ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            StreamLine(
                    StagnationPoint & aStagnationPoint,
                    Mesh * aMesh,
                    const id_t aID,
                    Cell< mesh::Element * > & aElements );


//----------------------------------------------------------------------------

            ~StreamLine();

//----------------------------------------------------------------------------

            // expose the mesh object
            Mesh *
            mesh();

//----------------------------------------------------------------------------

            // expose the gas object
            Gas &
            gas();

//----------------------------------------------------------------------------

            /**
             *
             * @param aAoA angle of attack in deg
             */
            void
            compute( const real & aAoA );

//----------------------------------------------------------------------------

            /*
             * return the freestream state
             */
            State &
            freestream();

//------------------------------------------------------------------------------

            /*
             * return the stagnation state
             */
            State &
            stagnation();

//----------------------------------------------------------------------------

            void
            print();

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------
// FUNCTIONS NEEDED FOR CONSTRUCTION
//----------------------------------------------------------------------------
            // make sure that the passed elements make sens
            void
            check_element_sanity( Cell< mesh::Element * > & aElements );

//----------------------------------------------------------------------------
            // make sure that the passed elements are continuous
            bool
            check_element_continuity( Cell< mesh::Element * > & aElements );

//----------------------------------------------------------------------------

            // make sure that the passed elements are all of the same type
            bool
            check_element_types( Cell< mesh::Element * > & aElements );

//----------------------------------------------------------------------------

            // create the panel objects for this streamline
            void
            create_panels( Cell< mesh::Element * > & aElements );

//----------------------------------------------------------------------------

            // compute data needed for the panels
            void
            compute_direction_vectors(
                    Cell< mesh::Element * > & aElements,
                    Matrix< real >          & aDirections );

//----------------------------------------------------------------------------

            // get the parameter coordinates for the element
            void
            populate_shape_derivative( Cell< Matrix< real > >  & adNdXi );

//----------------------------------------------------------------------------

            void
            collect_nodes(
                    Cell< mesh::Element * > & aElements,
                    Cell< mesh::Node * >    & aNodes ) ;

//----------------------------------------------------------------------------

            void
            collect_node_normals(
                    Cell< mesh::Node * >    & aNodes ,
                    Matrix< real >          & aNormals );

//----------------------------------------------------------------------------

            void
            compute_surface_coordinates(
                    Cell< mesh::Element * > & aElements,
                    Vector< real >          & aSurfaceCoordinates );

//----------------------------------------------------------------------------
// FUNCTIONS NEEDED FOR SURFACE INCLINATION METHOD
//----------------------------------------------------------------------------

            void
            compute_freestream_direction(
                    const real & aAoA,
                    Vector< real > & aFreestreamDirection );

//----------------------------------------------------------------------------

            void
            compute_modified_newton( const Vector< real > & aFreestreamDirection );

//----------------------------------------------------------------------------

            index_t
            find_stagnation_point();

//----------------------------------------------------------------------------

            real
            compute_surface_coordinate(
                    Vector< real > & aFreestreamDirection,
                    const index_t    aStagnationIndex );

//----------------------------------------------------------------------------

            // split the streamline into a lower and an upper part
            index_t
            split_streamline( Vector< real > & aFreestreamDirection );

//----------------------------------------------------------------------------

            // compute the prandtl meyer expansion.
            // argument is either mLowerPanels or mUpperPanels
            void
            compute_prandtl_meyer( Cell< Panel * > & aPanels ) ;

//----------------------------------------------------------------------------

            void
            print( Cell< Panel * > & aPanels  );

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------

        inline Mesh *
        StreamLine::mesh()
        {
            return mMesh ;
        }

//----------------------------------------------------------------------------

        inline Gas &
        StreamLine::gas()
        {
            return mGas ;
        }

//----------------------------------------------------------------------------

        /*
         * return the freestream state
         */
        inline State &
        StreamLine::freestream()
        {
            return mFreestream ;
        }

//------------------------------------------------------------------------------

        /*
         * return the stagnation state
         */
        inline State &
        StreamLine::stagnation()
        {
            return mStagnation ;
        }

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_BL_STREAMLINE_HPP
