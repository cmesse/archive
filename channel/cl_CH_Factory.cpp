//
// Created by Christian Messe on 21.11.20.
//
#include "constants.hpp"

#include "cl_CH_Factory.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Element_Factory.hpp"
#include "Mesh_Enums.hpp"
#include "fn_Mesh_compute_edge_lengths.hpp"
#include "CH_defines.hpp"
#include "cl_CH_Wall.hpp"
#include "cl_CH_GeometryCylinderCombustor.hpp"
#include "cl_CH_GeometryNozzle.hpp"

namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        Factory::Factory( HDF5 & aDatabase ) :
            mDatabase( aDatabase )
        {
        }

//------------------------------------------------------------------------------

        void
        Factory::create_channels( const string & aGroup,
                                  Mesh & aMesh,
                                  Cell< Segment * > & aSegments )
        {
            // select the HDF5 group
            mDatabase.select_group( aGroup );

            Matrix< id_t > tAllNodeIDs ;

            // get the coldgas nodes
            mDatabase.load_data( "ColdgasNodes", tAllNodeIDs );

            // get the number of segments
            uint tNumSegments = tAllNodeIDs.n_cols() ;

            // get the number of nodes per segment
            uint tNumNodesPerSegment = tAllNodeIDs.n_rows() ;

            // compute the reference coordinates for the surfaces
            Vector< real > tS ;
            this->compute_reference_coordinate("ChannelCenter", tS );

            // load cross sections
            Vector< real > tA ;
            mDatabase.load_data( "ChannelCrossSection", tA );

            // load perimeters
            Vector< real > tU ;
            mDatabase.load_data( "ChannelPerimeter", tU );
            //real tPerimeter ;

            // allocate segment container
            aSegments.set_size( tNumSegments, nullptr ) ;

            // temporary object for node IDs
            Vector< id_t > tNodeIDs( tNumNodesPerSegment ) ;

            // loop over all segments
            for( index_t k=0; k<tNumSegments; ++k )
            {

                // copy note IDs into temporary container
                tNodeIDs = tAllNodeIDs.col( k );

                // crate a new wall
                Wall * tWall = new Wall( aMesh, tNodeIDs );

                // we can uncomment this instead if tU( k )
                // to compute the perimeter from the mesh
                // the error is less then 1 %
                //tPerimeter = 2.0 * tWall->segment_length() ;

                // create the segment
                Segment * tSegment = new Segment( k+1, tS( k ), tA( k ), tU( k ), 1 );


                // add wall to segment
                tSegment->add_wall( 0, tWall );

                // add segment to output
                aSegments( k ) = tSegment ;
            }
        }

//------------------------------------------------------------------------------

        void
        Factory::create_cylinder_segments(
                Geometry          * aGeometry,
                Mesh              & aMesh,
                Cell< Segment * > & aSegments,
                bool                aReverse )
        {
            // select the HDF5 group
            mDatabase.select_group( "Chamber" );

            // get the type
            uint tType;
            mDatabase.load_data("Type", tType );

            BELFEM_ERROR( tType == 0, "Type must be cylinder chamber" );

            uint tNumElems;
            mDatabase.load_data( "NumElems", tNumElems );
            uint tNumSegments = 2 * tNumElems + 1;
            mDatabase.close_active_group() ;

            // allocate the segment container
            aSegments.set_size( tNumSegments, nullptr );

            mDatabase.select_group("Liner" );

            Matrix< id_t > tAllNodeIDs;
            mDatabase.load_data("HotgasNodes", tAllNodeIDs );

            Vector< id_t > tNodeIDs( tAllNodeIDs.n_rows() );

            real tLength = aGeometry->length() ;

            for( index_t k=0; k<tNumSegments; ++k )
            {
                // grab the first node from the mesh
                mesh::Node * tNode = aMesh.node( tAllNodeIDs( 0, k ) );

                // grab the X-Position
                real tX = tNode->x();

                // compute the cross section
                //real tA = constant::pi * tNode->y() * tNode->y() ;
                real tA = aGeometry->A( tX );

                // compute the perimeter
                //real tU = 2.0 * constant::pi * tNode->y() ;
                real tU = aGeometry->P( tX );

                if( aReverse )
                {
                    tX = tLength - tNode->x();
                }

                // create the segment
                Segment * tSegment = new Segment( k+1,
                                                  tX,
                                                  tA,
                                                  tU,
                                                  1 );

                // copy note IDs into temporary container
                tNodeIDs = tAllNodeIDs.col( k );

                // crate a new wall
                Wall * tWall = new Wall( aMesh, tNodeIDs );

                // add wall to segment
                tSegment->add_wall( 0, tWall );

                // add segment to output
                if( aReverse )
                {
                    aSegments( tNumSegments - k - 1 ) = tSegment;
                }
                else
                {
                    aSegments( k ) = tSegment ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Factory::create_nozzle_segments(
                Geometry          * aGeometry,
                Mesh              & aMesh,
                Cell< Segment * > & aSegments )
        {
            // select the HDF5 group
            mDatabase.select_group( "Nozzle" );

            // get the number of elements for the nozzle
            uint tNumElements ;
            mDatabase.load_data("NumElems", tNumElements );

            // cumpute the number of nodes
            index_t tNumNodes = 2 * tNumElements + 1 ;

            // close this group
            mDatabase.close_active_group() ;

            // select the liner group
            mDatabase.select_group( "Liner" );

            // get the node IDs
            Matrix < id_t > tAllNodeIDs ;
            mDatabase.load_data( "HotgasNodes", tAllNodeIDs );

            // total number of nodes
            index_t tNumAllNodes = tAllNodeIDs.n_cols() ;

            // node IDs per segment
            Vector< id_t > tNodeIDs( tAllNodeIDs.n_rows() );

            // id of first node
            index_t tNode0 = tNumAllNodes - tNumNodes ;

            // the output cell
            aSegments.set_size( tNumNodes, nullptr );

            // counter for segments
            index_t tCount = 0 ;


            for( index_t k=tNode0; k<tNumAllNodes; ++k )
            {
                // grab the first node from the mesh
                mesh::Node * tNode = aMesh.node( tAllNodeIDs( 0, k ) );

                // grab the X-Position
                real tX = tNode->x();

                // compute the cross section
                //real tA = constant::pi * tNode->y() * tNode->y() ;
                real tA = aGeometry->A( tX );

                // compute the perimeter
                //real tU = 2.0 * constant::pi * tNode->y() ;
                real tU = aGeometry->P( tX );

                // create the segment
                Segment * tSegment = new Segment( k+1,
                                                  tX,
                                                  tA,
                                                  tU,
                                                  1 );

                // copy note IDs into temporary container
                tNodeIDs = tAllNodeIDs.col( k );

                // crate a new wall
                Wall * tWall = new Wall( aMesh, tNodeIDs );

                // add wall to segment
                tSegment->add_wall( 0, tWall );

                // add segment to output cell
                aSegments( tCount++ ) = tSegment ;

            }
        }

//------------------------------------------------------------------------------

        channel::Geometry *
        Factory::create_cylinder_geometry()
        {
            // select group in database
            mDatabase.select_group( "Chamber" );

            // read geometry data
            real tChamberDiameter ;
            mDatabase.load_data( "ChamberDiameter", tChamberDiameter );

            real tChamberLength ;
            mDatabase.load_data( "ChamberLength", tChamberLength );

            real tCurvatureRadius ;
            mDatabase.load_data( "CurvatureRadius", tCurvatureRadius );

            real tCylinderLength ;
            mDatabase.load_data( "CylinderLength", tCylinderLength );

            real tKinkRadius ;
            mDatabase.load_data( "KinkRadius", tKinkRadius );

            real tThroatDiameter ;
            mDatabase.load_data( "ThroatDiameter", tThroatDiameter );

            // close this group
            mDatabase.close_active_group() ;

            // return the new object
            return new GeometryCylinderCombustor(
                    tThroatDiameter,
                    tChamberDiameter,
                    tCylinderLength,
                    tChamberLength,
                    tKinkRadius,
                    tCurvatureRadius );
        }

//------------------------------------------------------------------------------

        channel::Geometry *
        Factory::create_nozzle_geometry()
        {
            // the throat diameter is in another Groip, we must get that first
            mDatabase.select_group("Chamber");
            real tThroatDiameter ;
            mDatabase.load_data("ThroatDiameter", tThroatDiameter );

            // needed for the offset
            real tXoff ;
            mDatabase.load_data( "ChamberLength", tXoff );

            mDatabase.close_active_group() ;

            // select the nozzle group in database
            mDatabase.select_group( "Nozzle" );

            // the radius of the small circle at the beginning
            real tCircleRadius ;
            mDatabase.load_data( "CircleRadius", tCircleRadius );

            // the opening angle of the nozzle near the small circle,
            // need to covert to rad
            real tOpeningAngle ;
            mDatabase.load_data("OpeningAngle", tOpeningAngle );
            tOpeningAngle *= constant::deg ;

            // the angle at the end of the nozzle, need to convert to rad
            real tExhaustAngle ;
            mDatabase.load_data( "ExhaustAngle", tExhaustAngle );
            tExhaustAngle *= constant::deg ;

            // the ratio between the cross section at the exhaust and the throat
            real tExpansionRatio ;
            mDatabase.load_data("ExpansionRatio", tExpansionRatio );

            // the type of the nozzle. 0: Rao, 1: Bezier
            uint tType ;
            mDatabase.load_data("Type", tType );

            // these values depend on the type :

            // the length of the nozzle ( Bezier only )
            real tLength = BELFEM_QUIET_NAN ;

            // help parameters for the bezier curve
            real tXi = BELFEM_QUIET_NAN ;
            real tEta = BELFEM_QUIET_NAN ;

            if ( tType == 1 )
            {
                mDatabase.load_data( "Length", tLength );
                mDatabase.load_data( "xi", tXi );
                mDatabase.load_data( "eta", tEta );
            }

            // close this group
            mDatabase.close_active_group() ;

            channel::GeometryNozzle * aNozzle = nullptr ;

            // return the nozzle
            if( tType == 0 )
            {
                // create a new rao nozzle
                aNozzle = new GeometryNozzle(
                        tThroatDiameter,
                        tOpeningAngle,
                        tExhaustAngle,
                        tExpansionRatio,
                        tCircleRadius );
            }
            else if ( tType == 1 )
            {
                // create a new Bezier nozzle
                aNozzle = new GeometryNozzle(
                        tThroatDiameter,
                        tOpeningAngle,
                        tExhaustAngle,
                        tExpansionRatio,
                        tCircleRadius,
                        tLength,
                        tXi,
                        tEta );
            }
            else
            {
                BELFEM_ERROR( false, "Unknown Nozzle Type" );
            }

            // set the offset
            aNozzle->set_offset( tXoff );

            return aNozzle ;
        }

//------------------------------------------------------------------------------


        void
        Factory::compute_reference_coordinate(
                const string & aMatrixLabel,
                Vector< real > & aS )
        {
            // the reference coordinates
            Matrix< real > tX ;
            mDatabase.load_data( aMatrixLabel, tX ) ;

            uint tNumNodes = tX.n_cols() ;
            uint tNumElems = ( tNumNodes-1 ) /  2 ;

            // create some nodes
            Cell< mesh::Node * > tNodes( tNumNodes, nullptr ) ;

            for ( index_t k=0; k<tNumNodes; ++k )
            {
                tNodes( k ) = new mesh::Node( k+1, tX(0,k), tX( 1, k ) );
            }

            // create an element container
            Cell< mesh::Element * > tElements( tNumElems, nullptr );

            // create an element factory
            mesh::ElementFactory tFactory ;

            // reference offset
            uint tOff = 0 ;

            // create and link elements
            for ( index_t e=0; e<tNumElems; ++e )
            {
                mesh::Element * tElement = tFactory.create_element( ElementType::LINE3, e+1 ) ;

                tElement->insert_node( tNodes( tOff ), 0 ) ;
                tElement->insert_node( tNodes( tOff+2 ), 1 ) ;
                tElement->insert_node( tNodes( tOff+1 ), 2 ) ;
                tOff += 2 ;

                // add element to container
                tElements( e ) = tElement ;
            }

            Vector< real > tEdgeLengths ;
            mesh::compute_edge_lengths( 2, tElements, tEdgeLengths ) ;

            // from the edge lengths, we compute the surface coordinates

            // reset offset
            tOff = 0 ;

            aS.set_size( tNumNodes );

            // first entry
            aS( 0 ) = 0.0 ;

            for( index_t e=0; e<tNumElems; ++e )
            {
                aS( tOff+1 ) = aS( tOff ) + 0.5 * tEdgeLengths( e ) ;
                aS( tOff+2 ) = aS( tOff ) + tEdgeLengths( e ) ;
                tOff += 2 ;
            }

            // tidy up
            for ( mesh::Element * tElement : tElements )
            {
                delete tElement ;
            }
            for ( mesh::Node * tNode : tNodes )
            {
                delete tNode ;
            }
        }
//------------------------------------------------------------------------------

    }
}
