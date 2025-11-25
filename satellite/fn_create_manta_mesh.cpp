//
// Created by christian on 9/14/24.
//

#include "fn_create_manta_mesh.hpp"
#include "commtools.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_HDF5.hpp"
#include "fn_sum.hpp"
#include "fn_unique.hpp"
#include "fn_max.hpp"

namespace belfem
{
    Mesh *
    create_manta_mesh( const string & aPath )
    {
        Mesh * aMesh = nullptr ;

        if( comm_rank() == 0 )
        {
            HDF5 tFile( aPath, FileMode::OPEN_RDONLY );

            tFile.select_group("keyframe_00");
            Matrix< real > tX ;
            tFile.load_data("nodeCoords", tX );
            tFile.close_active_group();

            tFile.select_group("Mesh");

            // ids for the nodes
            Vector<  unsigned int > tNodeIDs ;
            tFile.load_data( "NodeIDs", tNodeIDs );

            Vector<  unsigned int > tGeometry ;
            tFile.load_data( "BlockIDs", tGeometry );

            Vector<  unsigned int > tBlockIDs( tGeometry.vector_data() );

            unique( tBlockIDs );

            // create the block map
            Map< index_t, id_t > tBlockMap ;
            index_t c = 0 ;
            for( id_t b: tBlockIDs )
            {
                tBlockMap[ b ] = c++ ;
            }

            // load node adjacency
            Matrix< unsigned int > tAdjacency;
            tFile.load_data("ElementTopology", tAdjacency);

            Vector< unsigned int > tElementIDs;
			tFile.load_data("ElementIDs", tElementIDs);

            aMesh = new Mesh( 3, 0, true );

            // now we can create the nodes
            index_t tNumNodes = tNodeIDs.length() ;

            Cell< mesh::Node * > & tNodes = aMesh->nodes() ;
            tNodes.set_size( tNumNodes, nullptr );

            Map< id_t, mesh::Node * > tNodeMap ;

            for( index_t k = 0; k<tNumNodes; ++k )
            {
                tNodes( k ) = new mesh::Node( tNodeIDs( k ), tX( k, 0 ), tX( k, 1 ), tX( k, 2 ) );
                tNodeMap[ tNodeIDs( k )] = tNodes( k ) ;
            }

            id_t tNumBlocks = tBlockIDs.length() ;

            // count elements per block
            Vector< id_t > tCount( tNumBlocks, 0 );
            for( luint b : tGeometry )
            {
                ++tCount( tBlockMap( b ) );
            }


            // create the blocks
            Cell< mesh::Block * > & tBlocks = aMesh->blocks() ;
            tBlocks.set_size( tNumBlocks, nullptr );
            c = 0 ;
            for( id_t b: tBlockIDs )
            {
                tBlocks( c ) = new mesh::Block( b, tCount( tBlockMap( b ) ) );
                ++c ;
            }

            // create the elements
            index_t tNumElems = sum( tCount );

            mesh::ElementFactory tFactory ;

            //Cell< mesh::Element * > tElements( tNumElems, nullptr );


            for(index_t e=0; e<tNumElems; ++e )
            {
                mesh::Element * tElement = tFactory.create_element( ElementType::TRI3, tElementIDs( e ) );

                for( uint k=0; k<3; ++k )
                {
                    tElement->insert_node( tNodes( tAdjacency( e, k ) ), k );
                }
                tElement->set_geometry_tag( tGeometry( e ) );
                tBlocks( tBlockMap( tGeometry( e ) ) )->insert_element( tElement );

                //tElements( e ) = tElement ;
            }

            aMesh->finalize() ;

            /* reorganize the elements in their original order
            c = 0 ;
            for( mesh::Element * tElement : tElements )
            {
              	tElement->set_index( c++ );
            }

            std::memmove(  aMesh->elements().data(), tElements.data(), tNumElems * sizeof( mesh::Element * ) )*/

            tFile.close_active_group() ;
            tFile.close();
        }
        else
        {
            aMesh = new Mesh( 3, 0, true ) ;

            comm_barrier() ;
        }


        return aMesh ;
    }

}
