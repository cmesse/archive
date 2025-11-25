//
// Created by Christian Messe on 23.11.20.
//
#include "constants.hpp"
#include "cl_CH_Wall.hpp"
#include "cl_Element_Factory.hpp"
#include "Mesh_Enums.hpp"
#include "fn_Mesh_compute_edge_lengths.hpp"
#include "fn_sum.hpp"
#include "fn_intpoints.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "fn_dot.hpp"
#include "assert.hpp"
namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        Wall::Wall( Mesh & aMesh, Vector< id_t > & aNodeIDs ) :
                mMesh( aMesh ),
                mNumNodes( aNodeIDs.length() ),
                mNumElements( ( aNodeIDs.length() - 1 ) / 2 )
        {
            this->collect_nodes_from_mesh( aNodeIDs );
            this->create_integration_elements();
            this->initialize_integration_weights();
        }


//------------------------------------------------------------------------------

        Wall::~Wall()
        {
            // delete the elements
            for( mesh::Element * tElement : mElements )
            {
                delete tElement;
            }
        }

//------------------------------------------------------------------------------

        real
        Wall::average_surface_temperature()
        {
            // grab the field from the mesh
            Vector< real > & tT = mMesh.field_data( "T" );

            Vector< real > tPhi( 3 );

            real tValue = 0.0 ;

            // element counter
            index_t e = 0;

            // loop over all elements
            for ( mesh::Element * tElement : mElements )
            {
                // collect field data
                for( uint i=0; i<3; ++i )
                {
                    tPhi( i ) = tT( tElement->node( i )->index() );
                }

                tValue += dot( mIntegrationWeights, tPhi ) * mElementLengths( e ) ;

                // integrate element counter
                e+=1 ;
            }

            return ( tValue /= mSegmentLength ) ;
        }

//------------------------------------------------------------------------------

        real
        Wall::average_heatload( const real & aAlpha, const real & aTinf )
        {
            // grab the field from the mesh
            Vector< real > & tTw  = mMesh.field_data( "T" );

            // heatload
            Vector< real > & tDotQ = mMesh.field_data( "dotQ" );

            // reference temperature for fluid
            Vector< real > & tTinf = mMesh.field_data( "Tinf" );

            // alpha value
            Vector< real > & tAlpha = mMesh.field_data( "alpha" );

            Vector< real > tPhi( 3 );


            // node / element counter
            index_t k ;

            // set nodal values
            for( mesh::Node * tNode : mNodes )
            {
                // get node index
                k = tNode->index() ;

                // set alpha condition
                tAlpha( k ) = aAlpha ;

                // set fluid condition
                tTinf( k ) = aTinf ;

                // compute heatload
                tDotQ( k ) = aAlpha * ( aTinf - tTw( k ) );

            }

            // reset counter
            k = 0 ;

            // integrated value
            real tValue = 0.0 ;

            // loop over all elements
            for ( mesh::Element * tElement : mElements )
            {
                // collect field data
                for( uint i=0; i<3; ++i )
                {
                    tPhi( i ) = tDotQ( tElement->node( i )->index() );
                }

                tValue += dot( mIntegrationWeights, tPhi ) * mElementLengths( k ) ;

                // integrate element counter
                k+=1 ;
            }

            // return averaged heatload
            return tValue / mSegmentLength ;
        }

//------------------------------------------------------------------------------

        void
        Wall::set_flowdata( const real & aT, const real & aP, const real & aMa  )
        {
            // Alpha-Parameter for heatload
            Vector< real > & tT = mMesh.field_data( "T_fluid" );
            Vector< real > & tP = mMesh.field_data( "p_fluid" );
            Vector< real > & tMa = mMesh.field_data( "Ma_fluid" );

            for( mesh::Node * tNode : mNodes )
            {
                tT( tNode->index() )  = aT ;
                tP( tNode->index() )  = aP ;
                tMa( tNode->index() ) = aMa ;
            }
        }

//------------------------------------------------------------------------------

        void
        Wall::collect_nodes_from_mesh( const Vector< id_t > & aNodeIDs )
        {
            // allocate memory
            mNodes.set_size( mNumNodes, nullptr ) ;

            // grab nodes from mesh
            for ( index_t k=0; k<mNumNodes; ++k )
            {
                mNodes( k ) = mMesh.node( aNodeIDs( k ) );
            }
        }

//------------------------------------------------------------------------------

        void
        Wall::create_integration_elements()
        {
            // the factory that builds these elements
            mesh::ElementFactory tFactory ;

            // allocate memory
            mElements.set_size( mNumElements, nullptr ) ;

            // offset for nodes
            index_t tOff = 0 ;

            // create and link elements
            for( index_t e=0; e<mNumElements; ++e )
            {
                // create a new element
                mesh::Element * tElement = tFactory.create_element( ElementType::LINE3, e+1 );

                // link element
                tElement->insert_node( mNodes( tOff ),  0 );
                tElement->insert_node( mNodes( tOff+2 ), 1 );
                tElement->insert_node( mNodes( tOff+1 ), 2 );

                // increment offset
                tOff += 2 ;

                // add element to container
                mElements( e ) = tElement ;
            }

            // compute the element lengths
            mesh::compute_edge_lengths( 3, mElements, mElementLengths );

            // compute the sum of the segment
            mSegmentLength = sum( mElementLengths );
        }

//------------------------------------------------------------------------------

        void
        Wall::initialize_integration_weights()
        {
            // since the shape funcition is always a line3, we can create
            // and integrate the shape analytically using the numbering scheme
            //
            // ( 0 ) ---- ( 2 ) ---- ( 1 )
            //
            // where 0 <= xi <= 1
            // so that
            // N0 = 2*xi^2 - 3*xi + 1
            // N1 = 2*xi^2 - xi
            // N2 = - 4*xi^2 + 4*xi

            //
            // Therefore, we end up with the Simpson rule
            mIntegrationWeights.set_size( 3 );
            mIntegrationWeights( 0 ) = 1.0 ;
            mIntegrationWeights( 1 ) = 1.0 ;
            mIntegrationWeights( 2 ) = 4.0 ;
            mIntegrationWeights /= 6.0 ;
        }

//------------------------------------------------------------------------------
    }
}