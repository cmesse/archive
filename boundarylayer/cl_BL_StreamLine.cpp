//
// Created by Christian Messe on 13.06.20.
//
#include "assert.hpp"
#include "cl_BL_StreamLine.hpp"
#include "meshtools.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"

#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_dot.hpp"

#include "fn_Mesh_compute_edge_lengths.hpp"
#include "fn_linspace.hpp"

namespace belfem
{
    namespace boundarylayer
    {
//----------------------------------------------------------------------------

        StreamLine::StreamLine(
                StagnationPoint & aStagnationPoint,
                Mesh * aMesh, const id_t aID, Cell< mesh::Element * > & aElements ) :
            mStagnationPoint( aStagnationPoint ),
            mGas( aStagnationPoint.gas() ),
            mFreestream( aStagnationPoint.freestream() ),
            mStagnation( aStagnationPoint.stagnation() ),
            mMesh( aMesh ),
            mID( aID ),
            mElementType( aElements( 0 )->type() )
        {
            // make sure that the passed mesh is OK
            this->check_element_sanity( aElements );

            this->create_panels( aElements );
        }

//----------------------------------------------------------------------------

        StreamLine::~StreamLine()
        {
            for( Panel * tPanel : mPanels )
            {
                delete tPanel ;
            }
        }

//----------------------------------------------------------------------------
// FUNCTIONS NEEDED FOR CONSTRUCTION
//----------------------------------------------------------------------------

        void
        StreamLine::check_element_sanity( Cell< mesh::Element * > & aElements )
        {
            // make sure that the type is correct
            BELFEM_ERROR( mesh::geometry_type( mElementType ) == GeometryType::LINE,
                         "Errow while creating streamline %lu: passed elements must all be of type LINE",
                         ( long unsigned int ) mID );

            // check continuity of mesh
            BELFEM_ERROR( this->check_element_continuity( aElements ),
                         "Errow while creating streamline %lu: passed elements must continuously connect.",
                         ( long unsigned int ) mID );

            // check element type
            BELFEM_ERROR( this->check_element_types( aElements ),
                         "Errow while creating streamline %lu: passed elements must all be of the same type.",
                         ( long unsigned int ) mID );

        }

//----------------------------------------------------------------------------
        bool
        StreamLine::check_element_continuity( Cell< mesh::Element * > & aElements )
        {
            // Node ID of first node
            id_t tNodeID = aElements( 0 )->node( 0 )->id();

            uint tNumNodes = mesh::number_of_nodes( mElementType );

            // Container for Nodes
            Cell< mesh::Node * > tNodes( tNumNodes, nullptr );

            for ( mesh::Element * tElement : aElements )
            {
                if ( tElement->node( 0 )->id() != tNodeID )
                {
                    // check if element can be flipped
                    // collect nodes
                    for ( uint k = 0; k < tNumNodes; ++k )
                    {
                        tNodes( k ) = tElement->node( k );
                    }
                }
                tNodeID = tElement->node( 1 )->id();
            }

            return true;
        }

//----------------------------------------------------------------------------

        bool
        StreamLine::check_element_types( Cell< mesh::Element * > & aElements )
        {
            for ( mesh::Element * tElement : aElements )
            {
                if ( tElement->type() != mElementType )
                {
                    return false;
                }
            }

            return true;
        }

//----------------------------------------------------------------------------

        // create the panel objects for this streamline
        void
        StreamLine::create_panels( Cell< mesh::Element * > & aElements )
        {
            // cell with nodes
            Cell< mesh::Node * > tNodes;

            // direction of streamline
            Matrix< real > tDirections ;
            Matrix< real > tNormals ;
            Vector< real > tSurfaceCoordinates ;

            this->collect_nodes( aElements, tNodes );

            this->compute_direction_vectors( aElements, tDirections );
            this->collect_node_normals( tNodes, tNormals );
            this->compute_surface_coordinates( aElements, tSurfaceCoordinates );

            // number of panels
            index_t tNumPanels = tSurfaceCoordinates.length() ;

            // allocate memory
            mPanels.set_size( tNumPanels, nullptr );

            for( index_t k=0; k<tNumPanels; ++k )
            {
                mPanels( k ) = new Panel(
                        mGas,
                        mFreestream,
                        mStagnation,
                        tNodes( k ),
                        tSurfaceCoordinates( k ),
                        tDirections.col( k ),
                        tNormals.col( k ) );
            }
        }

//----------------------------------------------------------------------------

        // compute data needed for the panels
        void
        StreamLine::compute_direction_vectors(
                Cell< mesh::Element * > & aElements,
                Matrix< real > & aDirections )
        {

            // get the parameter coordinates
            Cell< Matrix< real > > tdNdXi;
            this->populate_shape_derivative( tdNdXi );

            // compute number of nodes
            uint tNumNodesPerElement = tdNdXi.size();

            uint tNumNodes = ( tNumNodesPerElement - 1 ) * aElements.size() + 1;

            // populate directions vector
            aDirections.set_size( 3, tNumNodes, 0.0 );

            // offset
            index_t tOff = 0;

            // node coords per element
            Matrix< real > tNodeCoords( tNumNodesPerElement, 3 );

            // matrix with direction ( needed to make blaze work )
            Matrix< real > tDirectionMat( 3, 1 );

            // vector with direction
            Vector< real > tDirection( 3 );

            for ( mesh::Element * tElement : aElements )
            {
                for ( uint k = 0; k < tNumNodesPerElement; ++k )
                {
                    tNodeCoords( k, 0 ) = tElement->node( k )->x();
                    tNodeCoords( k, 1 ) = tElement->node( k )->y();
                    tNodeCoords( k, 2 ) = tElement->node( k )->z();
                }

                for ( uint k = 0; k < tNumNodesPerElement; ++k )
                {
                    tDirectionMat = trans( tdNdXi( k ) * tNodeCoords ) ;
                    tDirection = aDirections.col( tOff + k ) + tDirectionMat.col( 0 );

                    aDirections.set_col( tOff + k, tDirection );
                }

                tOff += tNumNodesPerElement - 1;
            }

            for ( uint k = 0; k < tNumNodes; ++k )
            {
                tDirection = aDirections.col( k );
                tDirection /= norm( tDirection );
                aDirections.set_col( k, tDirection );
            }
        }

//----------------------------------------------------------------------------

        void
        StreamLine::populate_shape_derivative( Cell< Matrix< real > > & adNdXi )
        {
            // how many nodes exist in one element
            uint tNumNodesPerElement = mesh::number_of_nodes( mElementType );

            // container for parameter coordinates
            Matrix< real > tXi;

            // allocate memory for container
            tXi.set_size( 1, tNumNodesPerElement );

            // populate parameter coordinates
            switch ( mElementType )
            {
                case ( ElementType::LINE2 ) :
                {
                    tXi( 0, 0 ) = -1.0;
                    tXi( 0, 1 ) = 1.0;
                    break;
                }
                case ( ElementType::LINE3 ) :
                {
                    tXi( 0, 0 ) = -1.0;
                    tXi( 0, 1 ) = 1.0;
                    tXi( 0, 2 ) = 0.0;
                    break;
                }
                case ( ElementType::LINE4 ) :
                {
                    tXi( 0, 0 ) = -1.0;
                    tXi( 0, 1 ) = 1.0;
                    tXi( 0, 2 ) = -1.0 / 3.0;
                    tXi( 0, 3 ) = 1.0 / 3.0;
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unsupported Element Type" );
                }
            }

            // factory for shape function interpolation
            fem::InterpolationFunctionFactory tFactory;

            // shape function for interpolation
            fem::InterpolationFunction * tShape = tFactory.create_lagrange_function( mElementType );

            // allocate memory for output
            adNdXi.set_size( tNumNodesPerElement, Matrix< real >( 1, tNumNodesPerElement ));

            // populate interpolation vectors
            for ( uint k = 0; k < tNumNodesPerElement; ++k )
            {
                tShape->dNdXi( tXi.col( k ), adNdXi( k ));
            }

            // delete the interpolation function
            delete tShape;
        }

//----------------------------------------------------------------------------

        void
        StreamLine::collect_nodes(
                Cell< mesh::Element * > & aElements,
                Cell< mesh::Node * > & aNodes )
        {
            // number of nodes per element
            uint tNumNodesPerElement = mesh::number_of_nodes( mElementType );

            // number of nodes
            index_t tNumNodes = aElements.size() * ( tNumNodesPerElement - 1 ) + 1;

            aNodes.set_size( tNumNodes, nullptr );

            index_t tCount = 0;

            switch ( mElementType )
            {
                case ( ElementType::LINE2 ) :
                {
                    for ( mesh::Element * tElement : aElements )
                    {
                        aNodes( tCount++ ) = tElement->node( 0 );
                    }
                    break;
                }
                case ( ElementType::LINE3 ) :
                {
                    for ( mesh::Element * tElement : aElements )
                    {
                        aNodes( tCount++ ) = tElement->node( 0 );
                        aNodes( tCount++ ) = tElement->node( 2 );
                    }
                    break;
                }
                case ( ElementType::LINE4 ) :
                {
                    for ( mesh::Element * tElement : aElements )
                    {
                        aNodes( tCount++ ) = tElement->node( 0 );
                        aNodes( tCount++ ) = tElement->node( 2 );
                        aNodes( tCount++ ) = tElement->node( 3 );
                    }
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal Element type." );
                }
            }

            // last node
            aNodes( tCount++ ) = aElements( aElements.size() - 1 )->node( 1 );
        }

//----------------------------------------------------------------------------

        void
        StreamLine::collect_node_normals(
                Cell< mesh::Node * > & aNodes,
                Matrix< real > & aNormals )
        {
            switch ( mMesh->number_of_dimensions())
            {
                case ( 2 ) :
                {
                    // make sure that normals exist
                    BELFEM_ERROR(
                            mMesh->field_exists( "SurfaceNormalsx" ) &&
                            mMesh->field_exists( "SurfaceNormalsy" ),
                            "Could not find surface normals for this mesh. Were they computed?" );

                    // grab normals from mesh
                    Vector< real > & tNormalX = mMesh->field_data( "SurfaceNormalsx" );
                    Vector< real > & tNormalY = mMesh->field_data( "SurfaceNormalsy" );

                    // allocate memory
                    aNormals.set_size( 3, aNodes.size());

                    // init counter
                    index_t tCount = 0;

                    // loop over all nodes
                    for ( mesh::Node * tNode : aNodes )
                    {
                        aNormals( 0, tCount ) = tNormalX( tNode->index());
                        aNormals( 1, tCount ) = tNormalY( tNode->index());
                        aNormals( 2, tCount ) = 0.0;
                        ++tCount;
                    }

                    break;
                }
                case ( 3 ) :
                {
                    // make sure that normals exist
                    BELFEM_ERROR(
                            mMesh->field_exists( "SurfaceeNormalsx" ) &&
                            mMesh->field_exists( "SurfaceeNormalsy" ) &&
                            mMesh->field_exists( "SurfaceeNormalsz" ),
                            "Could not find surface normals for this mesh. Were they computed?" );

                    Vector< real > & tNormalX = mMesh->field_data( "SurfaceeNormalsx" );
                    Vector< real > & tNormalY = mMesh->field_data( "SurfaceeNormalsy" );
                    Vector< real > & tNormalZ = mMesh->field_data( "SurfaceeNormalsz" );

                    // allocate memory
                    aNormals.set_size( 3, aNodes.size());

                    // init counter
                    index_t tCount = 0;

                    // loop over all nodes
                    for ( mesh::Node * tNode : aNodes )
                    {
                        aNormals( 0, tCount ) = tNormalX( tNode->index());
                        aNormals( 1, tCount ) = tNormalY( tNode->index());
                        aNormals( 2, tCount ) = tNormalZ( tNode->index());
                        ++tCount;
                    }

                    break;
                }
                default :
                {
                    BELFEM_ERROR( false, "Illegal mesh dimension %u.", ( uint ) mMesh->number_of_dimensions());
                }
            }

            // ceck sanity
            index_t tNumNodes = aNodes.size();

            // fix for first node
            Vector< real > tA( aNormals.col( 0 ));
            Vector< real > tB( aNormals.col( 1 ));

            // compute angle between both vectors
            if ( std::abs( dot( tA, tB )) < 0.99 )
            {
                // copy value of second node into first node
                aNormals.set_col( 0, tB );
            }

            // fix for last node
            tA = aNormals.col( tNumNodes - 1 );
            tB = aNormals.col( tNumNodes - 2 );

            // compute angle between both vectors
            if ( std::abs( dot( tA, tB )) < 0.99 )
            {
                // copy value of butlast node into last node
                aNormals.set_col( tNumNodes - 1, tB );
            }

            // check sanity of whole dataset
            for ( index_t k = 0; k < tNumNodes; ++k )
            {
                if ( std::abs( norm( aNormals.col( k )) - 1.0 ) > 1e-6 )
                {
                    Vector< real > tNorm = aNormals.col( k );
                    tNorm.print( "N" );
                    std::cout << norm( aNormals.col( k )) << std::endl;

                }
                BELFEM_ERROR( std::abs( norm( aNormals.col( k )) - 1.0 ) < 1e-6,
                             "Faulty node normal detected for node %lu at streamline %lu.",
                             ( long unsigned int ) aNodes( k )->id(),
                             ( long unsigned int ) mID );
            }
        }

//----------------------------------------------------------------------------

        void
        StreamLine::compute_surface_coordinates(
                Cell< mesh::Element * > & aElements,
                Vector< real > & aSurfaceCoordinates )
        {
            // compute the edge lengths of each element
            Vector< real > tEdgeLengths;
            mesh::compute_edge_lengths(
                    mMesh->number_of_dimensions(),
                    aElements,
                    tEdgeLengths );

           // compute number of nodes
           uint tNumNodesPerElement = mesh::number_of_nodes( mElementType );

           index_t tNumElems = aElements.size() ;
           index_t tNumNodes = tNumElems * ( tNumNodesPerElement - 1 ) + 1 ;

           // allocate memory for output
           aSurfaceCoordinates.set_size( tNumNodes );

           // local coordinate for each element
           Vector< real > tXi = linspace( 0.0, 1.0, tNumNodesPerElement );

           // initialize counter
           index_t tCount = 0 ;

           // initialize first node
           aSurfaceCoordinates( tCount++ ) = 0.0 ;

           // element offset
           real tOff = 0 ;

           // loop over all elements
           for( index_t e=0; e<tNumElems; ++e )
           {
                // get length of element
                real tLength = tEdgeLengths( e );

                // populate node coordiantes, starting with second node
                for( uint k=1; k<tNumNodesPerElement; ++k )
                {
                    aSurfaceCoordinates( tCount++ ) = tOff + tXi( k ) * tLength ;
                }

                // add length to offset
                tOff += tLength ;
           }
        }
//----------------------------------------------------------------------------
// FUNCTIONS NEEDED FOR SURFACE INCLINATION METHOD
//----------------------------------------------------------------------------

        void
        StreamLine::compute( const real & aAoA )
        {
            // compute the direction of the freestream flow
            Vector< real > tFreestreamDirection ;
            this->compute_freestream_direction( aAoA, tFreestreamDirection );

            // perform a modified newton
            this->compute_modified_newton( tFreestreamDirection );

            // split the streamline into lower and upper panels
            index_t tStagIndex = this->split_streamline( tFreestreamDirection );

            // perform the prandtl meyer expansion for both sides
            this->compute_prandtl_meyer( mLowerPanels ) ;
            this->compute_prandtl_meyer( mUpperPanels ) ;

            // compute the heatload at the stagnation point
            // compute the stagnation state
            // mStagnationPoint.compute_stagnation_heatload( mPanels( tStagIndex )->state()->Tw() );
            mStagnationPoint.compute_stagnation_heatload( mPanels( tStagIndex )->state()->Tw(),
                    mPanels( tStagIndex ) ->x() );

            this->print( mLowerPanels ) ;

        }

//----------------------------------------------------------------------------

        void
        StreamLine::compute_freestream_direction(
                const real     & aAoA,
                Vector< real > & aFreestreamDirection )
        {
            aFreestreamDirection.set_size( 3 );

            if( mMesh->number_of_dimensions() == 3 )
            {
                aFreestreamDirection( 0 ) = std::cos( aAoA * constant::deg );
                aFreestreamDirection( 1 ) = 0.0;
                aFreestreamDirection( 2 ) = std::sin( aAoA * constant::deg );
            }
            else
            {
                aFreestreamDirection( 0 ) = std::cos( aAoA * constant::deg );
                aFreestreamDirection( 1 ) = std::sin( aAoA * constant::deg );
                aFreestreamDirection( 2 ) = 0.0 ;
            }
        }

//----------------------------------------------------------------------------

        void
        StreamLine::compute_modified_newton(
                const Vector< real > & aFreestreamDirection )
        {

            // perform a modified newton method on all panels
            for ( Panel * tPanel : mPanels )
            {
                tPanel->unflag();
                tPanel->compute_aoa( aFreestreamDirection );
                tPanel->compute_newton();
            }
        }

//----------------------------------------------------------------------------

        // split the streamline into a lower and an upper part
        index_t
        StreamLine::split_streamline( Vector< real > & aFreestreamDirection )
        {
            // find the index of the stagnation point
            sint aStagIndex = this->find_stagnation_point() ;

            // surface coordinates for stagnation point
            real tX0 = this->compute_surface_coordinate( aFreestreamDirection, aStagIndex );
            real tS0 = mPanels( aStagIndex )->s() ;

            // with the surface coordinate at the reference point computed
            // we can compute all x-coordiates which are needed for the reference
            // reynolds number
            for ( Panel * tPanel : mPanels )
            {
                tPanel->set_x( std::abs( tPanel->s() - tS0 ) + tX0 );
            }


            // populate lower panels
            sint tNumPanels = mPanels.size() ;
            sint tNumLowerPanels = mPanels.size() - aStagIndex;
            mLowerPanels.set_size( tNumLowerPanels, nullptr );
            uint tCount = 0;
            for ( sint k = aStagIndex; k < tNumPanels; ++k )
            {
                mLowerPanels( tCount++ ) = mPanels( k );
            }

            // populate upper panels
            tCount = 0;
            sint tNumUpperPanels = aStagIndex + 1;
            mUpperPanels.set_size( tNumUpperPanels, nullptr );
            for ( int k = aStagIndex; k >= 0; --k )
            {
                mUpperPanels( tCount++ ) = mPanels( k );
            }

            return aStagIndex ;
        }

//----------------------------------------------------------------------------

        index_t
        StreamLine::find_stagnation_point()
        {
            // index of point
            index_t aIndex = gNoIndex ;

            // maximum pressure
            real    tPmax = 0.0 ;
            index_t tCount = 0 ;

            for ( Panel * tPanel : mPanels )
            {
                if ( tPanel->state()->p() > tPmax )
                {
                    tPmax = tPanel->state()->p();
                    aIndex = tCount;
                }
                tCount++;
            }

            return aIndex ;
        }

//----------------------------------------------------------------------------

        real
        StreamLine::compute_surface_coordinate(
                Vector< real > & aFreestreamDirection,
                const index_t    aStagnationIndex )
        {
            //  we take the stagnation point ant project it
            // on a unit sphere
            Vector< real > tX( 3 );
            tX( 0 ) = mPanels( aStagnationIndex )->node()->x() - mStagnationPoint.nose_radius() ;
            tX( 1 ) = mPanels( aStagnationIndex )->node()->y();
            tX( 2 ) = mPanels( aStagnationIndex )->node()->z();
            tX /= norm( tX );

            // this is the angle from this point to the actual stagnation point
            real tRefAngle = std::acos( dot( tX.vector_data(), -1 * aFreestreamDirection.vector_data() ));

            // this is the reference length at the stagnation point
            return mStagnationPoint.nose_radius() * tRefAngle;
        }

//----------------------------------------------------------------------------

        void
        StreamLine::compute_prandtl_meyer( Cell< Panel * > & aPanels )
        {
            index_t tNumPanels = aPanels.size();

            Panel * tPanel = aPanels( 0 );
            Panel * tPanel0;

            for ( index_t k = 1; k < tNumPanels; ++k )
            {
                // shift panels
                tPanel0 = tPanel;
                tPanel = aPanels( k );
                if ( tPanel0->state()->Ma() > 1.01 )
                {

                    // delta angle for prandtl meyer
                    real tNu = tPanel0->aoa() - tPanel->aoa();

                    if ( std::abs( tNu ) > 1e-7 )
                    {
                        tPanel->compute_prandtl_meyer(
                                tPanel0->state()->T(),
                                tPanel0->state()->p(),
                                tPanel0->state()->u(),
                                tNu );
                    }
                    else
                    {
                        // copy state
                        tPanel->state()->compute(
                                tPanel0->state()->T(),
                                tPanel0->state()->p(),
                                tPanel0->state()->u());
                    }

                    tPanel->state()->compute_Cp( mFreestream );
                }
            }
        }

//----------------------------------------------------------------------------

        void
        StreamLine::print()
        {
            std::cout << "numLowerPanels " <<mLowerPanels.size() << std::endl ;

            this->print( mLowerPanels );
        }

//----------------------------------------------------------------------------
        void
        StreamLine::print( Cell< Panel * > & aPanels )
        {
            Matrix< real > tData( aPanels.size(), 5 );

            // surface coordinate
            index_t tCount = 0 ;

            Vector< real > & tDotQ = mMesh->field_data("dotQ");

            for( Panel * tPanel : aPanels )
             {
                tData( tCount, 0 ) = tPanel->x() ;
                tData( tCount, 1 ) = tPanel->state()->T() ;
                tData( tCount, 2 ) = tPanel->state()->p() ;
                tData( tCount, 3 ) = tPanel->state()->u() ;
                tData( tCount, 4 ) = tDotQ( tPanel->node()->index() );

                ++tCount ;
            }

            tData.print("Data");

        }

//----------------------------------------------------------------------------
    }
}