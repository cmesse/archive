//
// Created by Christian Messe on 07.12.20.
//

#include "cl_IsotropicChannel.hpp"
#include "cl_CH_Factory.hpp"
#include "fn_gesv.hpp"
#include "intpoints.hpp"
#include "fn_dot.hpp"
#include "fn_BL_ReferenceTempertaure.hpp"
#include "cl_Element_Factory.hpp"
#include "fn_Mesh_compute_edge_lengths.hpp"
#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_gesv.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_linspace.hpp"
#include "cl_OneDMapper.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    IsotropicChannel::IsotropicChannel( const IsotropicChannelType aType,
                                        const channel::BoundaryLayerMethod aMethod,
                                        HDF5 & aDatabase,
                                        Gas  & aGas,
                                        Mesh * aMesh1,
                                        Mesh * aMesh2 ) :
        mType( aType ),
        mGas( aGas ),
        mMesh1( *aMesh1 )
    {
        // remember the molar fractions
        mInitialMolarFractions = mGas.molar_fractions() ;

        // create a factory
        channel::Factory tFactory( aDatabase );

        switch( aType )
        {
            case ( IsotropicChannelType::CylindricChamber ) :
            {
                // create the geometry object
                mGeometry = tFactory.create_cylinder_geometry();

                // create the cylinder object
                tFactory.create_cylinder_segments( mGeometry, mMesh1, mSegments, false );

                mStatesFuncition = & IsotropicChannel::compute_states_chamber;

                break ;
            }
            case( IsotropicChannelType::Nozzle ) :
            {
                // create the geometry object
                mGeometry = tFactory.create_nozzle_geometry() ;

                // crate the nozzle segments
                tFactory.create_nozzle_segments( mGeometry, mMesh1, mSegments );

                mStatesFuncition = & IsotropicChannel::compute_states_nozzle ;

                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "Invalid Channel Mode" );
            }
        }

        mNumberOfChannels = 1 ;

        mMolarFractions.set_size( mSegments.size(),
                                  Vector< real >( mGas.number_of_components() ) );

        // create the boundary layer object
        mBoundaryLayer = new channel::Boundarylayer( aGas,
                                                     aMethod,
                                                     channel::SigmaRecoveryMode::Petrukov,
                                                     100 , 1.05 ); // 60 1.25

        uint tNumSegments = mSegments.size() ;

        // allocate spline data
        mHeatData.set_size( tNumSegments, mBoundaryLayer->heat_spline()->matrix_data() );
        mViscosityData.set_size( tNumSegments, mBoundaryLayer->viscosity_spline()->matrix_data() );
        mConductivityData.set_size( tNumSegments, mBoundaryLayer->conductivity_spline()->matrix_data() );

        this->set_reverse_order_flag( false );
    }


//------------------------------------------------------------------------------

    IsotropicChannel::~IsotropicChannel()
    {
        if( mBoundaryLayer != nullptr )
        {
            delete mBoundaryLayer ;
        }

        for( channel::Segment * tSegment : mSegments )
        {
            delete tSegment ;
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::print()
    {
        for( channel::Segment * tSegment : mSegments )
        {
            tSegment->print();
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::pull_temperatures()
    {
        for( channel::Segment * tSegment : mSegments )
        {
            tSegment->pull_surface_temperatures();
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::push_heatloads()
    {

        for( channel::Segment * tSegment : mSegments )
        {
            tSegment->push_heatloads();
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::push_flowdata()
    {
        for( channel::Segment * tSegment : mSegments )
        {
            tSegment->push_flowdata();
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::compute_massflow(
            const real & aTt,
            const real & aPt,
                  real & aT,
                  real & aP,
                  real & aDotm )
    {
        // remix the gas
        mGas.remix( mInitialMolarFractions, true, false );

        // relaxation factor
        real tOmega = 0.5 ;

        // energy state
        real tHt = mGas.h( aTt, aPt );

        // entropy state
        real tS  = mGas.s( aTt, aPt );

        // temperature
        aT = aTt ;

        // pressure
        aP = aPt ;

        // flow velocity
        real tU = 0  ;

        // cross section of last channel
        Vector< real > tX( 2 );

        // the Jacobian
        Matrix< real > tJ( 2, 2 );

        // the pivot
        Vector< int > tPivot( 2 );

        // cancel criterion
        real tError = BELFEM_REAL_MAX ;

        // infinite loop counter
        uint tCount = 0 ;

        while( tError > 1e-6 )
        {
            // compute the new flow velocity, assuming critical flow
            tU = mGas.c( aT, aP ) ;

            // enthalpy error
            tX( 0 ) = mGas.h( aT, aP ) + 0.5 * tU * tU - tHt ;
            tX( 1 ) = mGas.s( aT, aP ) - tS ;

            // compute the residual error
            tError = std::sqrt(   tX( 0 ) * tX( 0 ) / ( tHt * tHt )
                                  + tX( 1 ) * tX( 1 ) / ( tS * tS ) );

            // compute the Jacobian
            tJ( 0, 0 ) = mGas.cp( aT, aP ) ;
            tJ( 1, 0 ) = mGas.dsdT( aT, aP );
            tJ( 0, 1 ) = mGas.dhdp( aT, aP );
            tJ( 1, 1 ) = mGas.dsdp( aT, aP );

            // solve the system
            gesv( tJ, tX, tPivot );

            // perform correction step
            aT -= tOmega * tX( 0 );
            aP -= tOmega * tX( 1 );


            if ( mIsReacting )
            {
                mGas.remix( mInitialMolarFractions, false, false );
                mGas.remix_to_equilibrium( aT, aP, true, false );
            }

            BELFEM_ERROR( tCount++ < 100, "too many iterations" );
        }


        // get the cross section
        const real & tA = mType == IsotropicChannelType::Nozzle ?
                    mSegments( 0 )->cross_section() :
                    mSegments( mSegments.size() -1 )->cross_section() ;

        aDotm = tU * tA * mGas.rho( aT, aP );
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::compute_states(
            const real & aTt,
            const real & aPt )
    {
        ( this->*mStatesFuncition ) ( aTt, aPt );
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::compute_heatloads()
    {
        ( this->*mHeatloadsFunction )();
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::set_surface_roughness( const real & aRa )
    {
        mBoundaryLayer->set_surface_roughness( aRa );
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::set_bartz_geometry_params( const real & aHydraulicDiameter,
                               const real & aNozzleCurvature )
    {
        mBoundaryLayer->set_bartz_geometry_params( aHydraulicDiameter, aNozzleCurvature );
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::set_friction_method( const channel::BoundaryLayerMethod aMethod )
    {
        mBoundaryLayer->set_friction_method( aMethod );
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::compute_heatloads_forward()
    {
        uint tCount = 0 ;


        mGas.remix( mMolarFractions( 0 ), true, true );

        mBoundaryLayer->set_flow_conditions(
                mSegments( 0 )->value( BELFEM_CHANNEL_TM ),
                mSegments( 0 )->value( BELFEM_CHANNEL_PM ),
                mSegments( 0 )->value( BELFEM_CHANNEL_UM ) );

        mBoundaryLayer->set_center_conditions(
                mSegments( 0 )->value( BELFEM_CHANNEL_TM ),
                mSegments( 0 )->value( BELFEM_CHANNEL_UM ) );

        mBoundaryLayer->set_wall_temperature( mSegments( 0 )->value( BELFEM_CHANNEL_TW1 ) );

        mBoundaryLayer->set_hydraulic_diameter( mSegments( 0 )->value( BELFEM_CHANNEL_DH ) );

        mBoundaryLayer->compute_initial_guesses() ;

        mBoundaryLayer->use_input_from_parameters( true );

        for( channel::Segment * tSegment : mSegments )
        {

            // update data
            mGas.remix( mMolarFractions( tCount ), true, false );

            // load lookup tables from memory
            mBoundaryLayer->heat_spline()->matrix_data() = mHeatData( tCount );
            mBoundaryLayer->viscosity_spline()->matrix_data() = mViscosityData( tCount );
            mBoundaryLayer->conductivity_spline()->matrix_data() = mConductivityData( tCount );

            // passing the false flag tells the boundary layer module to not
            // recompute the lookup tables. We know that they are up to date now
            mBoundaryLayer->compute( tSegment->data() , false );

            ++tCount;
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::compute_heatloads_backward()
    {
        int tLast = mSegments.size() -1 ;

        mGas.remix( mMolarFractions( tLast ), true, true );

        mBoundaryLayer->set_flow_conditions(
                mSegments( tLast  )->value( BELFEM_CHANNEL_TM ),
                mSegments( tLast  )->value( BELFEM_CHANNEL_PM ),
                mSegments( tLast  )->value( BELFEM_CHANNEL_UM ) );

        mBoundaryLayer->set_center_conditions(
                mSegments( tLast )->value( BELFEM_CHANNEL_TM ),
                mSegments( tLast )->value( BELFEM_CHANNEL_UM ) );

        mBoundaryLayer->set_wall_temperature( mSegments( tLast )->value( BELFEM_CHANNEL_TW1 ) );

        mBoundaryLayer->set_hydraulic_diameter( mSegments( tLast )->value( BELFEM_CHANNEL_DH ) );

        mBoundaryLayer->compute_initial_guesses() ;

        mBoundaryLayer->use_input_from_parameters( true );

        for( int k=tLast; k>=0; --k )
        {
            // grab segment
            channel::Segment * tSegment = mSegments ( k );

            // update data
            mGas.remix( mMolarFractions( k ), true, false );

            mBoundaryLayer->heat_spline()->matrix_data() = mHeatData( k );
            mBoundaryLayer->viscosity_spline()->matrix_data() = mViscosityData( k );
            mBoundaryLayer->conductivity_spline()->matrix_data() = mConductivityData( k );

            mBoundaryLayer->compute( tSegment->data(), false );
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::compute_states_chamber(
            const real & aTt,
            const real & aPt )
    {
        // remix the gas
        mGas.remix( mInitialMolarFractions, true, false );

        // relaxation factor
        real tOmega = 0.3 ;

        // entropy
        real tS = mGas.s( aTt, aPt );

        // total enthalpy
        real tHt = mGas.h( aTt, aPt );

        // static Temperature
        real tT ;

        // static pressure
        real tP ;

        // massflow
        real tDotM ;

        // velocity
        real tU = 0 ;

        // density
        real tRho ;

        // enthalpy
        real tH ;

        // infinite loop protection counter
        uint tCount ;

        // Entropy Error
        real tError ;

        this->compute_massflow( aTt, aPt, tT, tP, tDotM );

        // number of segments
        uint tN = mSegments.size() ;

        for( int k=tN-1; k>=0; --k )
        {
            channel::Segment * tSegment = mSegments( k );

            tError = BELFEM_REAL_MAX;

            const real & tA = tSegment->cross_section();
            tCount = 0 ;

            while ( tError > 1e-6 )
            {
                // compute density
                tRho = mGas.rho( tT, tP );

                tU = tDotM / ( tRho * tA ) ;

                tH = tHt - 0.5 * tU * tU ;

                tT *= 1.0 - tOmega ;
                tT += tOmega * mGas.T_from_h( tH, tP );

                tError = std::abs( ( mGas.s( tT, tP ) - tS ) / tS  );

                tP -= tOmega * ( mGas.s( tT, tP ) - tS ) / mGas.dsdp( tT, tP ) ;

                // perform reaction
                if ( mIsReacting )
                {
                    mGas.remix( mInitialMolarFractions, false, false );
                    mGas.remix_to_equilibrium( tT, tP, true, false );
                }
                BELFEM_ERROR( tCount++ < 1000 , "too many iterations" );

            }

            // remix and compute all transport data
            mGas.remix_to_equilibrium( tT, tP, true, true );

            // store result
            tSegment->set_value( BELFEM_CHANNEL_TM, tT );
            tSegment->set_value( BELFEM_CHANNEL_PM, tP );
            tSegment->set_value( BELFEM_CHANNEL_UM, tU );
            tSegment->set_value( BELFEM_CHANNEL_MAM, tU / mGas.c( tT, tP ));
            tSegment->set_value( BELFEM_CHANNEL_RM, mGas.R( tT, tP ) );
            tSegment->set_value( BELFEM_CHANNEL_HM, mGas.h( tT, tP ) );
            tSegment->set_value( BELFEM_CHANNEL_SM, mGas.s( tT, tP ) );

            // store mixture
            mMolarFractions( k ) = mGas.molar_fractions() ;

            // compute the splines using this mixture
            mBoundaryLayer->update_lookup_tables() ;

            // store coeffcients
            mHeatData( k ) = mBoundaryLayer->heat_spline()->matrix_data() ;
            mViscosityData( k ) = mBoundaryLayer->viscosity_spline()->matrix_data() ;
            mConductivityData( k ) = mBoundaryLayer->conductivity_spline()->matrix_data() ;
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::compute_states_nozzle(
            const real & aTt,
            const real & aPt )
    {
        // remix the gas
        mGas.remix( mInitialMolarFractions, true, false );

        // relaxation factor
        real tOmega = 0.3 ;

        // entropy
        real tS = mGas.s( aTt, aPt );

        // total enthalpy
        real tHt = mGas.h( aTt, aPt );

        // static Temperature
        real tT ;

        // static pressure
        real tP ;


        // infinite loop protection counter
        uint tCount ;

        // Entropy Error
        real tError ;

        // massflow
        real tDotM ;

        this->compute_massflow( aTt, aPt, tT, tP, tDotM );

        // velocity
        real tU = mGas.c( tT, tP );

        // inverse density
        real tV = mGas.v( tT, tP );

        // enthalpy
        real tH ;

        // number of segments
        uint tN = mSegments.size() ;

        for( uint k=0; k<tN; ++k )
        {
            channel::Segment * tSegment = mSegments( k );

            tError = BELFEM_REAL_MAX;

            const real & tA = tSegment->cross_section();
            tCount = 0 ;

            while ( tError > 1e-6 )
            {
                tH = mGas.h( tT, tP );
                tU = std::sqrt( 2.0 * ( tHt - tH ) );
                tV =  tA * tU / tDotM ;
                tP = mGas.p( tT, tV );

                tT -= tOmega * ( mGas.s( tT, tP ) - tS ) / mGas.dsdT( tT, tP ) ;

                tError = std::abs( ( mGas.s( tT, tP ) - tS ) / tS  );

                // perform reaction
                if ( mIsReacting )
                {
                    mGas.remix( mInitialMolarFractions, false, false );
                    mGas.remix_to_equilibrium( tT, tP, true, false );
                }
                BELFEM_ERROR( tCount++ < 1000 , "too many iterations" );

            }

            // remix and compute all transport data
            mGas.remix_to_equilibrium( tT, tP, true, true );

            // store result
            tSegment->set_value( BELFEM_CHANNEL_TM, tT );
            tSegment->set_value( BELFEM_CHANNEL_PM, tP );
            tSegment->set_value( BELFEM_CHANNEL_UM, tU );
            tSegment->set_value( BELFEM_CHANNEL_MAM, tU / mGas.c( tT, tP ));
            tSegment->set_value( BELFEM_CHANNEL_RM, mGas.R( tT, tP ) );
            tSegment->set_value( BELFEM_CHANNEL_HM, mGas.h( tT, tP ) );
            tSegment->set_value( BELFEM_CHANNEL_SM, mGas.s( tT, tP ) );

            // store mixture
            mMolarFractions( k ) = mGas.molar_fractions() ;

            // compute the splines using this mixture
            mBoundaryLayer->update_lookup_tables() ;

            // store coeffcients
            mHeatData( k ) = mBoundaryLayer->heat_spline()->matrix_data() ;
            mViscosityData( k ) = mBoundaryLayer->viscosity_spline()->matrix_data() ;
            mConductivityData( k ) = mBoundaryLayer->conductivity_spline()->matrix_data() ;
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::set_reverse_order_flag( const bool aFlag )
    {
        mComputeInReverseOrder = aFlag ;

        // set the heatload function
        mHeatloadsFunction = aFlag ?
             & IsotropicChannel::compute_heatloads_backward :
             & IsotropicChannel::compute_heatloads_forward ;
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::compute_surface_coordinates()
    {
        // get the number of nodes
        index_t tNumNodes = mSegments.size() ;

        // number of elements along this mesh
        index_t tNumElems = ( tNumNodes - 1 ) / 2 ;

        // initialize a counter
        index_t tCount = 0 ;

        // coordinate for first node
        real tX0 = mSegments( 0 )->x() ;

        // create a temporary list of nodes from the segment coordinates
        Cell< mesh::Node * > tNodes( tNumNodes, nullptr );
        for( channel::Segment * tSegment : mSegments )
        {
            tNodes( tCount ) = new mesh::Node(
                    tCount+1,
                    tSegment->x()-tX0,
                    0.5 * tSegment->value( BELFEM_CHANNEL_DH ) );

            // increment the counter
            ++tCount;
        }

        // initialize an offset
        index_t tOff = 0 ;

        // create an element factory
        mesh::ElementFactory tFactory ;

        // create a temporary list of elements
        Cell< mesh::Element * > tElements( tNumElems, nullptr );
        for( index_t e=0; e<tNumElems; ++e )
        {
            // create a new element
            mesh::Element * tElement = tFactory.create_element(
                    ElementType::LINE3,
                    e + 1 );

            // link element with nodes
            tElement->insert_node( tNodes( tOff   ), 0 );
            tElement->insert_node( tNodes( tOff+2 ), 1 );
            tElement->insert_node( tNodes( tOff+1 ), 2 );

            // add element to container
            tElements( e ) = tElement ;

            // increment offset
            tOff += 2;
        }

        // container for edge length
        Vector< real > tEdgeLengths ;

        // now we can compute the lengths of each edge
        mesh::compute_edge_lengths( 2, tElements, tEdgeLengths );

        Vector< real > tSurfaceCoords( tNumNodes );

        // reset offset
        tOff = 0.0 ;

        tSurfaceCoords( 0 ) = 0.0 ;

        // reset counter
        tCount = 0 ;

        for( mesh::Element * tElement : tElements )
        {
            // compute weighted distance of center point
            real tLa = std::sqrt(
                      std::pow( tElement->node( 2 )->x() - tElement->node( 0 )->x(), 2 )
                    + std::pow( tElement->node( 2 )->y() - tElement->node( 0 )->y(), 2 ) );
            real tLb = std::sqrt(
                    std::pow( tElement->node( 1 )->x() - tElement->node( 2 )->x(), 2 )
                    + std::pow( tElement->node( 1 )->y() - tElement->node( 2 )->y(), 2 ) );

            tSurfaceCoords( tOff + 1 ) = tSurfaceCoords( tOff ) + tLa / ( tLa + tLb ) * tEdgeLengths( tCount );
            tSurfaceCoords( tOff + 2 ) = tSurfaceCoords( tOff ) + tEdgeLengths( tCount ) ;
            ++tCount ;
            tOff += 2 ;
        }

        // write coordinates into segments
        tCount = 0 ;
        for( channel::Segment * tSegment : mSegments )
        {
            tSegment->set_value( BELFEM_CHANNEL_S, tSurfaceCoords( tCount++ ) );
        }

        // delete elements
        for( mesh::Element * tElement : tElements )
        {
            delete tElement ;
        }

        // delete nodes
        for( mesh::Node * tNode : tNodes )
        {
            delete tNode ;
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::load_moc_data( const string & aFilePath )
    {

        // open an HDF5 file
        HDF5 tDatabase( aFilePath, FileMode::OPEN_RDONLY );

        // load the data for X and Ma
        Vector< real > tX0 ;
        Vector< real > tMa0 ;
        tDatabase.select_group( "Characteristics" );
        tDatabase.load_data( "x", tX0 );
        tDatabase.load_data( "Ma", tMa0 );

        // now we must map these mach numbers onto the segments
        // first, we reallocate the vectors that we have declared
        index_t tN1 = mSegments.size() ;
        Vector< real > tX1( tN1 );
        Vector< real > tMa1( tN1 );

        // now we grab the x-cordinates from the segments
        index_t tCount = 0 ;
        for( channel::Segment * tSegment : mSegments )
        {
            tX1( tCount++ ) = tSegment->x();
        }
        tX1 -= mSegments( 0 )->x();

        // create the mapper
        OneDMapper tMapper( tX1, 1 );

        // project data onto mapper
        tMapper.project( tX0, tMa0, tMa1 );

        // write result into segments
        tCount = 0 ;
        for( channel::Segment * tSegment : mSegments )
        {
            tSegment->set_value( BELFEM_CHANNEL_MAM, tMa1( tCount++ ) ) ;
        }
    }
//------------------------------------------------------------------------------

    void
    IsotropicChannel::compute_states_nozzle_from_characteristics(
            const real & aTt,
            const real & aPt )
    {
        // remix the gas
        mGas.remix( mInitialMolarFractions, true, false );

        // relaxation factor
        real tOmega = 0.3 ;

        // entropy
        real tSt = mGas.s( aTt, aPt );

        // total enthalpy
        real tHt = mGas.h( aTt, aPt );

        real tError ;

        // counter
        uint tCount ;

        // initial guesses for the first state
        real tGamma = mGas.gamma( aTt, aPt ) ;
        real tMa0 = mSegments( 0 )->value( BELFEM_CHANNEL_MAM );

        real tT0 = aTt / ( 1.0 + 0.5 * ( tGamma - 1.0 ) * tMa0 * tMa0 );
        real tP0 = aPt * std::pow( tT0 / aTt, tGamma / ( tGamma - 1.0 ) ) ;

        // index for segment
        index_t k = 0 ;

        for ( channel::Segment * tSegment : mSegments )
        {
            // get the data container
            Vector< real > & tData = tSegment->data();

            // static temperature
            real & tT = tData( BELFEM_CHANNEL_TM );

            // static pressure
            real & tP = tData( BELFEM_CHANNEL_PM );

            // velocity at wall
            real & tU = tData( BELFEM_CHANNEL_UM );

            // enthalpy
            real & tH = tData( BELFEM_CHANNEL_HM );

            // entropy
            real & tS = tData( BELFEM_CHANNEL_SM );

            // get the mach number
            const real & tMa = tSegment->value( BELFEM_CHANNEL_MAM );

            // guess state
            tT = tT0;
            tP = tP0;

            // reset error
            tError = BELFEM_REAL_MAX;

            // reset counter
            tCount = 0;

            // The Jacobian
            Matrix< real > tJ( 2, 2 ) ;

            // the right hand side
            Vector< real > tX( 2 ) ;

            // pivot vectir
            Vector< int > tPivot( 2 );

            // molar fractions from the last step
            Vector< real > tLastMolarFractions = mInitialMolarFractions;

            while ( tError > 1e-6 )
            {
                tU = tMa * mGas.c( tT, tP );

                tH = mGas.h( tT, tP ) ;
                tS = mGas.s( tT, tP );

                tJ( 0, 0 ) = mGas.cp( tT, tP ) ; // dhdt
                tJ( 0, 1 ) = 0.0 ; // dhdp, since this is an ideal gas
                tJ( 1, 0 ) = mGas.dsdT( tT, tP );
                tJ( 1, 1 ) = mGas.dsdp( tT, tP );

                tX( 0 ) = tH + 0.5 * tU * tU - tHt ;
                tX( 1 ) = tS - tSt ;

                tError = std::sqrt( std::pow( tX( 0 ) / tHt, 2 ) + std::pow( tX( 1 ) / tSt , 2 ) );

                // std::cout << "   " << tCount << " " << tT << " " << tP*1e-5 << " " << tU << " " << tError << " " << std::endl ;

                // solve system
                gesv( tJ, tX, tPivot );

                // adapt values
                tT -= tOmega * tX( 0 );
                tP -= tOmega * tX( 1 );

                // perform reaction
                if ( mIsReacting )
                {
                    mGas.remix( tLastMolarFractions, false, false );
                    mGas.remix_to_equilibrium( tT, tP, true, false );
                }
                BELFEM_ERROR( tCount++ < 1000, "too many iterations" );
            }

            // std::cout << tSegment->x() << " " << tT << " " << tP*1e-5 << " " << tMa << " " << tCount << std::endl ;
            // save gas constant
            real & tR = tData( BELFEM_CHANNEL_RM );
            tR = mGas.R( tT, tP );

            // shift state
            tT0 = tT;
            tP0 = tP;

            // save data into boundary layer, also updates the splines
            mBoundaryLayer->set_flow_conditions( tT, tP, tU );

            // store mixture
            mMolarFractions( k ) = mGas.molar_fractions() ;

            // store coeffcients
            mHeatData( k ) = mBoundaryLayer->heat_spline()->matrix_data() ;
            mViscosityData( k ) = mBoundaryLayer->viscosity_spline()->matrix_data() ;
            mConductivityData( k ) = mBoundaryLayer->conductivity_spline()->matrix_data() ;

            // shift moler fractions
            tLastMolarFractions = mGas.molar_fractions() ;

            // increment counter
            ++k ;
        }
    }

//------------------------------------------------------------------------------

    void
    IsotropicChannel::compute_pressure_derivatives()
    {
        // number of segments
        index_t tN = mSegments.size() ;

        // vector for surface coordinates
        Vector< real > tX( tN );

        // vector for pressures
        Vector< real > tP( tN );

        // vector for pressure derivatives
        Vector< real > tdPdX( tN );

        // collect data
        index_t k=0 ;
        for( channel::Segment * tSegment : mSegments )
        {
            tX( k ) = tSegment->value( BELFEM_CHANNEL_S );
            tP( k ) = tSegment->value( BELFEM_CHANNEL_PM );
            ++k ;
        }

        // create a mapper
        OneDMapper tMapper( tX, 2 );

        // compute derivatives
        tMapper.derive( tP, tdPdX );

        // write data back
        k=0 ;
        for( channel::Segment * tSegment : mSegments )
        {
            tSegment->set_value( BELFEM_CHANNEL_DPDS, tdPdX( k ) );
        }
    }

//------------------------------------------------------------------------------
/*
    real
    IsotropicChannel::compute_first_segment( const real & aXoff )
    {
        // mBoundaryLayer->set_x_off( aXoff ) ;

        // get first segment
        channel::Segment * tSegment = mSegments( 0 );

        mBoundaryLayer->heat_spline()->matrix_data() = mHeatData( 0 );
        mBoundaryLayer->viscosity_spline()->matrix_data() = mViscosityData( 0 );
        mBoundaryLayer->conductivity_spline()->matrix_data() = mConductivityData( 0 );

        mBoundaryLayer->compute( tSegment->data(), false );

        return tSegment->value( BELFEM_CHANNEL_DOTQ );
    }

//------------------------------------------------------------------------------

    real
    IsotropicChannel::compute_segment( const index_t & aIndex )
    {

        // get first segment
        channel::Segment * tSegment = mSegments( aIndex );

        mBoundaryLayer->heat_spline()->matrix_data() = mHeatData( aIndex );
        mBoundaryLayer->viscosity_spline()->matrix_data() = mViscosityData( aIndex );
        mBoundaryLayer->conductivity_spline()->matrix_data() = mConductivityData( aIndex );

        mBoundaryLayer->compute( tSegment->data(), false );

        return tSegment->value( BELFEM_CHANNEL_DOTQ );
    }*/

//------------------------------------------------------------------------------
}