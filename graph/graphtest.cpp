/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <iostream>
#include <cassert>

#include "banner.hpp"
#include "cl_Communicator.hpp"
#include "cl_HDF5.hpp"
#include "cl_Logger.hpp"
#include "cl_Matrix.hpp"
#include "cl_Timer.hpp"
#include "cl_Vector.hpp"


#include "fn_Graph_max_cardinality_matching.hpp"

using namespace belfem;
using namespace belfem::graph;

/**
 * Helper function to create an undirected edge between two vertices.
 * Sets up bidirectional neighbor relationships.
 */
void
create_edge( Vertex * aA, Vertex * aB )
{
    aA->increment_vertex_counter();
    aB->increment_vertex_counter();
}

void
finalize_edges( Graph & aGraph )
{
    // Initialize vertex containers
    for ( Vertex * tVertex : aGraph )
    {
        tVertex->init_vertex_container();
    }
}

void
add_edge( Vertex * aA, Vertex * aB )
{
    aA->insert_vertex( aB );
    aB->insert_vertex( aA );
}

/**
 * Test 1: Simple path graph with 4 vertices
 *         0 -- 1 -- 2 -- 3
 * Expected: Maximum matching of size 2
 */
void
test_path_graph()
{
    std::cout << "Test 1: Path graph (4 vertices)..." << std::endl;

    Graph tGraph( 4, nullptr );

    // Create vertices
    for ( index_t k = 0; k < 4; ++k )
    {
        tGraph( k ) = new Vertex();
        tGraph( k )->set_index( k );
    }

    // Count edges first
    create_edge( tGraph( 0 ), tGraph( 1 ) );
    create_edge( tGraph( 1 ), tGraph( 2 ) );
    create_edge( tGraph( 2 ), tGraph( 3 ) );

    finalize_edges( tGraph );

    // Add edges
    add_edge( tGraph( 0 ), tGraph( 1 ) );
    add_edge( tGraph( 1 ), tGraph( 2 ) );
    add_edge( tGraph( 2 ), tGraph( 3 ) );

    // Run algorithm
    Cell< index_t > tMatch;
    index_t tCardinality = max_cardinality_matching( tGraph, tMatch );

    std::cout << "  Cardinality: " << tCardinality << " (expected: 2)" << std::endl;
    std::cout << "  Matching: ";
    for ( index_t k = 0; k < tMatch.size(); ++k )
    {
        if ( tMatch( k ) != gNoIndex )
        {
            std::cout << k << "->" << tMatch( k ) << " ";
        }
    }
    std::cout << std::endl;

    assert( tCardinality == 2 );

    // Cleanup
    for ( Vertex * tV : tGraph )
    {
        delete tV;
    }

    std::cout << "  PASSED" << std::endl;
}

/**
 * Test 2: Complete graph K4
 * Expected: Maximum matching of size 2 (perfect matching)
 */
void
test_complete_graph_k4()
{
    std::cout << "Test 2: Complete graph K4..." << std::endl;

    Graph tGraph( 4, nullptr );

    for ( index_t k = 0; k < 4; ++k )
    {
        tGraph( k ) = new Vertex();
        tGraph( k )->set_index( k );
    }

    // Count all edges
    for ( index_t i = 0; i < 4; ++i )
    {
        for ( index_t j = i + 1; j < 4; ++j )
        {
            create_edge( tGraph( i ), tGraph( j ) );
        }
    }

    finalize_edges( tGraph );

    // Add all edges
    for ( index_t i = 0; i < 4; ++i )
    {
        for ( index_t j = i + 1; j < 4; ++j )
        {
            add_edge( tGraph( i ), tGraph( j ) );
        }
    }

    Cell< index_t > tMatch;
    index_t tCardinality = max_cardinality_matching( tGraph, tMatch );

    std::cout << "  Cardinality: " << tCardinality << " (expected: 2)" << std::endl;
    assert( tCardinality == 2 );

    for ( Vertex * tV : tGraph )
    {
        delete tV;
    }

    std::cout << "  PASSED" << std::endl;
}

/**
 * Test 3: Triangle (3 vertices)
 * Expected: Maximum matching of size 1
 */
void
test_triangle()
{
    std::cout << "Test 3: Triangle (3 vertices)..." << std::endl;

    Graph tGraph( 3, nullptr );

    for ( index_t k = 0; k < 3; ++k )
    {
        tGraph( k ) = new Vertex();
        tGraph( k )->set_index( k );
    }

    create_edge( tGraph( 0 ), tGraph( 1 ) );
    create_edge( tGraph( 1 ), tGraph( 2 ) );
    create_edge( tGraph( 2 ), tGraph( 0 ) );

    finalize_edges( tGraph );

    add_edge( tGraph( 0 ), tGraph( 1 ) );
    add_edge( tGraph( 1 ), tGraph( 2 ) );
    add_edge( tGraph( 2 ), tGraph( 0 ) );

    Cell< index_t > tMatch;
    index_t tCardinality = max_cardinality_matching( tGraph, tMatch );

    std::cout << "  Cardinality: " << tCardinality << " (expected: 1)" << std::endl;
    assert( tCardinality == 1 );

    for ( Vertex * tV : tGraph )
    {
        delete tV;
    }

    std::cout << "  PASSED" << std::endl;
}

/**
 * Test 4: Bipartite graph (simple)
 *         0 -- 2
 *         |    |
 *         1 -- 3
 * Expected: Maximum matching of size 2
 */
void
test_bipartite()
{
    std::cout << "Test 4: Bipartite graph (4 vertices)..." << std::endl;

    Graph tGraph( 4, nullptr );

    for ( index_t k = 0; k < 4; ++k )
    {
        tGraph( k ) = new Vertex();
        tGraph( k )->set_index( k );
    }

    create_edge( tGraph( 0 ), tGraph( 2 ) );
    create_edge( tGraph( 0 ), tGraph( 3 ) );
    create_edge( tGraph( 1 ), tGraph( 2 ) );
    create_edge( tGraph( 1 ), tGraph( 3 ) );

    finalize_edges( tGraph );

    add_edge( tGraph( 0 ), tGraph( 2 ) );
    add_edge( tGraph( 0 ), tGraph( 3 ) );
    add_edge( tGraph( 1 ), tGraph( 2 ) );
    add_edge( tGraph( 1 ), tGraph( 3 ) );

    Cell< index_t > tMatch;
    index_t tCardinality = max_cardinality_matching( tGraph, tMatch );

    std::cout << "  Cardinality: " << tCardinality << " (expected: 2)" << std::endl;
    assert( tCardinality == 2 );

    for ( Vertex * tV : tGraph )
    {
        delete tV;
    }

    std::cout << "  PASSED" << std::endl;
}

/**
 * Test 5: Empty graph
 */
void
test_empty_graph()
{
    std::cout << "Test 5: Empty graph..." << std::endl;

    Graph tGraph;

    Cell< index_t > tMatch;
    index_t tCardinality = max_cardinality_matching( tGraph, tMatch );

    std::cout << "  Cardinality: " << tCardinality << " (expected: 0)" << std::endl;
    assert( tCardinality == 0 );
    assert( tMatch.size() == 0 );

    std::cout << "  PASSED" << std::endl;
}

/**
 * Test 6: Single vertex (no edges)
 */
void
test_single_vertex()
{
    std::cout << "Test 6: Single vertex..." << std::endl;

    Graph tGraph( 1, nullptr );
    tGraph( 0 ) = new Vertex();
    tGraph( 0 )->set_index( 0 );
    tGraph( 0 )->init_vertex_container( 0 );

    Cell< index_t > tMatch;
    index_t tCardinality = max_cardinality_matching( tGraph, tMatch );

    std::cout << "  Cardinality: " << tCardinality << " (expected: 0)" << std::endl;
    assert( tCardinality == 0 );

    delete tGraph( 0 );

    std::cout << "  PASSED" << std::endl;
}

/**
 * Test 7: Pentagon (odd cycle - tests blossom handling)
 *         0 -- 1
 *        /      \
 *       4        2
 *        \      /
 *         3 --
 * Expected: Maximum matching of size 2
 */
void
test_pentagon()
{
    std::cout << "Test 7: Pentagon (5 vertices)..." << std::endl;

    Graph tGraph( 5, nullptr );

    for ( index_t k = 0; k < 5; ++k )
    {
        tGraph( k ) = new Vertex();
        tGraph( k )->set_index( k );
    }

    create_edge( tGraph( 0 ), tGraph( 1 ) );
    create_edge( tGraph( 1 ), tGraph( 2 ) );
    create_edge( tGraph( 2 ), tGraph( 3 ) );
    create_edge( tGraph( 3 ), tGraph( 4 ) );
    create_edge( tGraph( 4 ), tGraph( 0 ) );

    finalize_edges( tGraph );

    add_edge( tGraph( 0 ), tGraph( 1 ) );
    add_edge( tGraph( 1 ), tGraph( 2 ) );
    add_edge( tGraph( 2 ), tGraph( 3 ) );
    add_edge( tGraph( 3 ), tGraph( 4 ) );
    add_edge( tGraph( 4 ), tGraph( 0 ) );

    Cell< index_t > tMatch;
    index_t tCardinality = max_cardinality_matching( tGraph, tMatch );

    std::cout << "  Cardinality: " << tCardinality << " (expected: 2)" << std::endl;
    assert( tCardinality == 2 );

    for ( Vertex * tV : tGraph )
    {
        delete tV;
    }

    std::cout << "  PASSED" << std::endl;
}

Communicator gComm;
Logger       gLog( 5 );


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    std::cout << "========================================" << std::endl;
    std::cout << "Maximum Cardinality Matching Tests" << std::endl;
    std::cout << "(Micali-Vazirani Algorithm)" << std::endl;
    std::cout << "========================================" << std::endl;

    test_empty_graph();
    test_single_vertex();
    test_triangle();
    test_path_graph();
    test_bipartite();
    test_complete_graph_k4();
    test_pentagon();

    std::cout << "========================================" << std::endl;
    std::cout << "All tests PASSED!" << std::endl;
    std::cout << "========================================" << std::endl;

    return gComm.finalize();
}
