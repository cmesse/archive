//
// Pocket Detection for Mesh Cleanup
// Implements Tarjan's algorithm with extensions for unicyclic component detection
//

#ifndef BELFEM_FN_GRAPH_POCKET_DETECTION_HPP
#define BELFEM_FN_GRAPH_POCKET_DETECTION_HPP
#include <algorithm>
#include <utility>       // For std::pair
#include <limits>

#include "cl_Graph_Vertex.hpp"
#include "cl_Cell.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace graph
    {
        struct ArticulationPointWork
        {
            DynamicBitset * mContains = nullptr ;
            DynamicBitset * mVisited = nullptr ;
            DynamicBitset * mInComponent = nullptr ;
            Vector< index_t > mParent ;
            Vector< index_t > mDepth ;

            ArticulationPointWork( const index_t aSize ) :
                mContains( new DynamicBitset( aSize ) ),
                mVisited( new DynamicBitset( aSize ) ),
                mInComponent( new DynamicBitset( aSize ) ),
                mParent( aSize, gNoIndex ),
                mDepth( aSize, gNoIndex )
            {

            }

            ~ArticulationPointWork()
            {
                delete mContains;
                delete mVisited;
                delete mInComponent;
            }
        };

        /**
         * Structure to hold detailed information about articulation points
         * and the components they separate
         */
        class ArticulationPointInfo
        {
            Vertex * mVertex ;
            Cell< Vertex * > & mGraph ;
            DynamicBitset * mIsUnicyclic = nullptr ;

            Vector< index_t > mComponentSizes;
            Vector< index_t > mComponentEdgeCounts;
            Vector< index_t > mCycleLength;


            Cell< Cell< Vertex * > > mSeparatedComponents;

            // work vectors
            DynamicBitset * mVisited ;
            DynamicBitset * mInComponent ;
            Vector< index_t > & mParent ;
            Vector< index_t > & mDepth ;

        public:

            ArticulationPointInfo(
                Vertex * aVertex,
                Cell< Vertex * > & aGraph,
                Cell< Cell< Vertex * > > aSeparatedComponents,
                ArticulationPointWork & aWork );

            ~ArticulationPointInfo()
            {
                delete mIsUnicyclic;
            }

            Vertex * vertex() const { return mVertex; }

            Cell< Cell< Vertex * > > & separated_components()
            {
                return mSeparatedComponents;
            }

            Cell< Vertex * > & separated_components( const index_t aIndex )
            {
                return mSeparatedComponents( aIndex );
            }

            bool
            unicyclic_test( const index_t aIndex ) const
            {
                return mIsUnicyclic->test( aIndex );
            }

            void
            unicylcic_set( const index_t aIndex )
            {
                mIsUnicyclic->set( aIndex );
            }

            void
            unicyclic_reset( const index_t aIndex )
            {
                mIsUnicyclic->reset( aIndex );
            }

            index_t &
            component_size( const index_t aIndex )
            {
                return mComponentSizes( aIndex );
            }

            index_t &
            component_edge_count( const index_t aIndex )
            {
                return mComponentEdgeCounts( aIndex );
            }

            index_t &
            cycle_length( const index_t aIndex )
            {
                return mCycleLength( aIndex );
            }

            index_t
            size() const
            {
                return mSeparatedComponents.size();
            }

        private:

            std::pair< bool, index_t >
            detect_cycle_properties( const index_t aIndex );
        };

        /**
         * Find articulation points using Tarjan's algorithm
         * @param aGraph The graph represented as Cell of Vertex pointers
         * @return Cell containing all articulation points
         */
        //Cell< Vertex * >
        //find_articulation_points( Cell< Vertex * > & aGraph );

        /**
         * Detect cycle properties in a component
         * @param aComponent The vertices in the component
         * @return Tuple of (has_cycle, cycle_length)
         */
        std::pair< bool, index_t >
        detect_cycle_properties( const Cell< Vertex * > & aComponent );

        /**
         * Extended articulation point finder with component analysis
         * @param aGraph The graph to analyze
         * @return Detailed information about each articulation point
         */
        Cell< ArticulationPointInfo * >
        find_articulation_points_with_components( Cell< Vertex * > & aGraph ) ;

        /**
         * Structure representing a detected pocket
         */
        struct PocketInfo
        {
            Vertex * mNeckVertex = nullptr;     // The articulation point
            Cell< Vertex * > mPocketVertices;   // Vertices in the pocket
            index_t mSize = 0;                  // Number of vertices
            index_t mEdgeCount = 0;             // Number of edges
            index_t mCycleLength = 0;           // Length of the cycle (if unicyclic)
            bool mIsUnicyclic = false;          // True if exactly one cycle

            // Default constructor
            PocketInfo() = default;

            // Move constructor
            PocketInfo( PocketInfo&& other ) noexcept
                : mNeckVertex( other.mNeckVertex )
                , mPocketVertices( std::move( other.mPocketVertices ) )
                , mSize( other.mSize )
                , mEdgeCount( other.mEdgeCount )
                , mCycleLength( other.mCycleLength )
                , mIsUnicyclic( other.mIsUnicyclic )
            {
                other.mNeckVertex = nullptr;
                other.mSize = 0;
                other.mEdgeCount = 0;
                other.mCycleLength = 0;
                other.mIsUnicyclic = false;
            }

            // Move assignment operator
            PocketInfo& operator=( PocketInfo&& other ) noexcept
            {
                if ( this != &other )
                {
                    mNeckVertex = other.mNeckVertex;
                    mPocketVertices = std::move( other.mPocketVertices );
                    mSize = other.mSize;
                    mEdgeCount = other.mEdgeCount;
                    mCycleLength = other.mCycleLength;
                    mIsUnicyclic = other.mIsUnicyclic;

                    other.mNeckVertex = nullptr;
                    other.mSize = 0;
                    other.mEdgeCount = 0;
                    other.mCycleLength = 0;
                    other.mIsUnicyclic = false;
                }
                return *this;
            }

            // Copy constructor and assignment could be defaulted if needed
            PocketInfo( const PocketInfo& ) = default;
            PocketInfo& operator=( const PocketInfo& ) = default;

            // Compute a quality score for pocket removal priority
            real quality_score() const
            {
                if ( !mIsUnicyclic ) return 0.0;

                // Smaller pockets with shorter cycles are higher priority
                real tSizeScore = 1.0 / ( 1.0 + mSize );
                real tCycleScore = 1.0 / ( 1.0 + mCycleLength );
                real tCompactnessScore = static_cast< real >( mEdgeCount ) /
                                          ( mSize * ( mSize - 1 ) / 2.0 );

                return tSizeScore + tCycleScore + tCompactnessScore;
            }
        };

        /**
         * Detect all pockets in the graph
         * @param aGraph The graph to analyze
         * @param aMaxPocketSize Maximum size for a valid pocket
         * @param aRequireUnicyclic Whether to require exactly one cycle
         * @return Collection of detected pockets with metadata
         */
        Cell< PocketInfo >
        detect_pockets_with_info( Cell< Vertex * > & aGraph,
                                 const index_t aMaxPocketSize = 20,
                                 const bool aRequireUnicyclic = true );

        /**
         * Remove detected pockets from the graph
         * @param aGraph The graph to modify
         * @param aPockets The pockets to remove
         * @return Number of vertices removed
         */
        index_t
        remove_pockets( Cell< Vertex * > & aGraph,
                       const Cell< PocketInfo > & aPockets );

    } // namespace graph
} // namespace belfem

#endif // BELFEM_FN_GRAPH_POCKET_DETECTION_HPP