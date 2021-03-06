#ifndef ADJACENCY_DI_GRAPH_TEST
#define ADJACENCY_DI_GRAPH_TEST

#include "testGroupInterface.hpp"
#include "adjacencyDiGraph.hpp"

#include <random>
#include <set>
#include <vector>

using namespace graph;

struct AdjacencyDiGraphTest : public ITestGroup
{
    // consts
    static const uint ITERATIONS = 100;
    static const uint TEST_SIZE = 100;

    // globals

    static std::default_random_engine mRandom;

    // helpers
    template< typename T
	      , template< typename K, typename V, typename CompT > typename VCT
	      , template< typename S > typename IECT
	      , template< typename S > typename OECT
	      , template< typename S > typename ValDistT
	      , typename ComparatorT = std::equal_to< T > >
    static AdjacencyDiGraph< T, VCT, IECT, OECT, ComparatorT > randomVertices( size_t number, ValDistT< T >& vdist )
    {
	AdjacencyDiGraph< T, VCT, IECT, OECT, ComparatorT > ag;
	for( uint i = 0; i < number; ++i )
	    addVertex( ag, vdist( mRandom ) );
	return ag;
    }

    template< typename T
	      , template< typename K, typename V, typename CompT > typename VCT
	      , template< typename S > typename IECT
	      , template< typename S > typename OECT
	      , typename ComparatorT = std::equal_to< T > >
    static AdjacencyDiGraph< T, VCT, IECT, OECT, ComparatorT > randomEdges( AdjacencyDiGraph< T, VCT, IECT, OECT, ComparatorT >& ag, size_t number )
    {
	typedef AdjacencyDiGraph< T, VCT, IECT, OECT, ComparatorT > G;
	typename DiGraphTraits< G >::VRangeT vs = vertices( ag );
	std::uniform_int_distribution<> jumpDist( 0, std::distance( vs.first, vs.second ) - 1 );
	for( uint i = 0; i < number; ++i )
	{
	    typename G::VIterT isrc = vs.first + jumpDist( mRandom );
	    typename G::VIterT itrg = vs.first + jumpDist( mRandom );
	    addEdge( ag, *isrc, *itrg );
	}
	return ag;
    }

    // test classes
    // vector acting as set
    template< typename T >
    using Gt = AdjacencyDiGraph< T, VecMap, InVec, InVec >;
    
    typedef AdjacencyDiGraph< int, VecMap, InVec, InVec > G;
    
    // test whether a value can be found in the graph after adding it
    class AddValueTest : public ITest
    {
	STATEFUL_TEST( AddValueTest );
      private:
	G mGraph;
	std::uniform_int_distribution<> mVdist;
	int mVal;
    };

    // test whether an edge with matching source and target can be found after adding it
    class AddEdgeTest : public ITest
    {
    	STATEFUL_TEST( AddEdgeTest );
      private:
    	G mGraph;
    	std::uniform_int_distribution<> mVdist;
    	int mSrcVal, mTrgVal;
    };

    // test whether removing vertex renders it absent from the graph
    class RemoveVertexTest : public ITest
    {
    	STATEFUL_TEST( RemoveVertexTest );
      private:
    	G mGraph;
    	std::uniform_int_distribution<> mVdist;
    	int mRemVal;
    };
    
    // test whether removing vertex does not leave invalid vertices
    class RemoveVertexNonCorruptionTest : public ITest
    {
    	STATEFUL_TEST( RemoveVertexNonCorruptionTest );
      private:
    	G mGraph;
    	std::vector< int > mRemain;
    	std::uniform_int_distribution<> mVdist;
    };

    // test whether removing vertex removes incoming connections of formerly adjacent vertices
    class RemoveVertexIncomingTest : public ITest
    {
    	STATEFUL_TEST( RemoveVertexIncomingTest );
      private:
    	G mGraph;
    	std::uniform_int_distribution<> mVdist;
    	int mRemVal;
    	std::vector< typename DiGraphTraits< G >::VertexT > mIns;
    };

    // test whether removing vertex removes incoming connections of formerly adjacent vertices
    class RemoveVertexOutgoingTest : public ITest
    {
    	STATEFUL_TEST( RemoveVertexOutgoingTest );
      private:
    	G mGraph;
    	std::uniform_int_distribution<> mVdist;
    	int mRemVal;
    	std::vector< typename DiGraphTraits< G >::VertexT > mIns;
    };

    // test whether removing vertex removes outgoing connectinos of formerly adjacent vertices
    class RemoveEdgeTest : public ITest
    {
    	STATEFUL_TEST( RemoveEdgeTest );
      private:
    	G mGraph;
    	std::uniform_int_distribution<> mJumpDist;
    	int mSrcVal, mTrgVal;
    };

    // test whether objects stored in graph are deleted as often as they are created
    class MemoryFreed : public ITest
    {
      public:
	class DummyValue
	{
	    uint& mCounterRef;
	    uint mId;
	  public:
	    DummyValue( uint& counterRef, const uint& id );
	    DummyValue( const DummyValue& orig );
	    DummyValue& operator =( const DummyValue& orig );
	    ~DummyValue();
	    bool operator ==( const DummyValue& ) const;
	};

	MemoryFreed( const uint&, const uint& );

	void iterate(); bool check() const;
      private:
	std::unique_ptr< Gt< DummyValue > > mpGraph;
	uint mObjectCounter = 0, mProduceCounter = 0;
	std::uniform_int_distribution<> mValDist, mSizeDist;
    };
    
    GROUP_CTOR_DECL( AdjacencyDiGraphTest );

    void init();
};

#endif
