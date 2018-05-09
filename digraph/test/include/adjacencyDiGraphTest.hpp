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
	      , template< typename K, typename V > typename VCT
	      , template< typename S > typename IECT
	      , template< typename S > typename OECT
	      , template< typename S > typename ValDistT >
    static AdjacencyDiGraph< T, VCT, IECT, OECT > randomVertices( size_t number, ValDistT< T >& vdist )
    {
	AdjacencyDiGraph< T, VCT, IECT, OECT > ag;
	for( uint i = 0; i < number; ++i )
	    addVertex( ag, vdist( mRandom ) );
	return ag;
    }

    template< typename T
	      , template< typename K, typename V > typename VCT
	      , template< typename S > typename IECT
	      , template< typename S > typename OECT >
    static AdjacencyDiGraph< T, VCT, IECT, OECT > randomEdges( AdjacencyDiGraph< T, VCT, IECT, OECT >& ag, size_t number )
    {
	typedef AdjacencyDiGraph< T, VCT, IECT, OECT > G;
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

    
    
    GROUP_CTOR_DECL( AdjacencyDiGraphTest );

    void init();
};

#endif
