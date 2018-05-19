#include "adjacencyDiGraphTest.hpp"
#include "testMacros.hpp"

#include <algorithm>
#include <iterator>


#ifndef DEBUG
#define DEBUG false
#endif

std::default_random_engine AdjacencyDiGraphTest::mRandom = std::default_random_engine( std::random_device()() );

AdjacencyDiGraphTest::TEST_CTOR( AddValueTest, "add values and find them again" );

void AdjacencyDiGraphTest::AddValueTest::init()
{
    mGraph = G();
    addVertex( mGraph, 0 );
    mVal = 0;
}

void AdjacencyDiGraphTest::AddValueTest::iterate()
{
    mVal = mVdist( mRandom );
    addVertex( mGraph, mVal );
}

bool AdjacencyDiGraphTest::AddValueTest::check() const
{
    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
    typename DiGraphTraits< G >::VIterT iv = findVertex( mGraph, mVal );

    if( iv < vs.first || iv >= vs.second )
	return false;		
    return true;
}

AdjacencyDiGraphTest::TEST_CTOR( AddEdgeTest, "add edges and find source and target" );

void AdjacencyDiGraphTest::AddEdgeTest::init()
{
    D( std::cout << "init add edge test" << std::endl; );
    mGraph = G();
    addVertex( mGraph, mVdist( mRandom ) );
    iterate();
    D( std::cout << "finished init add edge test" << std::endl; );
}

void AdjacencyDiGraphTest::AddEdgeTest::iterate()
{
    D( std::cout << "iterate add edge test" << std::endl; );
    addVertex( mGraph, mVdist( mRandom ) );
    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
    std::uniform_int_distribution<> edist( 0, std::distance( vs.first, vs.second ) - 1 );
    typename DiGraphTraits< G >::VIterT isrc = vs.first, itrg = vs.first;

    if( isrc == vs.second || itrg == vs.second )
	std::cout << "bad iterators, even before advancing " << std::endl;

    D( std::cout << "add edge in add edge test " << std::endl; );
    size_t srcSteps = edist( mRandom ), trgSteps = edist( mRandom );
    std::advance( isrc, srcSteps );
    std::advance( itrg, trgSteps );

    if( isrc == vs.second || itrg == vs.second )
	std::cout << "bad iterators, though advanced by " << srcSteps << " and " << trgSteps << std::endl;
    
    // ++itrg;
    graph::addEdge( mGraph, *isrc, *itrg );
    mSrcVal = value( mGraph, *isrc );
    mTrgVal = value( mGraph, *itrg );

    if( mSrcVal == 0 && mTrgVal == 0 )
	throw std::logic_error( "bad iterators" );
    
    D( std::cout << "finished iteration in add edge test" << std::endl; );
}

bool AdjacencyDiGraphTest::AddEdgeTest::check() const
{
    D( std::cout << "check add edge test" << std::endl; );
    typename DiGraphTraits< G >::VIterT src = findVertex( mGraph, mSrcVal )
	, trg = findVertex( mGraph, mTrgVal );

    typename DiGraphTraits< G >::OutRangeT srcOut = outEdges( mGraph, *src );
    typename DiGraphTraits< G >::InRangeT trgIn = inEdges( mGraph, *trg );

    typename DiGraphTraits< G >::InIterT fromSrc = findEdgeFrom( mGraph, *src, *trg );
    typename DiGraphTraits< G >::OutIterT toTrg = findEdgeTo( mGraph, *src, *trg );

    if( toTrg == srcOut.second ||
	fromSrc == trgIn.second )
    {
	D( std::cout << "failed - edge could not be found" << std::endl; );
	return false;
    }
    else if( value( mGraph, source( mGraph, *fromSrc ) ) != mSrcVal ||
	     value( mGraph, target( mGraph, *fromSrc ) ) != mTrgVal )
    {
	D( std::cout << "failed - edge retrieved at source should be ("
	   << mSrcVal << ", " << mTrgVal << ") but is ("
	   << value( mGraph, source( mGraph, *fromSrc ) ) << ", " << value( mGraph, target( mGraph, *fromSrc ) ) << ")" << std::endl; );
	for( typename DiGraphTraits< G >::InIterT ie = trgIn.first; ie != trgIn.second; ++ie )
	    std::cout << value( mGraph, source( mGraph, *ie ) ) << " to " << value( mGraph, target( mGraph, *ie ) ) << std::endl;
	
	return false;
    }
    else if( value( mGraph, source( mGraph, *toTrg ) ) != mSrcVal ||
	     value( mGraph, target( mGraph, *toTrg ) ) != mTrgVal )
    {
	D( std::cout << "failed - edge retrieved at target should be ("
	   << mSrcVal << "," << mTrgVal << ") but is ("
	   << value( mGraph, source( mGraph, *toTrg ) ) << "," << value( mGraph, target( mGraph, *toTrg ) ) << ")" << std::endl; );
	return false;
    }
    else
    {
	D( std::cout << "check add edge test succeeded" << std::endl; );
	return true;
    }
}

AdjacencyDiGraphTest::TEST_CTOR( RemoveVertexTest, "verify removed vertices vanish" );

void AdjacencyDiGraphTest::RemoveVertexTest::init()
{
    D( std::cout << "init remove vertex test" << std::endl; );
    mGraph = G();

    for( uint cInitVs = 0; cInitVs < mTestSize; ++cInitVs )
    {
	int val = mVdist( mRandom );
	addVertex( mGraph, val );
	D( std::cout << "inserted vertex " << val << std::endl; );
    }
    
    iterate();
    D( std::cout << "finished init remove vertex test" << std::endl; );
}

void AdjacencyDiGraphTest::RemoveVertexTest::iterate()
{
    D( std::cout << "iterate remove vertex test" << std::endl; );
    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
    uint noVertices = std::distance( vs.first, vs.second );
    if( noVertices > 0 )
    {
	uint noJumps = (std::uniform_int_distribution<>( 0, noVertices - 1 ) )( mRandom );
	typename DiGraphTraits< G >::VIterT irem = vs.first;
	std::advance( irem, noJumps );
	mRemVal = value( mGraph, *irem );
	D( std::cout << "attempting to remove vertex storing " << mRemVal << std::endl; );
	removeVertex( mGraph, *irem );
    }
    else
	D( std::cout << "nothing to remove, graph is empty" << std::endl; );
    D( std::cout << "done iterate remove vertex test" << std::endl; );
}

bool AdjacencyDiGraphTest::RemoveVertexTest::check() const
{
    D( std::cout << "check remove vertex test" << std::endl; );
    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
    // std::cout << "got range" << std::endl;
    if( findVertex( mGraph, mRemVal ) != vs.second )
    {
	D( std::cout << "check failed" << std::endl; );
	return false;
    }
    D( std::cout << "check succeeded" << std::endl; );
    return true;
}

AdjacencyDiGraphTest::TEST_CTOR( RemoveVertexNonCorruptionTest, "verify removal keeps all other vertices" );

void AdjacencyDiGraphTest::RemoveVertexNonCorruptionTest::init()
{
    D( std::cout << "init remove vertex test" << std::endl; );
    mGraph = G();

    for( uint cInitVs = 0; cInitVs < mTestSize + 1; ++cInitVs )
    {
	int val = mVdist( mRandom );
	addVertex( mGraph, val );
	D( std::cout << "inserted vertex " << val << std::endl; );
    }
    
    iterate();
    D( std::cout << "finished init remove vertex test" << std::endl; );
}

void AdjacencyDiGraphTest::RemoveVertexNonCorruptionTest::iterate()
{
    D( std::cout << "iterate remove vertex test" << std::endl; );
    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
    uint noVertices = std::distance( vs.first, vs.second );
    if( noVertices > 0 )
    {
	uint noJumps = (std::uniform_int_distribution<>( 0, noVertices - 1 ) )( mRandom );
	typename DiGraphTraits< G >::VIterT irem = vs.first;
	std::advance( irem, noJumps );
	int remVal = value( mGraph, *irem );
	mRemain.clear();
	for( typename DiGraphTraits< G >::VIterT iv = vs.first; iv != vs.second; ++iv )
	{
	    int ival = value( mGraph, *iv );
	    if( ival != remVal )
		mRemain.push_back( ival );
	}
	removeVertex( mGraph, *irem );
	D( std::cout << "removed" << remVal << std::endl; );
    }
    else
	D( std::cout << "nothing to remove, graph is empty" << std::endl; );
}

bool AdjacencyDiGraphTest::RemoveVertexNonCorruptionTest::check() const
{
    D( std::cout << "check remove vertex test" << std::endl; );
    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );

    for( int remainVal : mRemain )
    {
	if( findVertex( mGraph, remainVal ) == vs.second )
	{
	    D( std::cout << "check failed, could not find " << remainVal << std::endl; );
	    return false;
	}
    }
    D( std::cout << "check succeeded" << std::endl; );
    return true;
}

AdjacencyDiGraphTest::TEST_CTOR( RemoveVertexIncomingTest, "removed vertices are not found as sources" );

void AdjacencyDiGraphTest::RemoveVertexIncomingTest::init()
{
    D( std::cout << "init remove vertex incoming test" << std::endl; );
    mGraph = G();
    std::uniform_int_distribution<> eventDist( 0, 4 );

    for( uint i = 0; i < mTestSize; ++i )
    {
	// add vertex of random value
	int valAdd = mVdist( mRandom );
	addVertex( mGraph, valAdd );
	// would be nice if addVertex could return iterator
	typename DiGraphTraits< G >::VIterT ivadd = findVertex( mGraph, valAdd );
	D( std::cout << "add value " << valAdd << std::endl; );
	// add edges to preexisting vertices randomly in both directions
	typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
	for( DiGraphTraits< G >::VIterT iv = vs.first; iv != vs.second; ++iv )
	{
	    int event = eventDist( mRandom );
	    if( event == 0 )
	    {
		addEdge( mGraph, *ivadd, *iv );
		D( std::cout << "added (" << value( mGraph, *ivadd ) << ", " << value( mGraph, *iv ) << ")" << std::endl; );
	    }
	    else if( event == 1 )
	    {
		addEdge( mGraph, *iv, *ivadd );
		D( std::cout << "added (" << value( mGraph, *iv ) << ", " << value( mGraph, *ivadd ) << std::endl; );
	    }
	}
    }

    // std::cout << "graph constructed" << std::endl << mGraph;
    iterate();
}

void AdjacencyDiGraphTest::RemoveVertexIncomingTest::iterate()
{
    D( std::cout << "iterate remove vertex incoming test" << std::endl; );

    // D( std::cout << "on graph " << std::endl << mGraph << std::endl; );

    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
    uint noVertices = std::distance( vs.first, vs.second );
    if( noVertices > 0 )
    {
	uint noJumps = (std::uniform_int_distribution<>( 0, noVertices - 1 ) )( mRandom );
	typename DiGraphTraits< G >::VIterT irem = vs.first;
	std::advance( irem, noJumps );
	mRemVal = value( mGraph, *irem );
	outEdges( mGraph, *irem );

	D( std::cout << "attempting to remove vertex storing " << mRemVal << std::endl; );
	removeVertex( mGraph, *irem );
    }

    // D( std::cout << "changed graph " << std::endl << mGraph << std::endl; );
}

bool AdjacencyDiGraphTest::RemoveVertexIncomingTest::check() const
{
    D( std::cout << "check remove vertex incoming test" << std::endl; );
    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
    for( typename DiGraphTraits< G >::VIterT iv = vs.first; iv != vs.second; ++iv )
    {
	typename DiGraphTraits< G >::InRangeT ins = inEdges( mGraph, *iv );
	for( typename DiGraphTraits< G >::InIterT iie = ins.first; iie != ins.second; ++iie )
	{
	    if( value( mGraph, source( mGraph, *iie ) ) == mRemVal )
	    {
		D( std::cout << "found edge from " << value( mGraph, source( mGraph, *iie ) )
		   << " to " << value( mGraph, target( mGraph, *iie ) ) << std::endl; );
		return false;
	    }
	}
    }
    return true;
}

AdjacencyDiGraphTest::TEST_CTOR( RemoveVertexOutgoingTest, "removed vertices are not found as targets" );

void AdjacencyDiGraphTest::RemoveVertexOutgoingTest::init()
{
    D( std::cout << "init remove vertex incoming test" << std::endl; );
    mGraph = G();
    std::uniform_int_distribution<> boolDist( 0, 1 );

    for( uint i = 0; i < mTestSize; ++i )
    {
	// add vertex of random value
	int valAdd = mVdist( mRandom );
	addVertex( mGraph, valAdd );
	// would be nice if addVertex could return iterator
	typename DiGraphTraits< G >::VIterT ivadd = findVertex( mGraph, valAdd );
	D( std::cout << "add value " << valAdd << std::endl; );
	// add edges to preexisting vertices randomly in both directions
	typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
	for( DiGraphTraits< G >::VIterT iv = vs.first; iv != vs.second; ++iv )
	{
	    if( boolDist( mRandom ) == 1 )
		addEdge( mGraph, *ivadd, *iv );
	    if( boolDist( mRandom ) == 1 )
		addEdge( mGraph, *iv, *ivadd );
	}
    }
    iterate();
}

void AdjacencyDiGraphTest::RemoveVertexOutgoingTest::iterate()
{
    D( std::cout << "iterate remove vertex incoming test" << std::endl; );

    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
    uint noVertices = std::distance( vs.first, vs.second );
    if( noVertices > 0 )
    {
	uint noJumps = (std::uniform_int_distribution<>( 0, noVertices - 1 ) )( mRandom );
	typename DiGraphTraits< G >::VIterT irem = vs.first;
	std::advance( irem, noJumps );
	mRemVal = value( mGraph, *irem );
	outEdges( mGraph, *irem );

	D( std::cout << "attempting to remove vertex storing " << mRemVal << std::endl; );
	removeVertex( mGraph, *irem );
    }
}

bool AdjacencyDiGraphTest::RemoveVertexOutgoingTest::check() const
{
    D( std::cout << "check remove vertex incoming test" << std::endl; );
    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
    for( typename DiGraphTraits< G >::VIterT iv = vs.first; iv != vs.second; ++iv )
    {
	typename DiGraphTraits< G >::InRangeT ins = inEdges( mGraph, *iv );
	for( typename DiGraphTraits< G >::InIterT iie = ins.first; iie != ins.second; ++iie )
	{
	    if( value( mGraph, target( mGraph, *iie ) ) == mRemVal )
	    {
		D( std::cout << "found edge from " << value( mGraph, source( mGraph, *iie ) )
		   << " to " << value( mGraph, target( mGraph, *iie ) ) << std::endl; );
		return false;
	    }
	}
    }
    return true;
}

AdjacencyDiGraphTest::TEST_CTOR( RemoveEdgeTest, "verify removed edges vanish" );

void AdjacencyDiGraphTest::RemoveEdgeTest::init()
{
    D( std::cout << "init: remove edges " << std::endl; );
    std::uniform_int_distribution<> vdist;
    mGraph = randomVertices< int, VecMap, InVec, InVec, std::uniform_int_distribution >( mTestSize, vdist );
    randomEdges( mGraph, mTestSize + 1 );
    mJumpDist = std::uniform_int_distribution<>( 0, mTestSize - 1 );
    iterate();
}

void AdjacencyDiGraphTest::RemoveEdgeTest::iterate()
{
    D( std::cout << "iterate: remove edges" << std::endl; );
    typename DiGraphTraits< G >::VRangeT vs = vertices( mGraph );
    typename G::VIterT isrc = vs.first + mJumpDist( mRandom )
	, itrg = vs.first + mJumpDist( mRandom );
    mSrcVal = value( mGraph, *isrc );
    mTrgVal = value( mGraph, *itrg );

    removeEdge( mGraph, *isrc, *itrg );
}

bool AdjacencyDiGraphTest::RemoveEdgeTest::check() const
{
    D( std::cout << "check: remove edges" << std::endl; );
    typename DiGraphTraits< G >::VIterT isrc = findVertex( mGraph, mSrcVal )
	, itrg = findVertex( mGraph, mTrgVal );
    
    if( findEdgeTo( mGraph, *isrc, *itrg ) != outEdges( mGraph, *isrc ).second ||
	findEdgeFrom( mGraph, *isrc, *itrg ) != inEdges( mGraph, *itrg ).second )
	return false;

    return true;
}

GROUP_CTOR( AdjacencyDiGraphTest, "adjacency directed graph tests" );

void AdjacencyDiGraphTest::init()
{
    std::shared_ptr< InterleaveRunner > interleave = std::shared_ptr< InterleaveRunner >( new InterleaveRunner() );
    std::shared_ptr< InterleaveRandomRunner > rinterleave = std::shared_ptr< InterleaveRandomRunner >( new InterleaveRandomRunner );

    uint simpleTestSize = 5 * mTestSize
	, advancedSize = 0.5 * mTestSize;
    
    addTest( new AddValueTest( simpleTestSize, mRepetitions ), interleave );
    addTest( new AddEdgeTest( simpleTestSize, mRepetitions ), interleave );
    addTest( new RemoveVertexTest( simpleTestSize, mRepetitions ), rinterleave );
    addTest( new RemoveVertexNonCorruptionTest( advancedSize, mRepetitions ), rinterleave );
    addTest( new RemoveVertexIncomingTest( advancedSize, mRepetitions ), rinterleave );
    addTest( new RemoveVertexOutgoingTest( advancedSize, mRepetitions ), rinterleave );
    addTest( new RemoveEdgeTest( advancedSize, mRepetitions ), rinterleave );
}
