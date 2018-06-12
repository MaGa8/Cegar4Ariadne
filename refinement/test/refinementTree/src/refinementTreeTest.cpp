#include "refinementTreeTest.hpp"
#include "testMacros.hpp"

#include "geometry/interval.decl.hpp"
#include "geometry/box.decl.hpp"
#include "geometry/function_set.hpp"
#include "function/domain.hpp"
#include "function/constraint.hpp"

#include <algorithm>
#include <stack>
#include <cmath>

#ifndef DEBUG
#define DEBUG false
#endif


std::default_random_engine RefinementTreeTest::mRandom = std::default_random_engine( std::random_device()() );

RefinementTreeTest::TEST_CTOR( SizeTest, "number of states after expansion" );

void RefinementTreeTest::SizeTest::init()
{
    D( std::cout << "leaves test init" << std::endl; );
    Ariadne::BoundedConstraintSet safeSet( Ariadne::RealBox( { {0,10}, {0,10} } ) );

    Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::EffectiveScalarFunction constraintExpression = (a + b);
    // static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > vspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x, y} );
    mpRtree.reset( new ExactRefinementTree( safeSet, f, Ariadne::Effort( 5 ) ) );
    mExpansionCounter = 0;
    iterate();
}

void RefinementTreeTest::SizeTest::iterate()
{
    D( std::cout << "leaves test iterate" << std::endl; );
    refineRandomLeaf( *mpRtree, mRefiner );
    ++mExpansionCounter;
}

bool RefinementTreeTest::SizeTest::check() const
{
    D( std::cout << "leaves test check" << std::endl; );
    if( graph::size( mpRtree->graph() ) != 2 + (EXPANSION_SIZE - 1) * mExpansionCounter ) // outside and initial node 4 free
    {
	std::cout << "bad number of nodes " << graph::size( mpRtree->graph() ) << " after " << mExpansionCounter << " expansions" << std::endl;
	return false;
    }
    return true;
}

RefinementTreeTest::TEST_CTOR( IntersectionTest, "completeness of intersection with box" );

void RefinementTreeTest::IntersectionTest::init()
{
    Ariadne::BoundedConstraintSet safeSet( Ariadne::RealBox( { {0,10}, {0,10} } ) );

    Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::EffectiveScalarFunction constraintExpression = (a + b);
    // static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > vspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x, y} );
    mpRtree.reset( new ExactRefinementTree( safeSet, f, Ariadne::Effort( 5 ) ) );
    iterate();
}

void RefinementTreeTest::IntersectionTest::iterate()
{
    refineRandomLeaf( *mpRtree, mRefiner );
}

bool RefinementTreeTest::IntersectionTest::check() const
{
    // dimensions of root box
    Ariadne::ExactBoxType initAbs = mpRtree->initialEnclosure();
    Ariadne::Bounds< Ariadne::FloatDP > wid = initAbs[ 0 ].upper() - initAbs[ 0 ].lower()
	, hig = initAbs[ 1 ].upper() - initAbs[ 1 ].lower();
    // generate new random box
    std::uniform_real_distribution<> wdist( 0.0, wid.lower().get_d() ), hdist( 0.0, hig.lower().get_d() );
    double smallerWid = wdist( mRandom ), smallerHig = hdist( mRandom );
    Ariadne::ExactBoxType smallerBox( { { (wid.upper().get_d() - smallerWid) / 2, (wid.upper().get_d() + smallerWid) / 2 }
	    , { (hig.upper().get_d() - smallerHig) / 2, (hig.upper().get_d() + smallerHig) / 2 } } );

    // for each leaf, check: either does not intersect smaller box or is contained in image
    std::vector< typename ExactRefinementTree::NodeT > imageSmallerBox = mpRtree->intersection( smallerBox );
    auto vrange = graph::vertices( mpRtree->graph() );
    for( ; vrange.first != vrange.second; ++vrange.first )
    {
	// leaf should have value
	auto iImg = std::find_if( imageSmallerBox.begin(), imageSmallerBox.end()
				  , std::bind( &ExactRefinementTree::equal, *mpRtree, *vrange.first, std::placeholders::_1 ) );
	auto nval = mpRtree->nodeValue( *vrange.first );

	bool intersects = false;
	if( nval )
	    intersects = smallerBox.intersects( nval.value().get().getEnclosure() );
	else
	    intersects = possibly( !(intersection( smallerBox, mpRtree->initialEnclosure() ) == smallerBox) );

	if( intersects && iImg == imageSmallerBox.end() )
	{
	    printNodeValue( nval );
	    std::cout << " is not returned as intersection of " << smallerBox << " even though it is contained inside" << std::endl;
	    return false;
	}
	if( !intersects && iImg != imageSmallerBox.end() )
	{
	    printNodeValue( nval );
	    std::cout << " is returned as intersection with " << smallerBox << " even though it is not contained inside" << std::endl;
	    return false;
	}
    }
    return true;
}

RefinementTreeTest::TEST_CTOR( CSetIntersectionTest, "completeness of intersection with constraint set" );

void RefinementTreeTest::CSetIntersectionTest::init()
{
    Ariadne::BoundedConstraintSet safeSet( Ariadne::RealBox( { {0,10}, {0,10} } ) );

    // static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > vspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x, y} );
    mpRtree.reset( new ExactRefinementTree( safeSet, f, Ariadne::Effort( 5 ) ) );
    iterate();
}

void RefinementTreeTest::CSetIntersectionTest::iterate()
{
    refineRandomLeaf( *mpRtree, mRefiner );
}

bool RefinementTreeTest::CSetIntersectionTest::check() const
{
    Ariadne::ExactBoxType initAbs = mpRtree->initialEnclosure();
    // dimensions of root box
    Ariadne::Bounds< Ariadne::FloatDP > wid = initAbs[ 0 ].upper() - initAbs[ 0 ].lower()
	, hig = initAbs[ 1 ].upper() - initAbs[ 1 ].lower();

   // generate new random box
    std::uniform_real_distribution<> rdist( 0.0, std::min( wid.lower().get_d(), hig.lower().get_d() ) );
    Ariadne::EffectiveScalarFunction x = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, y = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );

    Ariadne::RealBox constraintBox = { { -rdist( mRandom ), rdist( mRandom ) }, { -rdist( mRandom ), rdist( mRandom ) } };
    Ariadne::BoundedConstraintSet constraints( constraintBox );

    // std::cout << "constraints " << constraints << std::endl;
    
    std::function< Ariadne::ValidatedUpperKleenean( const typename ExactRefinementTree::EnclosureT&, const Ariadne::BoundedConstraintSet& ) > intersect =
	[this] (auto& enc, auto& cs) { return !cs.separated( enc ).check( mpRtree->effort() ); };
	
    // for each leaf, check: either does not intersect smaller box or is contained in image
    std::vector< typename ExactRefinementTree::NodeT > constraintIntersection = mpRtree->intersection( constraints, intersect );

    auto vrange = graph::vertices( mpRtree->graph() );
    for( auto iv = vrange.first; iv != vrange.second; ++iv )
    {
	auto vval = mpRtree->nodeValue( *iv );
	if( vval )
	{
	    if( possibly( intersect( vval.value().get().getEnclosure(), constraints ) ) )
	    {
		typename std::vector< typename ExactRefinementTree::NodeT >::iterator iImg;
		iImg = std::find_if( constraintIntersection.begin(), constraintIntersection.end()
				     , [&vval, this] (const typename ExactRefinementTree::NodeT& n) {
					   auto nval = mpRtree->nodeValue( n );
					   if( !nval )
					       return false;
					   return vval.value().get() == nval.value().get(); } );
		// no equivalent box found in image
		if( iImg == constraintIntersection.end() )
		{
		    D( std::cout << "state " << vval.value().get() << " is not returned as image of " << constraints << " even though it is contained inside" << std::endl; );
		    return false;
		}
		constraintIntersection.erase( iImg );
	    }
	}
	// leave outside node be
    }
    if( constraintIntersection.size() > 1 ) // did not remove outside node
    {
	std::cout << "there are elements returned as intersection which could not be found as leaves" << std::endl;
	for( const typename ExactRefinementTree::NodeT& n : constraintIntersection )
	{
	    auto nval = mpRtree->nodeValue( n );
	    if( nval )
		std::cout << nval.value().get().getEnclosure();
	    else
		std::cout << "[outside] (okay, is not a leaf)";
	    std::cout << std::endl;
	}
	return false;
    }
    return true;
}

RefinementTreeTest::TEST_CTOR( RefinedNodesRemovalTest, "removes refined nodes from graph" );

void RefinementTreeTest::RefinedNodesRemovalTest::init()
{
    D( std::cout << "refined nodes removal test init" << std::endl; );
    Ariadne::ExactBoxType safeBox( { {0,10}, {0,9} } );
    mpRtree.reset( new ExactRefinementTree( getDefaultTree( safeBox, 10 ) ) );
    iterate();
}

void RefinementTreeTest::RefinedNodesRemovalTest::iterate()
{
    D( std::cout << "non leaf removal test iterate" << std::endl; );
    typename ExactRefinementTree::NodeT refined = refineRandomLeaf( *mpRtree, mRefiner );
}

bool RefinementTreeTest::RefinedNodesRemovalTest::check() const
{
    D( std::cout << "non leaf removal test check" << std::endl; );
    auto vrange = graph::vertices( mpRtree->graph() );
    for( auto iv = vrange.first; iv != vrange.second; ++iv )
    {
	auto vval = mpRtree->nodeValue( *iv );
	if( vval )
	{
	    auto iCover = std::find_if( vrange.first, vrange.second,
					[this, &vval] (const typename ExactRefinementTree::NodeT& n) {
					    auto nval = mpRtree->nodeValue( n );
					    if( !nval )
						return false;
					    if( nval.value().get() == vval.value().get() )
						return false;
					    Ariadne::BoundedConstraintSet vconstr( Ariadne::RealBox( vval.value().get().getEnclosure() ) );
					    return definitely( vconstr.covers( nval.value().get().getEnclosure() ) );
					} );
	    if( iCover != vrange.second )
	    {
		std::cout << "found " << mpRtree->nodeValue( *iCover ).value().get()
			  << " contained inside " << vval.value().get() << std::endl;
		return false;
	    }
	}
    }
    return true;
}

RefinementTreeTest::TEST_CTOR( PreimageTest, "preimage is complete and only complete" );

void RefinementTreeTest::PreimageTest::init()
{
    D( std::cout << "preimage test init" << std::endl; );
    Ariadne::BoundedConstraintSet safeSet( Ariadne::RealBox( { {0,1}, {0,2} } ) );
    // need refinement tree with non-static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > fspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( fspace, {x * x, y * y} );

    mpRtree.reset( new ExactRefinementTree( safeSet, f, Ariadne::Effort( 5 ) ) );
}

// expand
void RefinementTreeTest::PreimageTest::iterate()
{
    D( std::cout << "preimage test iterate" << std::endl; );

    GraphVertexPrintConverter conv( *mpRtree );
    mRefined = refineRandomLeaf( *mpRtree, mRefiner );
}

// check for each child of refined node r':
//     for each leaf l: if l maps in to r' then l is in preimage of r'
bool RefinementTreeTest::PreimageTest::check() const
{
    D( std::cout << "preimgae test check" << std::endl; );

    auto vrange = graph::vertices( mpRtree->graph() );
    for( auto iv = vrange.first; iv != vrange.second; ++iv )
    {
	auto preimg = mpRtree->preimage( *iv );

	for( auto iu = vrange.first; iu != vrange.second; ++iu )
	{
	    auto iFound = std::find_if( preimg.begin(), preimg.end()
					, std::bind( &ExactRefinementTree::equal, &*mpRtree, *iu, std::placeholders::_1 ) );

	    bool isReach = possibly( mpRtree->isReachable( *iu, *iv ) );
	    // false negative
	    if( (isReach && iFound == preimg.end() ) ||
		(!isReach && iFound != preimg.end() ) )
	    {
		std::cout << "failed! ";
		printNodeValue( mpRtree->nodeValue( *iu ) );
		std::cout << " maps into ";
		printNodeValue( mpRtree->nodeValue( *iv ) );
		std::cout << " but is " << (isReach ? "not" : "") << " contained in preimage" << std::endl;
		return false;
	    }
	}
    }
    return true;
}

RefinementTreeTest::TEST_CTOR( PostimageTest, "postimage is complete and only complete" );

void RefinementTreeTest::PostimageTest::init()
{
    D( std::cout << "preimage test init" << std::endl; );
    Ariadne::BoundedConstraintSet safeSet( Ariadne::RealBox( { {0,1}, {0,2} } ) );
    // need refinement tree with non-static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > fspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( fspace, {x * x, y * y} );

    mpRtree.reset( new ExactRefinementTree( safeSet, f, Ariadne::Effort( 5 ) ) );
}

// expand
void RefinementTreeTest::PostimageTest::iterate()
{
    D( std::cout << "preimage test iterate" << std::endl; );
    mRefined = refineRandomLeaf( *mpRtree, mRefiner );
}

// check for each child of refined node r':
//     for each leaf l: if r' maps in to l then l is in postimage of r'
bool RefinementTreeTest::PostimageTest::check() const
{
    D( std::cout << "postimage test check" << std::endl; );

    auto vrange = graph::vertices( mpRtree->graph() );
    for( auto iv = vrange.first; iv != vrange.second; ++iv )
    {
	auto postimg = mpRtree->postimage( *iv );

	for( auto iu = vrange.first; iu != vrange.second; ++iu )
	{
	    auto iFound = std::find_if( postimg.begin(), postimg.end()
					, std::bind( &ExactRefinementTree::equal, &*mpRtree, *iu, std::placeholders::_1 ) );

	    bool isReach = possibly( mpRtree->isReachable( *iv, *iu ) );
	    // false negative
	    if( (isReach && iFound == postimg.end() ) ||
		(!isReach && iFound != postimg.end() ) )
	    {
		std::cout << "failed! ";
		printNodeValue( mpRtree->nodeValue( *iu ) );
		std::cout << " is being mapped into from ";
		printNodeValue( mpRtree->nodeValue( *iv ) );
		std::cout << " but is " << (isReach ? "not" : "") << " contained in postimage" << std::endl;
		return false;
	    }
	}
    }
    return true;
}

RefinementTreeTest::TEST_CTOR( AlwaysUnsafeTest, "always unsafe node is being mapped to in unsafe, static system" );

void RefinementTreeTest::AlwaysUnsafeTest::init()
{
    D( std::cout << "always unsafe test init" << std::endl; );

    Ariadne::BoundedConstraintSet safeSet( Ariadne::RealBox( { {0,10}, {0,10} } ) );
    // static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > vspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x + 1, y + 1} );
    mpRtree.reset( new ExactRefinementTree( safeSet, f, Ariadne::Effort( 5 ) ) );
}

void RefinementTreeTest::AlwaysUnsafeTest::iterate()
{
    refineRandomLeaf( *mpRtree, mRefiner );
}

// some leaf node exists that maps to always unsafe, because the system is static and unsafe by construction
bool RefinementTreeTest::AlwaysUnsafeTest::check() const
{
    std::function< Ariadne::ExactBoxType( const typename ExactRefinementTree::MappingT::ValueT& ) > gvalCon = GraphVertexPrintConverter( *mpRtree );
    auto unsafePreimg = mpRtree->preimage( mpRtree->outside() );
    if( unsafePreimg.empty() )
    {
	std::cout << "failure! always unsafe is not being mapped to" << std::endl;
	graph::print( std::cout, mpRtree->graph(), gvalCon );

	return false;
    }

    const auto vrange = graph::vertices( mpRtree->graph() );
    for( auto iv = vrange.first; iv != vrange.second; ++iv )
    {
	bool mapsToUnsafe = possibly( mpRtree->isReachable( *iv, mpRtree->outside() ) )
	    , inUnsafePreimg = std::any_of( unsafePreimg.begin(), unsafePreimg.end()
					    , std::bind( &ExactRefinementTree::equal, &*mpRtree, *iv, std::placeholders::_1 ) );

	if( mapsToUnsafe != inUnsafePreimg )
	{
	    std::cout << "node ";
	    printNodeValue( mpRtree->nodeValue( *iv ) );
	    std::cout << ", maps to unsafe " << mapsToUnsafe << ", contained in unsafe preimage " << inUnsafePreimg << std::endl;
	    return false;
	}
    }
    return true;
}

RefinementTreeTest::GROUP_CTOR( RefinementTreeTest, "refinement tree" );

void RefinementTreeTest::init()
{
    std::shared_ptr< InterleaveRandomRunner > pRinterleave( new InterleaveRandomRunner() );
    std::shared_ptr< ContinuousRandomRunner > pRcontinuous( new ContinuousRandomRunner() );
    std::shared_ptr< OnlyOnceRunner > pOnce( new OnlyOnceRunner() );

    addTest( new SizeTest( 1 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new IntersectionTest( 1 * mTestSize, 0.1 * mRepetitions ), pRcontinuous );
    addTest( new CSetIntersectionTest( 1*mTestSize, 0.1 * mRepetitions ), pRcontinuous );
    addTest( new RefinedNodesRemovalTest( 1 * mTestSize, 0.1 * mRepetitions ), pRcontinuous );
    addTest( new PreimageTest( 1 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new PostimageTest( 1 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new AlwaysUnsafeTest( 1 * mTestSize, 0.1 * mRepetitions ), pRinterleave );    
}
