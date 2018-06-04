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

RefinementTreeTest::ExpansionTest::ExpansionTest( uint testSize, uint reps)
    : ITest( "expansion", testSize, reps )
    , EXPANSION_SIZE( 2 )
{}

void RefinementTreeTest::ExpansionTest::init()
{
    D( std::cout << "expansion test init" << std::endl; );
    
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

void RefinementTreeTest::ExpansionTest::iterate()
{
    D( std::cout << "expansion test iterate " << std::endl; );
    // find leaf randomly
    std::vector< typename ExactRefinementTree::NodeT > ls = mpRtree->leaves(); // should select all leaves

    std::uniform_int_distribution<> pickDist( 0, ls.size() - 1 );
    typename std::vector< typename ExactRefinementTree::NodeT >::iterator iexp = ls.begin() + pickDist( mRandom );
    mPreviousNoNodes = mpRtree->tree().size();
    mPreviousHeight = tree::subtreeHeight( mpRtree->tree(), tree::root( mpRtree->tree() ) );
    mExpandNodeDepth = tree::depth( mpRtree->tree()
				    , static_cast< InsideGraphValue< typename ExactRefinementTree::RefinementT::NodeT >& >( *graph::value( mpRtree->leafMapping(), *iexp ) ).treeNode() );

    mpRtree->refine( *iexp, mRefiner );
}

bool RefinementTreeTest::ExpansionTest::check() const
{
    D( std::cout << "expansion test check" << std::endl; );
    size_t noNodesExpanded = mpRtree->tree().size();
    size_t newHeight = subtreeHeight( mpRtree->tree(), root( mpRtree->tree() ) );
    if( noNodesExpanded - mPreviousNoNodes != EXPANSION_SIZE )
    {
	D( std::cout << "numer of nodes test failed: previously " << mPreviousNoNodes << " now having " << noNodesExpanded << " expected the difference to be " << EXPANSION_SIZE << std::endl; );
	return false;
    }
    if( newHeight < mPreviousHeight )
    {
	D( std::cout << "depth test failed" << std::endl; );
    	return false;
    }
    // newHeight should be strictly higher if leaf node was expanded
    // expanded leaf node iff mPreviousHeight - 1 == mExpandNodeDepth
    if( newHeight <= mPreviousHeight && mPreviousHeight == mExpandNodeDepth )
    {
	D( std::cout << "stricter depth test failed - refined node at deepest level" << std::endl; );
	D( std::cout << "previous height " << mPreviousHeight << ", new height " << newHeight << " expanded node at depth " << mExpandNodeDepth << std::endl; );
	return false;
    }

    return true;
}

RefinementTreeTest::TEST_CTOR( LeavesTest, "number of leaves after expansion" );

void RefinementTreeTest::LeavesTest::init()
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

void RefinementTreeTest::LeavesTest::iterate()
{
    D( std::cout << "leaves test iterate" << std::endl; );
    LargestSideRefiner refiner;
    refineRandomLeaf( *mpRtree, refiner );
    ++mExpansionCounter;
}

bool RefinementTreeTest::LeavesTest::check() const
{
    D( std::cout << "leaves test check" << std::endl; );
    std::vector< typename ExactRefinementTree::NodeT > ls = mpRtree->leaves();
    if( ls.size() != 1 + (EXPANSION_SIZE - 1) * mExpansionCounter )
    {
	D( std::cout << "bad number of nodes " << ls.size() << " after " << mExpansionCounter << " expansions" << std::endl; );
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
    for( typename ExactRefinementTree::NodeT leaf : mpRtree->leaves( tree::root( mpRtree->tree() ) ) )
    {
	// leaf should have value
	const typename ExactRefinementTree::EnclosureT& leafBox = mpRtree->nodeValue( leaf ).value().get().getEnclosure();
	if( smallerBox.intersects( leafBox ) )
	{
	    typename std::vector< typename ExactRefinementTree::NodeT >::iterator iImg;
	    iImg = std::find_if( imageSmallerBox.begin(), imageSmallerBox.end()
				 , [&leaf, this] (const typename ExactRefinementTree::NodeT& n) {
				     std::optional< std::reference_wrapper< const InteriorTreeValue< typename ExactRefinementTree::EnclosureT > > > oN = mpRtree->nodeValue( n );
				     if( !oN )
					 return false;
				     return mpRtree->nodeValue( leaf ).value().get() == oN.value().get(); } );
	    // no equivalent box found in image
	    if( iImg == imageSmallerBox.end() )
	    {
		D( std::cout << "leaf box " << leafBox << " is not returned as image of " << smallerBox << " even though it is contained inside" << std::endl; );
		return false;
	    }
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
    
    // Ariadne::ConstraintSet constraints ( { -Ariadne::Real( rdist( mRandom ) ) <= x + y <= Ariadne::Real( rdist( mRandom ) )
    // 		, -Ariadne::Real( rdist( mRandom ) ) <= x <= Ariadne::Real( rdist( mRandom ) )
    // 		} );

    Ariadne::ConstraintSet constraints( { -Ariadne::Real( rdist( mRandom ) ) <= x <= Ariadne::Real( rdist( mRandom ) )
		, -Ariadne::Real( rdist( mRandom ) ) <= y <= Ariadne::Real( rdist( mRandom ) ) } );

    // std::cout << "constraints " << constraints << std::endl;
    
    std::function< Ariadne::ValidatedUpperKleenean( const typename ExactRefinementTree::EnclosureT&, const Ariadne::ConstraintSet& ) > intersect =
	[this] (auto& enc, auto& cs) { return !cs.separated( enc ).check( mpRtree->effort() ); };
	
    // for each leaf, check: either does not intersect smaller box or is contained in image
    std::vector< typename ExactRefinementTree::NodeT > constraintIntersection = mpRtree->intersection( constraints, intersect );
    for( typename ExactRefinementTree::NodeT leaf : mpRtree->leaves() )
    {
	// leaf should have value
	const typename ExactRefinementTree::EnclosureT& leafBox = mpRtree->nodeValue( leaf ).value().get().getEnclosure();
	if( possibly( intersect( leafBox, constraints ) ) )
	{
	    typename std::vector< typename ExactRefinementTree::NodeT >::iterator iImg;
	    iImg = std::find_if( constraintIntersection.begin(), constraintIntersection.end()
				 , [&leaf, this] (const typename ExactRefinementTree::NodeT& n) {
				     auto nval = mpRtree->nodeValue( n );
				     if( !nval )
					 return false;
				     return mpRtree->nodeValue( leaf ).value().get() == nval.value().get(); } );
	    // no equivalent box found in image
	    if( iImg == constraintIntersection.end() )
	    {
		D( std::cout << "leaf box " << leafBox << " is not returned as image of " << constraints << " even though it is contained inside" << std::endl; );
		return false;
	    }
	    constraintIntersection.erase( iImg );
	}
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

RefinementTreeTest::TEST_CTOR( NonLeafRemovalTest, "removes nodes from graph that are not leaves" );

void RefinementTreeTest::NonLeafRemovalTest::init()
{
    D( std::cout << "non leaf removal test init" << std::endl; );
    Ariadne::ExactBoxType safeBox( { {0,10}, {0,9} } );
    mpRtree.reset( new ExactRefinementTree( getDefaultTree( safeBox, 10 ) ) );
    iterate();
}

void RefinementTreeTest::NonLeafRemovalTest::iterate()
{
    D( std::cout << "non leaf removal test iterate" << std::endl; );
    typename ExactRefinementTree::NodeT refined = refineRandomLeaf( *mpRtree, mRefiner );
}

bool RefinementTreeTest::NonLeafRemovalTest::check() const
{
    D( std::cout << "non leaf removal test check" << std::endl; );
    std::stack< typename ExactRefinementTree::RefinementT::NodeT > toCheck;
    toCheck.push( tree::root( mpRtree->tree() ) );

    while( !toCheck.empty() )
    {
	typename ExactRefinementTree::RefinementT::NodeT& nex = toCheck.top();
	if( !tree::isLeaf( mpRtree->tree(), nex ) )
	{
	    typename ExactRefinementTree::MappingT::VIterT ivfound = graph::findVertex( mpRtree->leafMapping()
											, typename ExactRefinementTree::MappingT::ValueT( new InsideGraphValue< typename ExactRefinementTree::RefinementT::NodeT >( nex ) ) );
	    if( ivfound != vertices( mpRtree->leafMapping() ).second )
	    {
		D( std::cout << "box " << tree::value( mpRtree->tree(), nex )->getEnclosure()
		   << " is not a leaf but found in graph as "; ); 
		D( printNodeValue< ExactRefinementTree >( mpRtree->nodeValue( *ivfound ) ); );
		D( std::cout << std::endl; );
		return false;
	    }
	}

	toCheck.pop();
	typename tree::FixedBranchTreeTraits< typename ExactRefinementTree::RefinementT >::CRangeT cs = children( mpRtree->tree(), nex );
	for( ; cs.first != cs.second; ++cs.first )
	    toCheck.push( *cs.first );
	return true;
    }
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
    mRefined = refineRandomLeaf( *mpRtree, mRefiner );
}

// check for each child of refined node r':
//     for each leaf l: if l maps in to r' then l is in preimage of r'
bool RefinementTreeTest::PreimageTest::check() const
{
    D( std::cout << "preimgae test check" << std::endl; );
    auto allLeaves = mpRtree->leaves();
    for( typename ExactRefinementTree::NodeT& refinement : mpRtree->leaves( mRefined ) )
    {
	auto preimg = mpRtree->preimage( refinement );
	
	for( typename ExactRefinementTree::NodeT& leaf : allLeaves )
	{
	    Ariadne::ExactBoxType leafBox = mpRtree->nodeValue( leaf ).value().get().getEnclosure()
		, refinementBox = mpRtree->nodeValue( refinement ).value().get().getEnclosure();
	    auto iFound = std::find_if( preimg.begin(), preimg.end()
					, [this, &leaf] (const typename ExactRefinementTree::NodeT& pn) {
					    std::optional< std::reference_wrapper< const InteriorTreeValue< typename ExactRefinementTree::EnclosureT > > > opn = mpRtree->nodeValue( pn );
					    if( !opn )
						return false;
					    return mpRtree->nodeValue( leaf ).value().get() == opn.value().get(); } );
	    bool isReach = possibly( mpRtree->isReachable( mpRtree->nodeValue( leaf ).value(), refinement ) );
	    // false negative
	    if( isReach && iFound == preimg.end() )
	    {
		D( std::cout << "failed! " << leafBox << " maps into " << refinementBox << " but is not contained in preimage" << std::endl; );
		return false;
	    }
	    // false positive
	    if( !isReach && iFound != preimg.end() )
	    {
		D( std::cout << "failed! " << leafBox << " does not map into " << refinementBox << " but is contained in preimage" << std::endl; );
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
    D( std::cout << "preimgae test check" << std::endl; );
    auto allLeaves = mpRtree->leaves();
    for( typename ExactRefinementTree::NodeT& refinement : mpRtree->leaves( mRefined ) )
    {
	auto postimg = mpRtree->postimage( refinement );
	
	for( typename ExactRefinementTree::NodeT& leaf : allLeaves )
	{
	    Ariadne::ExactBoxType leafBox = mpRtree->nodeValue( leaf ).value().get().getEnclosure()
		, refinementBox = mpRtree->nodeValue( refinement ).value().get().getEnclosure();
	    auto iFound = std::find_if( postimg.begin(), postimg.end()
					, [this, &leaf] (const typename ExactRefinementTree::NodeT& pn) {
					    std::optional< std::reference_wrapper< const InteriorTreeValue< typename ExactRefinementTree::EnclosureT > > > opn = mpRtree->nodeValue( pn );
					    if( !opn )
						return false;
					    return mpRtree->nodeValue( leaf ).value().get() == opn.value().get(); } );
	    
	    bool isReach = possibly( mpRtree->isReachable( mpRtree->nodeValue( refinement ).value(), leaf ) );
	    // false negative
	    if( isReach && iFound == postimg.end() )
	    {
		D( std::cout << "failed! " << refinementBox << " maps into " << leafBox << " but is not contained in postimage" << std::endl; );
		return false;
	    }
	    // false positive
	    if( !isReach && iFound != postimg.end() )
	    {
		D( std::cout << "failed! " << refinementBox << " does not map into " << leafBox << " but is contained in postimage" << std::endl; );
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
	D( std::cout << "failure! always unsafe is not being mapped to" << std::endl; );
	D( graph::print( std::cout, mpRtree->leafMapping(), gvalCon ); );

	return false;
    }
    
    for( typename ExactRefinementTree::NodeT& n : mpRtree->leaves() )
    {
	bool  mapsToUnsafe = possibly( mpRtree->isReachable( mpRtree->nodeValue( n ).value().get(), mpRtree->outside() ) )
	    , inUnsafePreimg = std::any_of( unsafePreimg.begin(), unsafePreimg.end()
					    , [this, &n] (const typename ExactRefinementTree::NodeT& pn) {
						return mpRtree->nodeValue( pn ).value().get() == mpRtree->nodeValue( n ).value().get(); } );
	if( mapsToUnsafe != inUnsafePreimg )
	{
	    D( std::cout << "node "; );
	    D( printNodeValue< ExactRefinementTree >( mpRtree->nodeValue( n ) ); );
	    D( std::cout << ", maps to unsafe " << mapsToUnsafe << ", contained in unsafe preimage " << inUnsafePreimg << std::endl; );
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

    addTest( new ExpansionTest( 1 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new LeavesTest( 1 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new IntersectionTest( 1 * mTestSize, 0.1 * mRepetitions ), pRcontinuous );
    addTest( new CSetIntersectionTest( 1*mTestSize, 0.1 * mRepetitions ), pRcontinuous );
    addTest( new NonLeafRemovalTest( 1 * mTestSize, 0.1 * mRepetitions ), pRcontinuous );
    addTest( new PreimageTest( 1 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new PostimageTest( 1 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new AlwaysUnsafeTest( 1 * mTestSize, 0.1 * mRepetitions ), pRinterleave );    
}
