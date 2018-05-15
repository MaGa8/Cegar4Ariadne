#include "refinementTreeTest.hpp"
#include "testMacros.hpp"
#include "refinementStrategy.hpp"

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
    , mpRefiner( new LargestSideRefiner< Ariadne::ExactIntervalType >() )
    , EXPANSION_SIZE( 2 )
{}

void RefinementTreeTest::ExpansionTest::init()
{
    D( std::cout << "expansion test init" << std::endl; );
    mpBox.reset( new Ariadne::ExactBoxType( { {0,10}, {0,10} } ) );

    Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::EffectiveScalarFunction constraintExpression = (a + b);
    // only upper bound on + values
    Ariadne::EffectiveConstraint c1 = (constraintExpression <= 10 );
    Ariadne::EffectiveConstraintSet cs = { c1 };
    // static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > vspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x, y} );
    mpRtree.reset( new ExactRefinementTree( *mpBox, cs, f, Ariadne::Effort( 5 ) ) );
    iterate();
}

void RefinementTreeTest::ExpansionTest::iterate()
{
    D( std::cout << "expansion test iterate " << std::endl; );
    // find leaf randomly
    std::vector< typename ExactRefinementTree::NodeT > ls = mpRtree->image( *mpBox ); // should select all leaves

    std::uniform_int_distribution<> pickDist( 0, ls.size() - 1 );
    typename std::vector< typename ExactRefinementTree::NodeT >::iterator iexp = ls.begin() + pickDist( mRandom );
    mPreviousNoNodes = mpRtree->tree().size();
    mPreviousHeight = tree::subtreeHeight( mpRtree->tree(), tree::root( mpRtree->tree() ) );
    mExpandNodeDepth = tree::depth( mpRtree->tree()
				    , static_cast< typename ExactRefinementTree::InsideGraphValue& >( *graph::value( mpRtree->leafMapping(), *iexp ) ).treeNode() );

    mpRtree->refine( *iexp, *mpRefiner );
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
    mpBox.reset( new Ariadne::ExactBoxType( { {0,10}, {0,10} } ) );

    Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::EffectiveScalarFunction constraintExpression = (a + b);
    // only upper bound on + values
    Ariadne::EffectiveConstraint c1 = (constraintExpression <= 10 );
    Ariadne::EffectiveConstraintSet cs = { c1 };
    // static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > vspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x, y} );
    mpRtree.reset( new ExactRefinementTree( *mpBox, cs, f, Ariadne::Effort( 5 ) ) );
    mExpansionCounter = 0;
    iterate();
}

void RefinementTreeTest::LeavesTest::iterate()
{
    D( std::cout << "leaves test iterate" << std::endl; );
    LargestSideRefiner< Ariadne::ExactIntervalType > refiner;
    refineRandomLeaf( *mpRtree, refiner );
    ++mExpansionCounter;
}

bool RefinementTreeTest::LeavesTest::check() const
{
    D( std::cout << "leaves test check" << std::endl; );
    std::vector< typename ExactRefinementTree::NodeT > ls = mpRtree->image( *mpBox );
    if( ls.size() != 1 + (EXPANSION_SIZE - 1) * mExpansionCounter )
    {
	D( std::cout << "bad number of nodes " << ls.size() << " after " << mExpansionCounter << " expansions" << std::endl; );
	return false;
    }
    return true;
}

RefinementTreeTest::TEST_CTOR( ImageTest, "completeness of components of image" );

void RefinementTreeTest::ImageTest::init()
{
    mpRootBox.reset( new Ariadne::ExactBoxType( { {0,10}, {0,10} } ) );

    Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::EffectiveScalarFunction constraintExpression = (a + b);
    // only upper bound on + values
    Ariadne::EffectiveConstraint c1 = (constraintExpression <= 10 );
    Ariadne::EffectiveConstraintSet cs = { c1 };
    // static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > vspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x, y} );
    mpRtree.reset( new ExactRefinementTree( *mpRootBox, cs, f, Ariadne::Effort( 5 ) ) );
    iterate();
}

void RefinementTreeTest::ImageTest::iterate()
{
    LargestSideRefiner< typename ExactRefinementTree::IntervalT > refiner;
    refineRandomLeaf( *mpRtree, refiner );
}

bool RefinementTreeTest::ImageTest::check() const
{
    // dimensions of root box
    Ariadne::Bounds< Ariadne::FloatDP > wid = (*mpRootBox)[ 0 ].upper() - (*mpRootBox)[ 0 ].lower()
	, hig = (*mpRootBox)[ 1 ].upper() - (*mpRootBox)[ 1 ].lower();
    // generate new random box
    std::uniform_real_distribution<> wdist( 0.0, wid.lower().get_d() ), hdist( 0.0, hig.lower().get_d() );
    double smallerWid = wdist( mRandom ), smallerHig = hdist( mRandom );
    Ariadne::ExactBoxType smallerBox( { { (wid.upper().get_d() - smallerWid) / 2, (wid.upper().get_d() + smallerWid) / 2 }
	    , { (hig.upper().get_d() - smallerHig) / 2, (hig.upper().get_d() + smallerHig) / 2 } } );

    // for each leaf, check: either does not intersect smaller box or is contained in image
    std::vector< typename ExactRefinementTree::NodeT > imageSmallerBox = mpRtree->image( smallerBox );
    for( typename ExactRefinementTree::NodeT leaf : mpRtree->leaves( tree::root( mpRtree->tree() ) ) )
    {
	// leaf should have value
	const typename ExactRefinementTree::EnclosureT& leafBox = mpRtree->nodeValue( leaf ).value().get().getEnclosure();
	if( smallerBox.intersects( leafBox ) )
	{
	    typename std::vector< typename ExactRefinementTree::NodeT >::iterator iImg;
	    iImg = std::find_if( imageSmallerBox.begin(), imageSmallerBox.end()
				 , [&leaf, this] (const typename ExactRefinementTree::NodeT& n) {
				     return mpRtree->nodeValue( leaf ).value().get() == mpRtree->nodeValue( n ).value().get(); } );
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

RefinementTreeTest::TEST_CTOR( NonLeafRemovalTest, "removes nodes from graph that are not leaves" );

void RefinementTreeTest::NonLeafRemovalTest::init()
{
    D( std::cout << "non leaf removal test init" << std::endl; );
    mpRootBox.reset( new Ariadne::ExactBoxType( { {0,10}, {0,9} } ) );
    mpRtree.reset( new ExactRefinementTree( getDefaultTree( *mpRootBox, 10 ) ) );
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
											, typename ExactRefinementTree::MappingT::ValueT( new ExactRefinementTree::InsideGraphValue( nex ) ) );
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
    mpRootBox.reset( new Ariadne::ExactBoxType( { {0,1}, {0,2} } ) );
    // need refinement tree with non-static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > fspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( fspace, {x * x, y * y} );

    Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::EffectiveConstraint c1 = ( (a + b) <= 5 );
    Ariadne::EffectiveConstraintSet cs = {c1};

    mpRtree.reset( new ExactRefinementTree( *mpRootBox, cs, f, Ariadne::Effort( 5 ) ) );
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
					    std::optional< std::reference_wrapper< const typename ExactRefinementTree::InteriorTreeValue > > opn = mpRtree->nodeValue( pn );
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
    mpRootBox.reset( new Ariadne::ExactBoxType( { {0,1}, {0,2} } ) );
    // need refinement tree with non-static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > fspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( fspace, {x * x, y * y} );

    Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::EffectiveConstraint c1 = ( (a + b) <= 5 );
    Ariadne::EffectiveConstraintSet cs = {c1};

    mpRtree.reset( new ExactRefinementTree( *mpRootBox, cs, f, Ariadne::Effort( 5 ) ) );
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
					    std::optional< std::reference_wrapper< const typename ExactRefinementTree::InteriorTreeValue > > opn = mpRtree->nodeValue( pn );
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

    typename ExactRefinementTree::EnclosureT rtAbs = Ariadne::ExactBoxType( { {0,10}, {0,10} } );
    Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::EffectiveScalarFunction constraintExpression = (a + b);
    // only upper bound on + values
    Ariadne::EffectiveConstraint c1 = (constraintExpression <= 9 );
    Ariadne::EffectiveConstraintSet cs = { c1 };
    // static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > vspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x, y} );
    mpRtree.reset( new ExactRefinementTree( rtAbs, cs, f, Ariadne::Effort( 5 ) ) );
}

void RefinementTreeTest::AlwaysUnsafeTest::iterate()
{
    refineRandomLeaf( *mpRtree, mRefiner );
}

// some leaf node exists that maps to always unsafe, because the system is static and unsafe by construction
bool RefinementTreeTest::AlwaysUnsafeTest::check() const
{
    auto unsafePreimg = mpRtree->preimage( mpRtree->alwaysUnsafe() );
    if( unsafePreimg.empty() )
    {
	D( std::cout << "failure! always unsafe is not being mapped to" << std::endl; );
	return false;
    }
    
    for( typename ExactRefinementTree::NodeT& n : mpRtree->leaves() )
    {
	bool unsafe = definitely( !mpRtree->isSafe( n ) )
	    , mapsToUnsafe = possibly( mpRtree->isReachable( mpRtree->nodeValue( n ).value().get(), mpRtree->alwaysUnsafe() ) )
	    , inUnsafePreimg = std::any_of( unsafePreimg.begin(), unsafePreimg.end()
					    , [this, &n] (const typename ExactRefinementTree::NodeT& pn) {
						return mpRtree->nodeValue( pn ).value().get() == mpRtree->nodeValue( n ).value().get(); } );
	std::array< bool, 3 > conds = { unsafe, mapsToUnsafe, inUnsafePreimg };
	if( !std::all_of( conds.begin(), conds.end(), [] (const bool b) { return b; } )
	    && !std::all_of( conds.begin(), conds.end(), [] (const bool b) { return !b; } ) )
	{
	    D( std::cout << "node "; );
	    D( printNodeValue< ExactRefinementTree >( mpRtree->nodeValue( n ) ); );
	    D( std::cout << " is unsafe " << unsafe << ", maps to unsafe " << mapsToUnsafe << ", contained in unsafe preimage " << inUnsafePreimg << std::endl; );
	    return false;
	}
	
    }
}

RefinementTreeTest::TEST_CTOR( PositiveCounterexampleTest, "test whether counterexample can be found" );

void RefinementTreeTest::PositiveCounterexampleTest::iterate() {}

bool RefinementTreeTest::PositiveCounterexampleTest::check() const
{
    /*
      l = nonSafetyLevel, c = constraintBoundary
      need: 1 - 1 / 2^x + (1 + c) - (1 + c) / 2^x > l
      => - (2 + c) / 2^x > l - 2 - c
      => - (2 + c) / (l - 2 - c) > 2^x
      => x > log2( (2 + c) / (2 + c - l) )
      in one dimension, so assuming equal number of splits along each dimension, multiply by 2
    */
    const double constraintBoundary = 2;
    const uint refinementDepth = 8;
    Ariadne::Effort effort( 5 );
    
    D( std::cout << "preimage test init" << std::endl; );
    Ariadne::ExactBoxType rootBox( { {0,1}, {0,3} } )
	, initialBox( { {0,0.75}, {0,1.25} } );
    // need refinement tree with non-static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > fspace = {x, y};
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( fspace, {x * x, y * y} );

    Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::EffectiveConstraint c1 = ( (a + b) <= constraintBoundary );
    Ariadne::EffectiveConstraintSet cs = {c1};

    ExactRefinementTree rtree( rootBox, cs, f, effort );

    // refine range (0,1.1) to less than 0.1, so 1.1 / 2^x <= 0.1
    std::cout << "need to refine " << refinementDepth << " levels " << std::endl;
    refineEqualDepth( rtree, refinementDepth );

    // root box is initial state
    std::vector< typename ExactRefinementTree::NodeT > initImage = rtree.image( initialBox );
    auto cexPath = findCounterexample( rtree, initImage.begin(), initImage.end() );

    std::cout << "is counterexample spurious? "
	      << isSpurious( rtree, cexPath.begin(), cexPath.end(), initImage.begin(), initImage.end(), effort )
	      << std::endl;

    if( cexPath.empty() )
	return false;
    
    return true;
}


RefinementTreeTest::GROUP_CTOR( RefinementTreeTest, "refinement tree" );

void RefinementTreeTest::init()
{
    std::shared_ptr< InterleaveRandomRunner > pRinterleave( new InterleaveRandomRunner() );
    std::shared_ptr< ContinuousRandomRunner > pRcontinuous( new ContinuousRandomRunner() );
    std::shared_ptr< OnlyOnceRunner > pOnce( new OnlyOnceRunner() );

    addTest( new ExpansionTest( mTestSize, mRepetitions ), pRinterleave );
    addTest( new LeavesTest( mTestSize, mRepetitions ), pRinterleave );
    addTest( new ImageTest( mTestSize, mRepetitions ), pRcontinuous );
    addTest( new NonLeafRemovalTest( mTestSize, mRepetitions ), pRcontinuous );
    addTest( new PreimageTest( mTestSize, mRepetitions ), pRinterleave );
    addTest( new PostimageTest( mTestSize, mRepetitions ), pRinterleave );
    addTest( new AlwaysUnsafeTest( mTestSize, mRepetitions ), pRinterleave );
    addTest( new PositiveCounterexampleTest( mTestSize, mRepetitions ), pOnce );
}
