#include "cegarTest.hpp"
#include "testMacros.hpp"

#include <limits>

#ifndef DEBUG
#define DEBUG true
#endif

std::default_random_engine CegarTest::mRandom = std::default_random_engine( std::random_device()() );

std::function< Ariadne::ValidatedUpperKleenean( const typename CegarTest::ExactRefinementTree::EnclosureT&, const Ariadne::BoundedConstraintSet& ) > CegarTest::mIntersectConstraints =
    [] (auto& enc, auto& cset) {return !(cset.separated( enc ).check( Ariadne::Effort( 10 ) ) ); };

CegarTest::TEST_CTOR( FindCounterexampleTest, "finds counterexample" )

void CegarTest::FindCounterexampleTest::init()
{
    double boundary = 1000, delta = 1;

    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( {x, y}, {(1 + x)*(1 + x), (1 + y)*(1 + y)} );

    Ariadne::EffectiveScalarFunction cx = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, cy = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    mpInitial.reset( new Ariadne::BoundedConstraintSet( Ariadne::RealBox( { {-1, 1}, {-1, 1} } ) ) );
    
    Ariadne::RealBox initialAbs = { {-boundary - delta, boundary + delta}               // initial box wide enough so counterexamples inside exist
					 , {-boundary - delta, boundary + delta} };          // cs.constraint_bounds();

    mpRtree.reset( new ExactRefinementTree( Ariadne::BoundedConstraintSet( initialAbs ), f, Ariadne::Effort( 10 ) ) );
}

void CegarTest::FindCounterexampleTest::iterate()
{
    refineRandomLeaf( *mpRtree, mRefiner );
}

bool CegarTest::FindCounterexampleTest::check() const
{
    auto initialNodes = mpRtree->intersection( *mpInitial, mIntersectConstraints );

    CounterexampleStore< typename ExactRefinementTree::EnclosureT, RandomStateValue, GreatestState > store( mStateH, mCexH );
    findCounterexample( *mpRtree, initialNodes.begin(), initialNodes.end(), store );
    // system radially diverges
    if( !store.hasCounterexample() )
    {
	std::cout << "no counterexample obtainable" << std::endl;
	return false;
    }

    auto counterex = store.obtain();
    if( definitely( !mpRtree->overlapsConstraints( *mpInitial, counterex.first.front() ) ) )
    {
	std::cout << "counterexample does not begin in initial set" << std::endl;
	std::cout << "first node ";
	printNode( *mpRtree, counterex.first.front() );
	return false;
    }
    return true;
}


CegarTest::TEST_CTOR( FindNoCounterexampleTest, "not finds nonexistent counterexamples" )

void CegarTest::FindNoCounterexampleTest::init()
{
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::RealConstant b1( "b1", Ariadne::Real( 3.7 ) ), b2( "b2", Ariadne::Real( 3.8 ) );
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( {x, y}, { b1*x*(1 - x), b2*y*(1 - y)} );

    Ariadne::EffectiveScalarFunction cx = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, cy = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    mpInitial.reset( new Ariadne::BoundedConstraintSet( Ariadne::RealBox( { {-1, 1}, {-1, 1 } } ) ) );

    Ariadne::RealBox initialAbs = { {0, 1}               // initial box wide enough so counterexamples inside exist
					 , {0, 1} };          // cs.constraint_bounds();

    mpRtree.reset( new ExactRefinementTree( Ariadne::BoundedConstraintSet( initialAbs ), f, Ariadne::Effort( 10 ) ) );
}

void CegarTest::FindNoCounterexampleTest::iterate()
{
    refineRandomLeaf( *mpRtree, mRefiner );
}


bool CegarTest::FindNoCounterexampleTest::check() const
{
    auto initialNodes = mpRtree->intersection( *mpInitial, mIntersectConstraints );
    CounterexampleStore< typename ExactRefinementTree::EnclosureT, RandomStateValue, GreatestState > store( mStateH, mCexH );
    findCounterexample( *mpRtree, initialNodes.begin(), initialNodes.end(), store );
    // system radially diverges
    if( store.hasCounterexample() )
    {
	auto counterex = store.obtain();
	if( definitely( !mpRtree->isSafe( counterex.first.back() ) ) &&
	    definitely( !isSpurious( *mpRtree, counterex.first.begin(), counterex.first.end(), *mpInitial, Ariadne::Effort( 10 ) ) ) )
	{
	    printCounterexample( *mpRtree, counterex.first.begin(), counterex.first.end() );
	    std::cout << "is spurious " << isSpurious( *mpRtree, counterex.first.begin(), counterex.first.end(), *mpInitial, Ariadne::Effort( 10 ) ) << std::endl;
	    return false;
	}
    }
    return true;
}

CegarTest::InitialAbstraction::InitialAbstraction( uint size, uint repetitions )
    : ITest( "initial abstractions are complete and only complete", size, repetitions )
    , mTerm( size * mMaxNodesFactor )
{}

void CegarTest::InitialAbstraction::iterate()
{
    double wid = mInitialBoxLengthDist( mRandom )
	, hig = mInitialBoxLengthDist( mRandom )
	, delta = mDeltaDist( mRandom );

    mpRtree.reset( henonMap< typename ExactRefinementTree::EnclosureT >( wid + delta, hig + delta
									 , 1.4, 0.3, Ariadne::Effort( 10 ) ) );

    Ariadne::EffectiveScalarFunction cx = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	    , cy = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::RealConstant w = constant( "w", wid )
	, h = constant( "h", hig )
	, zero = constant( "zero", 0.0 );
    
    mpInitialSet.reset( new Ariadne::BoundedConstraintSet( { {0, wid}, {0, hig} },
							   { zero <= cx <= w, zero <= cy <= h } ) );
}

bool CegarTest::InitialAbstraction::check() const
{
    Ariadne::Effort effort( 10 );
    CheckInitialSet checkObs( *mpRtree, Ariadne::ConstraintSet( mpInitialSet->constraints() ), effort );
    KeepInitialSet< ExactRefinementTree > keeper;

    cegar( *mpRtree, *mpInitialSet, effort, mRefinement, mStateH, mCexH, mTerm, checkObs, keeper );
    if( checkObs.mpIssue )
    {
	std::cout << "found " << *checkObs.mpIssue << " in initial set, even though it is not in " << *mpInitialSet << std::endl;
	return false;
    }
    
    // verify that all final leaf nodes inside the initial set are contained
    ExactRefinementTree::NodeComparator ncomp( *mpRtree );
    for( auto leafRange = mpRtree->graph().vertices(); leafRange.first != leafRange.second; ++leafRange.first )
    {
	if( definitely( mpRtree->overlapsConstraints( *mpInitialSet, *leafRange.first ) ) )
	{
	    auto ifound = std::find_if( keeper.mNodes.begin(), keeper.mNodes.end()
					, std::bind( &ExactRefinementTree::NodeComparator::operator()
						     , &ncomp, *leafRange.first, std::placeholders::_1 ) );
	    if( ifound == keeper.mNodes.end() )
	    {
		std::cout << "node in initial set ";
		printNode( *mpRtree, *leafRange.first );
		std::cout << " is not contained in the initial set of the last iteration " << std::endl;
		return false;
	    }
	}
    }
    
    return true;
}

CegarTest::VerifySafety::VerifySafety( uint size, uint repetitions )
    : ITest( "verify safe by construction system", size, repetitions )
    , mTerm( size * mMaxNodesFactor )
{}

void CegarTest::VerifySafety::iterate()
{
    double wid = mInitialBoxLengthDist( mRandom )
	, hig = mInitialBoxLengthDist( mRandom )
	, delta = mDeltaDist( mRandom )
	, boundary = std::max( wid + delta, hig + delta )
	, attraction = mAttractionDist( mRandom );

    mpInitialSet.reset( new Ariadne::BoundedConstraintSet( Ariadne::RealBox( { {-wid, wid}, {-hig, hig} } ) ) );

    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::RealConstant a( "a", Ariadne::Real( attraction ) );
    Ariadne::EffectiveVectorFunction insideMap = Ariadne::make_function( {x, y}, {a*x, a*y} );

    // figure out how to use bounds on safe region to do this
    Ariadne::RealBox initialAbs = { {-boundary, boundary}, {-boundary, boundary} };

    mpRtree.reset( new ExactRefinementTree( Ariadne::BoundedConstraintSet( initialAbs ), insideMap, Ariadne::Effort( 10 ) ) );
}

bool CegarTest::VerifySafety::check() const
{
    auto cegarResult = cegar( *mpRtree, *mpInitialSet, Ariadne::Effort( 10 ), mRefinement, mStateH, mCexH, mTerm
			      // , iniPrint, cexPrint
			      );
    if( !definitely( cegarResult.first ) )
    {
	std::cout << "could not verify safety of system safe by construction, safety " << cegarResult.first
		  << " with " << graph::size( mpRtree->graph() ) << " nodes "
		  << std::endl;
	return false;
    }
    return true;
}

CegarTest::VerifyCounterexamples::VerifyCounterexamples( uint size, uint reps )
    : ITest( "verify consecutive nodes in counterexample are reachable", size, reps )
    , mTerm( mMaxNodesFactor * size )
{}

void CegarTest::VerifyCounterexamples::iterate()
{
    D( std::cout << "verify counterexamples iteration" << std::endl; );
    double wi = mInitialBoxLengthDist( mRandom )
	, hi = mInitialBoxLengthDist( mRandom )
	, ws = mSafeBoxLengthDist( mRandom )
	, hs = mSafeBoxLengthDist( mRandom );
    
    mpRtree.reset( henonMap< typename ExactRefinementTree::EnclosureT >( ws, hs, 1.4, 0.3, Ariadne::Effort( 10 ) ) );
    mpInitialSet.reset( new Ariadne::BoundedConstraintSet( { {0, wi}, {0, hi} } ) );
}

bool CegarTest::VerifyCounterexamples::check() const
{
    CounterexampleVerifier verifier;
    cegar( *mpRtree, *mpInitialSet, Ariadne::Effort( 10 ), mRefinement, mStateH, mCexH, mTerm, verifier );

    if( !verifier.mBadCounterexample.empty() )
    {
	std::cout << "found counterexample with bad link " << std::endl;
	printCounterexample( *mpRtree, verifier.mBadCounterexample.begin(), verifier.mBadCounterexample.end() );
	std::cout << "bad link being ";
	printNode( *mpRtree, *(verifier.mBadLink - 1) );
	std::cout << " -> ";
	printNode( *mpRtree, *verifier.mBadLink );
	std::cout << std::endl;
	return false;
    }
    return true;
}

CegarTest::LoopTest::LoopTest( uint size, uint reps )
    : ITest( "whole cegar loop", size, reps )
    , mTerm( mMaxNodesFactor * size )
{}

void CegarTest::LoopTest::iterate()
{
    double wi = mInitialBoxLengthDist( mRandom )
	, hi = mInitialBoxLengthDist( mRandom )
	, ws = mSafeSetBoxLengthDist( mRandom )
	, hs = mSafeSetBoxLengthDist( mRandom );
    
    mpRtree.reset( henonMap< typename ExactRefinementTree::EnclosureT >( ws, hs, 1.4, 0.3, Ariadne::Effort( 10 ) ) );
    mpInitialSet.reset( new Ariadne::BoundedConstraintSet( { {0, wi}, {0, hi} } ) );
}

bool CegarTest::LoopTest::check() const
{
    Ariadne::UpperBoxType initialBb = mpInitialSet->bounding_box();
    std::uniform_real_distribution<> widthDist = std::uniform_real_distribution<>( initialBb[ 0 ].lower().get_d()
										   , initialBb[ 0 ].upper().get_d() )
	, heightDist = std::uniform_real_distribution<>( initialBb[ 1 ].lower().get_d()
							 , initialBb[ 1 ].upper().get_d() );

    auto unsafeTrajectory = findUnsafeTrajectory( widthDist, heightDist, mpRtree->dynamics(), mpRtree->constraints(), mTestSize, mTestSize );

    auto cegarCounterex = cegar( *mpRtree, *mpInitialSet, Ariadne::Effort( 10 ), mRefinement, mStateH, mCexH, mTerm );

    if( !unsafeTrajectory.empty() && definitely( cegarCounterex.first ) )
    {
	std::cout << "found ";
	for_each( unsafeTrajectory.begin(), unsafeTrajectory.end()
		  , [] (const Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > >& p) {std::cout << p << " -> "; } );
	std::cout << std::endl << "but cegar found nothing" << std::endl << std::endl;
	std::cout << "initial set " << *mpInitialSet << std::endl;
	std::cout << "safe set " << mpRtree->constraints() << std::endl << std::endl;
	std::cout << "graph " << std::endl;
	std::function< Ariadne::ExactBoxType( const typename ExactRefinementTree::MappingT::ValueT& ) > printConv =
	    std::bind( &CegarTest::extractEnclosure< Ariadne::ExactBoxType >, *mpRtree, std::placeholders::_1 );

	graph::print( std::cout, mpRtree->graph(), printConv );

	for( auto vrange = graph::vertices( mpRtree->graph() ); vrange.first != vrange.second; ++vrange.first )
	{
	    auto nval = mpRtree->nodeValue( *vrange.first );
	    if( nval )
	    {
		auto mapped = Ariadne::image( nval.value().get().getEnclosure(), mpRtree->dynamics() );
		printNode( *mpRtree, *vrange.first ); std::cout << " -> " << mapped << std::endl;
	    }
	}
	
	return false;
    }
    else if( unsafeTrajectory.empty() )
    {
	if( definitely( !cegarCounterex.first ) )
	{                                                                 // because mapping random points comes with no guarantees
	    std::cout << "cegar found ";
	    printCounterexample( *mpRtree, cegarCounterex.second.begin(), cegarCounterex.second.end() ); // does endl
	    std::cout << "but could not find a point that maps there, " << std::endl
		      << "maybe the number of points sampled should be increased?" << std::endl;
	}
	else if( !definitely( cegarCounterex.first ) ) // case: indeterminate
	{
	    std::cout << "no points sampled left the safe set but cegar could not verify safety" << std::endl
		      << "maybe the maximum number of nodes of cegar should be increased?" << std::endl;
	}
	std::cout << std::endl;
    }

    return true;
}

CegarTest::GROUP_CTOR( CegarTest, "top level cegar functions" )

void CegarTest::init()
{
    std::shared_ptr< ITestRunner > pRinterleave( new InterleaveRandomRunner() )
	, pStateless( new StatelessRunner() );
    addTest( new FindCounterexampleTest( 0.5 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new FindNoCounterexampleTest( 0.5 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new InitialAbstraction( mTestSize, 0.1 * mRepetitions ), pStateless );
    addTest( new VerifySafety( mTestSize, 0.1 * mRepetitions ), pStateless );
    addTest( new VerifyCounterexamples( mTestSize, 0.05 * mRepetitions ), pStateless );
    addTest( new LoopTest( mTestSize, 0.05 * mRepetitions ), pStateless );
}
