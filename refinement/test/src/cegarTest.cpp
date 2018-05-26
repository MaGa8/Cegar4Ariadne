#include "cegarTest.hpp"
#include "testMacros.hpp"

#include "expression/space.hpp"
#include "expression/expression.hpp"
#include "function/function.hpp"

#include <limits>

#ifndef DEBUG
#define DEBUG false
#endif

std::default_random_engine CegarTest::mRandom = std::default_random_engine( std::random_device()() );

std::function< Ariadne::ValidatedLowerKleenean( const typename CegarTest::ExactRefinementTree::EnclosureT&, const Ariadne::ConstraintSet& ) > CegarTest::mIntersectConstraints =
    [] (auto& enc, auto& cset ) {return cset.overlaps( enc ).check( Ariadne::Effort( 10 ) ); };

CegarTest::TEST_CTOR( FindCounterexampleTest, "finds counterexample" )

void CegarTest::FindCounterexampleTest::init()
{
    double boundary = 1000, delta = 1;

    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( {x, y}, {(1 + x)*(1 + x), (1 + y)*(1 + y)} );

    Ariadne::EffectiveScalarFunction cx = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, cy = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::ConstraintSet cs = { cx*cx + cy*cy <= boundary };          // large value so it will take some time to diverge enough
    
    mpInitial.reset( new Ariadne::ConstraintSet( { -1 <= cx <= 1, -1 <= cy <= 1 } ) );

    Ariadne::ExactBoxType initialAbs = { {-boundary - delta, boundary + delta}               // initial box wide enough so counterexamples inside exist
					 , {-boundary - delta, boundary + delta} };          // cs.constraint_bounds();

    mpRtree.reset( new ExactRefinementTree( initialAbs, cs, f, Ariadne::Effort( 10 ) ) );
}

void CegarTest::FindCounterexampleTest::iterate()
{
    refineRandomLeaf( *mpRtree, mRefiner );
}

bool CegarTest::FindCounterexampleTest::check() const
{
    auto initialNodes = mpRtree->intersection( *mpInitial, mIntersectConstraints );
    
    auto counterex = findCounterexample( *mpRtree, initialNodes.begin(), initialNodes.end() );
    // system radially diverges
    if( graph::value( mpRtree->leafMapping(), counterex.back() )->isInside() )
	D( std::cout << "inside" << std::endl; );
    else
	D( std::cout << "outside" << std::endl; );
    if( counterex.empty() )
	return false;

    if( definitely( !mpRtree->relates( counterex.front(), *mpInitial, mIntersectConstraints ) ) )
    {
	D( std::cout << "counterexample does not begin in initial set" << std::endl; );
	return false;
    }
    return true;
}


CegarTest::TEST_CTOR( FindNoCounterexampleTest, "not finds nonexistent counterexamples" )

void CegarTest::FindNoCounterexampleTest::init()
{
    double boundary = 1.005, delta = -0.006; // cheating a bit, so initial abstraction is within constraint boundary

Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::RealConstant b1( "b1", Ariadne::Real( 3.7 ) ), b2( "b2", Ariadne::Real( 3.8 ) );
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( {x, y}, { b1*x*(1 - x), b2*y*(1 - y)} );

    Ariadne::EffectiveScalarFunction cx = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, cy = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
    Ariadne::ConstraintSet cs = { cx*cx + cy*cy <= boundary };          // large value so it will take some time to diverge enough

mpInitial.reset( new Ariadne::ConstraintSet( { -1 <= cx <= 1, -1 <= cy <= 1 } ) );

    Ariadne::ExactBoxType initialAbs = { {0, 1}               // initial box wide enough so counterexamples inside exist
					 , {0, 1} };          // cs.constraint_bounds();

    mpRtree.reset( new ExactRefinementTree( initialAbs, cs, f, Ariadne::Effort( 10 ) ) );
}

void CegarTest::FindNoCounterexampleTest::iterate()
{
    refineRandomLeaf( *mpRtree, mRefiner );
}


bool CegarTest::FindNoCounterexampleTest::check() const
{
    auto initialNodes = mpRtree->intersection( *mpInitial, mIntersectConstraints );
    auto counterex  = findCounterexample( *mpRtree, initialNodes.begin(), initialNodes.end() );
    // system radially diverges
    if( !counterex.empty() &&
	definitely( !mpRtree->isSafe( counterex.back() ) ) &&
	definitely( !isSpurious( *mpRtree, counterex.begin(), counterex.end(), *mpInitial, Ariadne::Effort( 10 ) ) ) )
    {
	printCounterexample( *mpRtree, counterex.begin(), counterex.end() );
	std::cout << "is spurious " << isSpurious( *mpRtree, counterex.begin(), counterex.end(), *mpInitial, Ariadne::Effort( 10 ) ) << std::endl;
	return false;
    }
    return true;
}

CegarTest::TEST_CTOR( LoopTest, "whole cegar loop" )

void CegarTest::LoopTest::iterate()
{
    double wid = mInitialBoxLengthDist( mRandom ),
	hig = mInitialBoxLengthDist( mRandom ),
	delta = mDeltaDist( mRandom )
	, boundary = std::max( wid + delta, hig + delta );

    std::cout << "initial " << wid << " x " << hig << " within square bound " << boundary << std::endl;

    Ariadne::RealConstant w( "w", Ariadne::Real( wid ) )
	, h( "h", Ariadne::Real( hig ) ); // must be the same as wid
    Ariadne::EffectiveScalarFunction cx = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, cy = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );

    mpInitialSet.reset( new Ariadne::ConstraintSet( { cx >= -w, cx <= w, cy >= -h, cy <= h } ) );
    Ariadne::ExactBoxType initialSetBox = { {-wid, wid}, {-hig, hig} };
    mpInitialSetBox.reset( new Ariadne::ExactBoxType( initialSetBox ) );

    // mWidthDist = std::uniform_real_distribution<>( initial[ 0 ].lower().get_d(), initial[ 0 ].upper().get_d() );
    // mHeightDist = std::uniform_real_distribution<>( initial[ 1 ].lower().get_d(), initial[ 1 ].upper().get_d() );

    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::RealConstant a( "a", Ariadne::Real( 1.4 ) ), b( "b", Ariadne::Real( 0.3 ) );
    Ariadne::EffectiveVectorFunction henon = Ariadne::make_function( {x, y}, {1 - a*x*x + y, b*x} );

    Ariadne::ConstraintSet cs = { cx*cx + cy*cy <= boundary };          // large value so it will take some time to diverge enough

    // figure out how to use bounds on safe region to do this
    Ariadne::ExactBoxType initialAbs = { {-boundary, boundary}, {-boundary, boundary} };

    mpRtree.reset( new ExactRefinementTree( initialAbs, cs, henon, Ariadne::Effort( 10 ) ) );
}

bool CegarTest::LoopTest::check() const
{
    bool foundCounterexample = false;

    std::uniform_real_distribution<> widthDist = std::uniform_real_distribution<>( (*mpInitialSetBox)[ 0 ].lower().get_d()
										   , (*mpInitialSetBox)[ 0 ].upper().get_d() )
	, heightDist = std::uniform_real_distribution<>( (*mpInitialSetBox)[ 1 ].lower().get_d()
							 , (*mpInitialSetBox)[ 1 ].upper().get_d() );

    std::vector< Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > > counterex;
    for( uint npt = 0; npt < mTestSize && !foundCounterexample; ++npt )
    {
	double x = widthDist( mRandom ), y = heightDist( mRandom );

    	Ariadne::ExactPoint initialPt( { x, y } );
	// refactor this into its own function
	Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > mappedPt = Ariadne::EffectiveVectorFunction::identity( initialPt.dimension() ).evaluate( initialPt );
	counterex = {initialPt};

    	for( uint nmap = 0; nmap < mTestSize && !foundCounterexample; ++nmap )
    	{
    	    mappedPt = mpRtree->dynamics().evaluate( mappedPt );
	    counterex.push_back( mappedPt );

	    Ariadne::ExactBoxType mappedBox( mappedPt.dimension() );
	    for( uint i = 0; i < mappedPt.dimension(); ++i )
		mappedBox[ i ] = {Ariadne::cast_exact( mappedPt[ i ].lower() ), Ariadne::cast_exact( mappedPt[ i ].upper() )};
	    if( definitely( mpRtree->constraints().separated( mappedBox ) ) )
		foundCounterexample = true;
    	}
    }

    auto cegarCounterex = cegar( *mpRtree, *mpInitialSet, Ariadne::Effort( 10 ), mRefinement, mLocator, MAX_NODES_FACTOR * mTestSize );
    std::cout << "safety " << cegarCounterex.first << std::endl;
      
    if( foundCounterexample && definitely( cegarCounterex.first ) )
    {
	std::cout << "found ";
	for_each( counterex.begin(), counterex.end()
		  , [] (const Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > >& p) {std::cout << p << " -> "; } );
	std::cout << std::endl << "but cegar found nothing" << std::endl;
	return false;
    }
    if( !foundCounterexample && !cegarCounterex.second.empty() )
    {                                                                 // because mapping random points comes with no guarantees
	std::cout << "cegar found ";
	printCounterexample( *mpRtree, cegarCounterex.second.begin(), cegarCounterex.second.end() ); // does endl
	std::cout << "but could not find a point that maps there" << std::endl;
    }
    // if( definitely( cegarCounterex.first ) ) std::cout << "cegar says: safe system" << std::endl;
    // else if( cegarCounterex.second.empty() ) std::cout << "cegar presumes: safe system" << std::endl;

    return true;
}




CegarTest::GROUP_CTOR( CegarTest, "top level cegar functions" )

void CegarTest::init()
{
    std::shared_ptr< ITestRunner > pRinterleave( new InterleaveRandomRunner() )
	, pStateless( new StatelessRunner() );
    addTest( new FindCounterexampleTest( 0.5 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new FindNoCounterexampleTest( 0.5 * mTestSize, 0.1 * mRepetitions ), pRinterleave );
    addTest( new LoopTest( 0.5 * mTestSize, 0.1 * mRepetitions ), pStateless );                        // alternatingly calling iterate then check is okay
}
