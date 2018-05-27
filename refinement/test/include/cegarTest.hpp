#ifndef CEGAR_TEST_HPP
#define CEGAR_TEST_HPP

#include "testGroupInterface.hpp"
#include "cegar.hpp"
#include "guide.hpp"

#include <random>

struct CegarTest : public ITestGroup
{
    typedef RefinementTree< Ariadne::ExactBoxType > ExactRefinementTree;

    static std::default_random_engine mRandom;

    static std::function< Ariadne::ValidatedLowerKleenean( const typename ExactRefinementTree::EnclosureT&
							   , const Ariadne::ConstraintSet& ) > mIntersectConstraints;

    template< typename E, typename R >
    static typename RefinementTree< E >::NodeT refineRandomLeaf( RefinementTree< E >& rt, const R& refiner )
    {
	auto ls = rt.leaves();
	// need to store n otherwise graph part will be removed from memory (will be removed from graph)
	typename RefinementTree< E >::NodeT n = *(ls.begin() + (std::uniform_int_distribution<>( 0, ls.size() - 1 )( mRandom ) ) );
	rt.refine( n, refiner );
	return n;
    }

    template< typename Rtree, typename Iterator >
    static void printCounterexample( const Rtree& rtree, Iterator startCex, const Iterator& endCex )
    {
	for( ; startCex != endCex; ++startCex )
	{
	    auto nVal = rtree.nodeValue( *startCex );
	    if( nVal )
		std::cout << nVal.value().get().getEnclosure();
	    else
		std::cout << "[outside]";
	    std::cout << " (" << rtree.isSafe( *startCex ) << ") -> ";
	}
	std::cout << std::endl;
    }

    typedef std::vector< Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > > BoundsCounterexampleT;
    
    static BoundsCounterexampleT findUnsafeTrajectoryFrom( const Ariadne::ExactPoint& initial
							   , const Ariadne::EffectiveVectorFunction& dynamics
							   , const Ariadne::ConstraintSet& safeSet
							   , const uint& noMappings )
    {
	BoundsCounterexampleT path;
	Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > mapped = Ariadne::EffectiveVectorFunction::identity( initial.dimension() ).evaluate( initial );
	path.push_back( mapped );
	
	for( int nmap = 0; nmap < noMappings; ++nmap )
	{
	    mapped = dynamics.evaluate( mapped );
	    path.push_back( mapped );
	    Ariadne::ExactBoxType ptBox = boundsPoint2Box( mapped );
	    if( definitely( safeSet.separated( ptBox ) ) )
		return path;
	}
	return {};
    }

    //! \return true if every point sampled did not leave the safe set when mapped noMappings times forward
    static BoundsCounterexampleT findUnsafeTrajectory ( std::uniform_real_distribution<>& widthDist
							, std::uniform_real_distribution<>& heightDist
							, const Ariadne::EffectiveVectorFunction& dynamics
							, const Ariadne::ConstraintSet& safeSet
							, const uint& noSamples
							, const uint& noMappings )
    {
	bool foundUnsafe = false;
	for( uint npt = 0; npt < noSamples && !foundUnsafe; ++npt )
	{
	    double x = widthDist( mRandom ), y = heightDist( mRandom );
	    Ariadne::ExactPoint initialPt( { x, y } );

	    auto unsafeTrajectory = findUnsafeTrajectoryFrom( initialPt, dynamics, safeSet, noMappings );
	    if( !unsafeTrajectory.empty() )
		return unsafeTrajectory;
	}
	return {};
    }

    // find a counterexample in an obviously unsafe system
    class FindCounterexampleTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ConstraintSet > mpInitial;
	LargestSideRefiner mRefiner;
	KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT > mGuide = KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT >( 1 );
	STATEFUL_TEST( FindCounterexampleTest );
    };

    // do not find a counterexample in a system that is safe
    // counterexamples might still be reported, by they are either spurious or terminate in a state whose safety cannot be determined
    // implicitly tests isSpurious
    class FindNoCounterexampleTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ConstraintSet > mpInitial;
	LargestSideRefiner mRefiner;
	KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT > mGuide = KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT > ( 1 );
	STATEFUL_TEST( FindNoCounterexampleTest );
    };

    // no explicit test for isSpurious as it is hard to construct cases where a counterexample is definitely deemed spurious

    // test that counterexample with single broken link is detected
    class LoopTest : public ITest
    {
	static const uint MAX_NODES_FACTOR = 20; // because: 50 * 200 = b10,000
	std::exponential_distribution<> mInitialBoxLengthDist = std::exponential_distribution<>( 100 )
	    , mDeltaDist = std::exponential_distribution<>( 0.5 );
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ConstraintSet > mpInitialSet;
	std::unique_ptr< Ariadne::ExactBoxType > mpInitialSetBox;
	LargestSideRefiner mRefinement;
	CompleteCounterexample mLocator;
	KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT > mGuide = KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT >( 0.1 );

	STATELESS_TEST( LoopTest );
    };

    // class VerifySafety : public ITest
    // {
    // 	std::unique_ptr< ExactRefinementTree > mpInitialSet;
    // 	std::unique_ptr< Ariadne::ConstraintSet > mpInitialSet;
    // 	LargestSideRefiner mRefinement;
    // 	CompleteCounterexample mLocator;
    // 	KeepRandomCounterexamples mGuide( 0.1 );

    // 	STATELESS_TEST( VerifySafety );
    // };

    GROUP_CTOR_DECL( CegarTest );
    
    void init();
};


#endif CEGAR_TEST_HPP
