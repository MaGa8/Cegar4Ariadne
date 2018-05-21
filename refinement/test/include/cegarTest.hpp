#ifndef CEGAR_TEST_HPP
#define CEGAR_TEST_HPP

#include "testGroupInterface.hpp"
#include "cegar.hpp"

#include <random>

struct CegarTest : public ITestGroup
{
    typedef RefinementTree< Ariadne::ExactBoxType > ExactRefinementTree;

    static std::default_random_engine mRandom;

    template< typename E >
    static typename RefinementTree< E >::NodeT refineRandomLeaf( RefinementTree< E >& rt, const IRefinement< E >& refiner )
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

    // find a counterexample in an obviously unsafe system
    class FindCounterexampleTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	Ariadne::ExactBoxType mInitial;
	LargestSideRefiner< Ariadne::ExactIntervalType > mRefiner;
	STATEFUL_TEST( FindCounterexampleTest );
    };

    // do not find a counterexample in a system that is safe
    // counterexamples might still be reported, by they are either spurious or terminate in a state whose safety cannot be determined
    // implicitly tests isSpurious
    class FindNoCounterexampleTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	Ariadne::ExactBoxType mInitial;
	LargestSideRefiner< Ariadne::ExactIntervalType > mRefiner;
	STATEFUL_TEST( FindNoCounterexampleTest );
    };

    // no explicit test for isSpurious as it is hard to construct cases where a counterexample is definitely deemed spurious

    // test that counterexample with single broken link is detected
    class LoopTest : public ITest
    {
	static const uint MAX_NODES_FACTOR = 20; // because: 50 * 200 = b10,000
	std::unique_ptr< ExactRefinementTree > mpRtree;
	Ariadne::ExactBoxType mInitial;
	LargestSideRefiner< Ariadne::ExactIntervalType > mRefinement;
	CompleteCounterexample< Ariadne::ExactBoxType > mLocator;

	STATELESS_TEST( LoopTest );
    };

    GROUP_CTOR_DECL( CegarTest );
    
    void init();
};


#endif CEGAR_TEST_HPP
