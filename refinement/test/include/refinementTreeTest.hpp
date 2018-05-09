#ifndef REFINEMENT_TREE_TEST
#define REFINEMENT_TREE_TEST

#include "testGroupInterface.hpp"
#include "refinementTree.hpp"
#include "refinementStrategy.hpp"

#include "expression/space.hpp"
#include "expression/expression.hpp"
#include "function/function.hpp"
#include "function/function.decl.hpp"

#include <memory>
#include <random>

struct RefinementTreeTest : public ITestGroup
{
    typedef RefinementTree< Ariadne::EffectiveIntervalType > EffectiveRefinementTree;
    typedef RefinementTree< Ariadne::ExactIntervalType > ExactRefinementTree;

    // typedef RefinementTree< int > EffectiveRefinementTree;
    // typedef RefinementTree< int > ExactRefinementTree;
    
    static std::default_random_engine mRandom;

    template< typename IntervalT >
    static typename RefinementTree< IntervalT >::NodeT refineRandomLeaf( RefinementTree< IntervalT >& rt, const Ariadne::Box< IntervalT >& rootBox, const IRefinementStrategy< IntervalT >& refiner )
    {
	auto ls = rt.image( rootBox ); // should select all leaves, right?
	// need to store n otherwise graph part will be removed from memory (will be removed from graph)
	typename RefinementTree< IntervalT >::NodeT n = *(ls.begin() + (std::uniform_int_distribution<>( 0, ls.size() - 1 )( mRandom ) ) );
	rt.refine( n, refiner );
	return n;
    }

    template< typename IntervalT >
    static void refineEqualDepth( RefinementTree< IntervalT >& rt, const uint depth )
    {
	LargestSideRefiner< IntervalT > refiner;
	for( uint i = 0; i < depth; ++i )
	{
	    auto lvs = rt.leaves();
	    for( typename RefinementTree< IntervalT >::NodeT& lf : lvs )
		rt.refine( lf, refiner );
	}
    }

    template< typename IntervalT >
    static bool nodeEquals( const RefinementTree< IntervalT >& rt, const typename RefinementTree< IntervalT >::NodeT& n1, const typename RefinementTree< IntervalT >::NodeT& n2 )
    {
	return tree::value( rt.getTree()
			    , graph::value( rt.getLeafMapping(), n1 ) ).getEnclosure() == tree::value( rt.getTree()
													, graph::value( rt.getLeafMapping(), n2 ) ).getEnclosure();
    }

    template< typename IntervalT >
    static RefinementTree< IntervalT > getDefaultTree( const Ariadne::Box< IntervalT > rootBox, uint cap )
    {
	Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
	Ariadne::EffectiveScalarFunction constraintExpression = (a + b);
	// only upper bound on + values
	Ariadne::EffectiveConstraint c1 = (constraintExpression <= cap );
	Ariadne::EffectiveConstraintSet cs = { c1 };
	// static dynamics
	Ariadne::RealVariable x( "x" ), y( "y" );
	Ariadne::Space< Ariadne::Real > vspace = {x, y};
	Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x, y} );
	return RefinementTree< IntervalT >( rootBox, cs, f, Ariadne::Effort( 5 ) );
    }
    
    // test expansion:
    // proper number of nodes in graph: number of refinements - 1
    // proper depth: increasing if refining deepest level node otherwise stagnant
    class ExpansionTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ExactBoxType > mpBox;
	std::unique_ptr< IRefinementStrategy< Ariadne::ExactIntervalType > > mpRefiner;
	const uint EXPANSION_SIZE;
	uint mPreviousNoNodes, mPreviousHeight, mExpandNodeDepth;
	STATEFUL_TEST( ExpansionTest );
    };
    
    // test leaves after couple of expansions
    class LeavesTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ExactBoxType > mpBox;
	const uint EXPANSION_SIZE = 2;
	uint mExpansionCounter;
	STATEFUL_TEST( LeavesTest );
    };

    // test image: brute force all leaves and test for intersection
    class ImageTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ExactBoxType > mpRootBox;
	STATEFUL_TEST( ImageTest );
    };
    
    // test whether non-leafs are absent from tree
    class NonLeafRemovalTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ExactBoxType > mpRootBox;
	LargestSideRefiner< Ariadne::ExactIntervalType > mRefiner;
	STATEFUL_TEST( NonLeafRemovalTest );
	
    };

    // test preimage after some expansions: 
    class PreimageTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ExactBoxType > mpRootBox;
	typename ExactRefinementTree::NodeT mRefined;
	LargestSideRefiner< Ariadne::ExactIntervalType > mRefiner;
	STATEFUL_TEST( PreimageTest );
    };

    // test postimage after some expansions
    class PostimageTest : public ITest
    {
    	std::unique_ptr< ExactRefinementTree > mpRtree;
    	std::unique_ptr< Ariadne::ExactBoxType > mpRootBox;
    	typename ExactRefinementTree::NodeT mRefined;
    	LargestSideRefiner< Ariadne::ExactIntervalType > mRefiner;
    	STATEFUL_TEST( PostimageTest );
    };

    //positve test for finding counterexamples
    class PositiveCounterexampleTest : public ITest
    {
    	STATELESS_TEST( PositiveCounterexampleTest );
    };
    
    // negative test for finding counterexamples

    GROUP_CTOR_DECL( RefinementTreeTest );

    void init();
};

#endif
