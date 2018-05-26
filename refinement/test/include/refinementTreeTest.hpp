#ifndef REFINEMENT_TREE_TEST
#define REFINEMENT_TREE_TEST

#include "testGroupInterface.hpp"
#include "cegar.hpp"

#include "expression/space.hpp"
#include "expression/expression.hpp"
#include "function/function.hpp"
#include "function/function.decl.hpp"

#include <memory>
#include <random>

struct RefinementTreeTest : public ITestGroup
{
    // typedef RefinementTree< Ariadne::EffectiveIntervalType > EffectiveRefinementTree;
    typedef RefinementTree< Ariadne::Box< Ariadne::ExactIntervalType > > ExactRefinementTree;

    // typedef RefinementTree< int > EffectiveRefinementTree;
    // typedef RefinementTree< int > ExactRefinementTree;
    
    static std::default_random_engine mRandom;
    
    template< typename E >
    class GraphVertexPrintConverter
    {
      public:
	GraphVertexPrintConverter( const RefinementTree< E >& rtree ) : mRtree( rtree ) {}
	
	E operator ()( const typename RefinementTree< E >::MappingT::ValueT& val )
	{
	    if( val->isInside() )
	    {
		auto& inVal = static_cast< InsideGraphValue< typename RefinementTree< E >::RefinementT::NodeT >& >( *val );
		return tree::value( mRtree.tree(), inVal.treeNode() )->getEnclosure();
	    }
	    else
		return E::zero( 2 ); // hard coded cause using 2d tests only
	}
	
      private:
	const RefinementTree< E >& mRtree;
    };

    template< typename E, typename R >
    static typename RefinementTree< E >::NodeT refineRandomLeaf( RefinementTree< E >& rt, const R& refiner )
    {
	auto ls = rt.leaves();
	// need to store n otherwise graph part will be removed from memory (will be removed from graph)
	typename RefinementTree< E >::NodeT n = *(ls.begin() + (std::uniform_int_distribution<>( 0, ls.size() - 1 )( mRandom ) ) );
	rt.refine( n, refiner );
	return n;
    }

    template< typename E >
    static void refineEqualDepth( RefinementTree< E >& rt, const uint depth )
    {
	LargestSideRefiner refiner;
	for( uint i = 0; i < depth; ++i )
	{
	    auto lvs = rt.leaves();
	    for( typename RefinementTree< E >::NodeT& lf : lvs )
		rt.refine( lf, refiner );
	}
    }

    template< typename E >
    static bool nodeEquals( const RefinementTree< E >& rt, const typename RefinementTree< E >::NodeT& n1, const typename RefinementTree< E >::NodeT& n2 )
    {
	return rt.nodeValue( n1 ) == rt.nodeValue( n2 );
    }

    template< typename IntervalT >
    static RefinementTree< Ariadne::Box< IntervalT > > getDefaultTree( const Ariadne::Box< IntervalT > rootBox, uint cap )
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
	return RefinementTree< Ariadne::Box< IntervalT > >( rootBox, cs, f, Ariadne::Effort( 5 ) );
    }

    template< typename RefTree >
    static void printNodeValue( const std::optional< std::reference_wrapper< const InteriorTreeValue< typename RefTree::EnclosureT > > > otn )
    {
	if( otn )
	    std::cout << otn.value().get().getEnclosure() << " ";
	else
	    std::cout << "[unsafe] ";
    }
    
    // test expansion:
    // proper number of nodes in graph: number of refinements - 1
    // proper depth: increasing if refining deepest level node otherwise stagnant
    class ExpansionTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ExactBoxType > mpBox;
	LargestSideRefiner mRefiner;
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

    // test intersection: brute force all leaves and test for intersection
    class IntersectionTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ExactBoxType > mpRootBox;
	LargestSideRefiner mRefiner;
	STATEFUL_TEST( IntersectionTest );
    };

    // test intersection with random circular constraint set
    class CSetIntersectionTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ExactBoxType > mpRootBox;
	LargestSideRefiner mRefiner;
	STATEFUL_TEST( CSetIntersectionTest );
    };
    
    // test whether non-leafs are absent from graph
    class NonLeafRemovalTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ExactBoxType > mpRootBox;
	LargestSideRefiner mRefiner;
	STATEFUL_TEST( NonLeafRemovalTest );
	
    };

    // test preimage after some expansions: 
    class PreimageTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::ExactBoxType > mpRootBox;
	typename ExactRefinementTree::NodeT mRefined;
	LargestSideRefiner mRefiner;
	STATEFUL_TEST( PreimageTest );
    };

    // test postimage after some expansions
    class PostimageTest : public ITest
    {
    	std::unique_ptr< ExactRefinementTree > mpRtree;
    	std::unique_ptr< Ariadne::ExactBoxType > mpRootBox;
    	typename ExactRefinementTree::NodeT mRefined;
    	LargestSideRefiner mRefiner;
    	STATEFUL_TEST( PostimageTest );
    };
    
    // test that there exists a node in static unsafe system that maps to always unsafe
    class AlwaysUnsafeTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	LargestSideRefiner mRefiner;
	STATEFUL_TEST( AlwaysUnsafeTest );
    };

    // move away to cegar test
    // //positve test for finding counterexamples
    // class PositiveCounterexampleTest : public ITest
    // {
    // 	STATELESS_TEST( PositiveCounterexampleTest );
    // };
    
    // negative test for finding counterexamples

    GROUP_CTOR_DECL( RefinementTreeTest );

    void init();
};

#endif

