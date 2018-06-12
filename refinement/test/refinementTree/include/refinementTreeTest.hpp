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
		auto& inVal = static_cast< InsideGraphValue< E >& >( *val );
		return inVal.getEnclosure();
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
	auto vrange = graph::vertices( rt.graph() );
	// need to store n otherwise graph part will be removed from memory (will be removed from graph)
	typename RefinementTree< E >::NodeT n;
	do
	{
	    uint jump = std::uniform_int_distribution<>( 0, std::distance( vrange.first, vrange.second ) - 1 )( mRandom );
	    auto irefine = vrange.first;
	    std::advance( irefine, jump );
	    n = *irefine;
	} while( rt.equal( rt.outside(), n ) );
							
	rt.refine( n, refiner );
	return n;
    }

    template< typename E >
    static void refineEqualDepth( RefinementTree< E >& rt, const uint depth )
    {
	LargestSideRefiner refiner;
	for( uint i = 0; i < depth; ++i )
	{
	    auto vrange = graph::vertices( rt.graph() );
	    for( auto inode = vrange.first; inode != vrange.second; ++inode )
		rt.refine( *inode, refiner );
	}
    }

    template< typename IntervalT >
    static RefinementTree< Ariadne::Box< IntervalT > > getDefaultTree( const Ariadne::Box< IntervalT > rootBox, uint cap )
    {
	Ariadne::EffectiveScalarFunction a = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	, b = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );
	Ariadne::EffectiveScalarFunction constraintExpression = (a + b);
	// static dynamics
	Ariadne::RealVariable x( "x" ), y( "y" );
	Ariadne::Space< Ariadne::Real > vspace = {x, y};
	Ariadne::EffectiveVectorFunction f = Ariadne::make_function( vspace, {x, y} );
	return RefinementTree< Ariadne::Box< IntervalT > >( Ariadne::BoundedConstraintSet( Ariadne::RealBox( rootBox ) )
							    , f, Ariadne::Effort( 5 ) );
    }

    template< typename E >
    static void printNodeValue( const std::optional< std::reference_wrapper< const InsideGraphValue< E > > >& otn )
    {
	if( otn )
	    std::cout << otn.value().get().getEnclosure() << " ";
	else
	    std::cout << "[unsafe] ";
    }
    
    // test leaves after couple of expansions
    class SizeTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	const uint EXPANSION_SIZE = 2;
	LargestSideRefiner mRefiner;
	uint mExpansionCounter;
	STATEFUL_TEST( SizeTest );
    };

    // test intersection: brute force all leaves and test for intersection
    class IntersectionTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	LargestSideRefiner mRefiner;
	STATEFUL_TEST( IntersectionTest );
    };

    // test intersection with random circular constraint set
    class CSetIntersectionTest : public ITest
    {
    	std::unique_ptr< ExactRefinementTree > mpRtree;
    	LargestSideRefiner mRefiner;
    	STATEFUL_TEST( CSetIntersectionTest );
    };
    
    // test whether non-leafs are absent from graph
    class RefinedNodesRemovalTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	LargestSideRefiner mRefiner;
	STATEFUL_TEST( RefinedNodesRemovalTest );
	
    };

    // test preimage after some expansions: 
    class PreimageTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	typename ExactRefinementTree::NodeT mRefined;
	LargestSideRefiner mRefiner;
	STATEFUL_TEST( PreimageTest );
    };

    // test postimage after some expansions
    class PostimageTest : public ITest
    {
    	std::unique_ptr< ExactRefinementTree > mpRtree;
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

