#ifndef CEGAR_TEST_HPP
#define CEGAR_TEST_HPP

#include "testGroupInterface.hpp"
#include "cegar.hpp"
#include "guide.hpp"

#include "expression/space.hpp"
#include "expression/expression.hpp"
#include "function/function.hpp"
// #include "function/function_set.hpp"

#include <random>

struct CegarTest : public ITestGroup
{
    typedef RefinementTree< Ariadne::ExactBoxType > ExactRefinementTree;

    static std::default_random_engine mRandom;

    static std::function< Ariadne::ValidatedUpperKleenean( const typename ExactRefinementTree::EnclosureT&
							   , const Ariadne::BoundedConstraintSet& ) > mIntersectConstraints;

    template< typename E >
    static E dummy( const E& e ) {return e; }

    template< typename E >
    static E extractEnclosure( const RefinementTree< E >& rtree, const typename RefinementTree< E >::MappingT::ValueT& val )
    {
    	if( val->isInside() )
    	    return tree::value( rtree.tree()
				, static_cast< const InsideGraphValue< typename RefinementTree< E >::RefinementT::NodeT >& >( *val ).treeNode() )->getEnclosure();
    	else
    	    return E( E::zero( 0 ) );
    }
    
    template< typename E, typename R >
    static typename RefinementTree< E >::NodeT refineRandomLeaf( RefinementTree< E >& rt, const R& refiner )
    {
	auto ls = rt.leaves();
	// need to store n otherwise graph part will be removed from memory (will be removed from graph)
	typename RefinementTree< E >::NodeT n = *(ls.begin() + (std::uniform_int_distribution<>( 0, ls.size() - 1 )( mRandom ) ) );
	rt.refine( n, refiner );
	return n;
    }

    template< typename Rtree >
    static void printNode( const Rtree& rtree, const typename Rtree::NodeT& n )
    {
	auto nVal = rtree.nodeValue( n );
	if( nVal )
	    std::cout << nVal.value().get().getEnclosure();
	else
	    std::cout << "[outside]";
	std::cout << " (" << rtree.isSafe( n ) << ") ";
    }
    
    template< typename Rtree, typename Iterator >
    static void printCounterexample( const Rtree& rtree, Iterator startCex, const Iterator& endCex )
    {
	for( ; startCex != endCex; ++startCex )
	{
	    printNode( rtree, *startCex );
	    std::cout << "-> ";
	}
	std::cout << std::endl;
    }

    typedef std::vector< Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > > BoundsCounterexampleT;
    
    static BoundsCounterexampleT findUnsafeTrajectoryFrom( const Ariadne::ExactPoint& initial
							   , const Ariadne::EffectiveVectorFunction& dynamics
							   , const Ariadne::BoundedConstraintSet& safeSet
							   , const uint& noMappings )
    {
	BoundsCounterexampleT path;
	Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > mapped = Ariadne::EffectiveVectorFunction::identity( initial.dimension() ).evaluate( initial );
	path.push_back( mapped );
	
	for( uint nmap = 0; nmap < noMappings; ++nmap )
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
							, const Ariadne::BoundedConstraintSet& safeSet
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

    static Ariadne::RealConstant constant( std::string name, double value )
    {
	return Ariadne::RealConstant( name, Ariadne::Real( value ) );
    }

    //! \return refinement tree for henon map with strictly positive initial and safe sets (i.e. left hand bottom corner is origin
    template< typename BoxType >
    static RefinementTree< BoxType >* henonMap( double safeWidth, double safeHeight
						, double a, double b, const Ariadne::Effort& e )
    {
	Ariadne::EffectiveScalarFunction cx = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 0 )
	    , cy = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( 2 ), 1 );

	Ariadne::BoundedConstraintSet safeSet( { {-safeWidth, safeWidth}, {-safeHeight, safeHeight} } );

	Ariadne::RealVariable x( "x" ), y( "y" );
	Ariadne::RealConstant fa( "a", Ariadne::Real( a ) ), fb( "b", Ariadne::Real( b ) );
	Ariadne::EffectiveVectorFunction henon = Ariadne::make_function( {x, y}, {1 - fa*x*x + y, fb*x} );

	return new RefinementTree< BoxType >( safeSet, henon, e );
    }

    // find a counterexample in an obviously unsafe system
    class FindCounterexampleTest : public ITest
    {
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::BoundedConstraintSet > mpInitial;
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
	std::unique_ptr< Ariadne::BoundedConstraintSet > mpInitial;
	LargestSideRefiner mRefiner;
	KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT > mGuide = KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT > ( 1 );
	STATEFUL_TEST( FindNoCounterexampleTest );
    };

    // no explicit test for isSpurious as it is hard to construct cases where a counterexample is definitely deemed spurious

    struct PrintInitialSet : public CegarObserver
    {
	template< typename Rtree, typename IterT >
	void searchCounterexample( const Rtree& rtree, IterT iBegin, const IterT& iEnd )
	{
	    std::cout << "current initial abstractions " << std::endl;
	    for( ; iBegin != iEnd; ++iBegin )
	    {
		printNode( rtree, *iBegin );
		std::cout << std::endl;
	    }
	    std::cout << std::endl;
	}
    };

    struct PrintCounterexample : public CegarObserver
    {
	template< typename Rtree, typename IterT >
	void processCounterexample( const Rtree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd )
	{
	    std::cout << "counterexample found " << std::endl;
	    for( ; iCounterexBegin != iCounterexEnd; ++iCounterexBegin )
	    {
		printNode( rtree, *iCounterexBegin );
		std::cout << " -> " << std::endl;
	    }
	    std::cout << std::endl;
	}
    };

    struct CheckInitialSet : public CegarObserver
    {
	typedef ExactRefinementTree Rtree;

	CheckInitialSet( const Rtree& rtree
			 , const Ariadne::ConstraintSet& initialSet
			 , const Ariadne::Effort& effort)
	    : mRtree( rtree )
	    , mInitialSet( initialSet )
	    , mEffort( effort )
	{}

	template< typename IterT >
	void searchCounterexample( const Rtree& rtree, IterT iBegin, const IterT& iEnd )
	{
	    if( !mpIssue )
	    {
		for( ; iBegin != iEnd; ++iBegin )
		{
		    auto nval = mRtree.nodeValue( *iBegin );
		    if( nval)
		    {
			const typename Rtree::EnclosureT& enc = nval.value().get().getEnclosure();
			if( definitely( mInitialSet.separated( enc ) ) ) 
			    mpIssue.reset( new typename Rtree::EnclosureT( enc ) );
		    }
		}
	    }
	}

	const Rtree& mRtree;
	Ariadne::ConstraintSet mInitialSet;
	Ariadne::Effort mEffort;
	std::unique_ptr< typename Rtree::EnclosureT > mpIssue;
    };

    template< typename Rtree >
    struct KeepInitialSet : public CegarObserver
    {
	template< typename IterT >
	void searchCounterexample( const Rtree& rtree, IterT iBegin, const IterT& iEnd )
	{
	    mNodes.clear();
	    mNodes.reserve( std::distance( iBegin, iEnd ) );
	    for( ; iBegin != iEnd; ++iBegin )
		mNodes.push_back( *iBegin );
	}

	std::vector< typename Rtree::NodeT > mNodes;
    };
    
    // test that initial abstractions are updated correctly
    class InitialAbstraction : public ITest
    {
	const double safeSetWidth = 10;
	const uint mMaxNodesFactor = 10;

    	std::unique_ptr< ExactRefinementTree > mpRtree;
    	std::unique_ptr< Ariadne::BoundedConstraintSet > mpInitialSet;
    	LargestSideRefiner mRefinement;
    	CompleteCounterexample mLocator;
    	KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT > mGuide = KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT >( 0.1 );
    	std::uniform_real_distribution<> mInitialBoxLengthDist = std::uniform_real_distribution<>( 0, 0.25 );
	std::exponential_distribution<> mDeltaDist = std::exponential_distribution<>( 1 );

    	STATELESS_TEST( InitialAbstraction );
    };

    // verify safety of trivially safe system
    class VerifySafety : public ITest
    {
	static const uint MAX_NODES_FACTOR = 100;
    	std::unique_ptr< ExactRefinementTree > mpRtree;
    	std::unique_ptr< Ariadne::BoundedConstraintSet > mpInitialSet;
    	LargestSideRefiner mRefinement;
    	CompleteCounterexample mLocator;
	KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT > mGuide = KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT >( 0.1 );
	std::exponential_distribution<> mInitialBoxLengthDist = std::exponential_distribution<>( 10 )
	    , mDeltaDist = std::exponential_distribution<>( 0.5 );
	std::uniform_real_distribution<> mAttractionDist = std::uniform_real_distribution( 0.1, 0.9 );
	
    	STATELESS_TEST( VerifySafety );
    };

    struct CounterexampleVerifier : public CegarObserver
    {
	template< typename IterT >
	void processCounterexample( const ExactRefinementTree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd )
	{
	    if( !mBadCounterexample.empty() )
		return;

	    auto icex = iCounterexBegin;
	    auto currentVal = rtree.nodeValue( *iCounterexBegin ), lastVal = currentVal;

	    while( currentVal && mBadCounterexample.empty() && icex != iCounterexEnd )
	    {
		++icex;
		lastVal = currentVal;

		if( icex != iCounterexEnd )
		{
		    currentVal = rtree.nodeValue( *icex );
		    if( currentVal )
		    {
			Ariadne::UpperBoxType mapped = Ariadne::image( lastVal.value().get().getEnclosure(), rtree.dynamics() );
			if( definitely( Ariadne::intersection( mapped, currentVal.value().get().getEnclosure() ).is_empty() ) )
			{
			    mBadCounterexample = std::vector< typename ExactRefinementTree::NodeT >( iCounterexBegin, iCounterexEnd );
			    mBadLink = mBadCounterexample.begin();
			    std::advance( mBadLink, std::distance( iCounterexBegin, icex ) );
			}
		    }
		}
	    }
	}

	std::vector< typename ExactRefinementTree::NodeT > mBadCounterexample;
	std::vector< typename ExactRefinementTree::NodeT >::const_iterator mBadLink;
    };
    
    //! \class tests that successive states in counterexamples are reachable
    class VerifyCounterexamples : public ITest
    {
	static const uint MAX_NODES_FACTOR = 25;
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::BoundedConstraintSet > mpInitialSet;
	LargestSideRefiner mRefinement;
	CompleteCounterexample mLocator;
	KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT > mGuide = KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT >( 0.1 );
	std::exponential_distribution<> mInitialBoxLengthDist = std::exponential_distribution<>( 25 )    
	    , mSafeBoxLengthDist = std::exponential_distribution<>( 0.01 );                              // narrow initial set, wide safe set => many counterexamples

	STATELESS_TEST( VerifyCounterexamples );
    };

    // test that counterexample with single broken link is detected
    class LoopTest : public ITest
    {
	static const uint MAX_NODES_FACTOR = 100; // because: 50 * 200 = b10,000
	std::exponential_distribution<> mInitialBoxLengthDist = std::exponential_distribution<>( 3 ) // so mean is 0.33, quite close to safe region
	    , mSafeSetBoxLengthDist = std::exponential_distribution<>( 0.1 );
	std::unique_ptr< ExactRefinementTree > mpRtree;
	std::unique_ptr< Ariadne::BoundedConstraintSet > mpInitialSet;
	LargestSideRefiner mRefinement;
	CompleteCounterexample mLocator;
	KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT > mGuide = KeepRandomCounterexamples< typename ExactRefinementTree::EnclosureT >( 0.1 );

	STATELESS_TEST( LoopTest );
    };

    GROUP_CTOR_DECL( CegarTest );
    
    void init();
};


#endif //CEGAR_TEST_HPP
