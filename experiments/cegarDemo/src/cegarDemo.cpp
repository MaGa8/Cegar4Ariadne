#include "cegar.hpp"
#include "guide.hpp"
#include "locator.hpp"
#include "cegarObserver.hpp"
#include "visualization.hpp"

#include "numeric/real.hpp"
#include "expression/expression.hpp"
#include "expression/space.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "geometry/interval.hpp"
#include "output/graphics.hpp"


struct TreeSizePrinter : public CegarObserver
{
    template< typename Rtree >
    void startIteration( const Rtree& rtree )
    {
	std::cout << "new iteration, number of nodes " << rtree.tree().size() << std::endl;
    }
};

int main()
{
    typedef RefinementTree< Ariadne::ExactBoxType > ExactRefinementTree;
    
    Ariadne::Effort effort( 10 );
    uint maxNodes = 100000;

    Ariadne::BoundedConstraintSet hugeSafeSet( { {-20, 20}, {-20, 20} } )
	, wideSafeSet( { {-10,10}, {-10,10} } )
	, medSafeSet( { {-5,5}, {-5,5} } )
	, medSmallSafeSet( { {-3.5, 3.5}, {-3.5, 3.5} } )
	, tinkerbellSafeSet( { {-3, 3}, {-3,3} } )
	, smallSafeSet( { {-2.5,2.5}, {-2.5,2.5} } )
	, tinySafeSet( { {-1.25,1.25}, {-1.25,1.25} } )
	, initialSet( { {0, 0.4}, {0, 0.4} } );

    // need refinement tree with non-static dynamics
    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::Space< Ariadne::Real > fspace = {x, y};
    	
    Ariadne::RealConstant a( "a", Ariadne::Real( 1.4 ) ), b( "b", Ariadne::Real( 0.3 ) );
    // Henon map - needs less nodes than "henonSimple"
    // can do with 5k: (1.3, 0.3) or (1.4, 0.25)
    // with 7.5k (1.35, 0.3)
    Ariadne::EffectiveVectorFunction henon = Ariadne::make_function( fspace, {1 - a*x*x + y, b*x} );

    // can do with 5k: (1.3, 0.3) or (1.4, 0.25)
    Ariadne::EffectiveVectorFunction henonSimple = Ariadne::make_function( fspace, {1 - a*x*x + b*y, x} );
    // logistic map
    // within 0.5k (0.75, 0.65)
    // within 2k (0.9, 0.9)
    // within 5k (0.99, 0.95)
    // within 10k (0.99, 0.98)
    // within 15k (1.1, 0.98)
    Ariadne::RealConstant rx( "rx", Ariadne::Real( 0.99 ) ), ry( "ry", Ariadne::Real( 0.99 ) );
    Ariadne::EffectiveVectorFunction logistic = Ariadne::make_function( fspace, {rx*x*(1 - x), ry*y*(1 - y)} );
    // bogdanov map
    Ariadne::RealConstant eps( "e", Ariadne::Real( 0.005 ) ), k( "k", Ariadne::Real( 0.01 ) ), mu( "mu", Ariadne::Real( 0.0001 ) );
    Ariadne::EffectiveVectorFunction bogdanov = Ariadne::make_function( fspace, {x + y + eps*y + k*x*(x - 1) + mu*x*y, y + eps*y + k*x*(x - 1) + mu*x*y } );

    Ariadne::RealConstant at( "a", Ariadne::Real( 0.9 ) )
	, bt( "b", Ariadne::Real( -0.6013 ) )
	, ct( "c", Ariadne::Real( 2 ) )
	, dt( "d", Ariadne::Real( 0.5 ) );
    Ariadne::EffectiveVectorFunction tinkerbell = Ariadne::make_function( {x,y}, {x*x - y*y + at*x + bt*y, 2*x*y + ct*x + dt*y} );

    ExactRefinementTree rtree( tinySafeSet, logistic, effort );

    LargestSideRefiner refiner;
    // CompleteCounterexample locator;
    // RandomStates locator;
    StateVolume stateH;
    GreatestState counterexH;
    TreeSizePrinter sizePrinter;
    DebugOutput dbgout;

    LimitedIterations termination( maxNodes );
    
    auto safety = cegar( rtree, initialSet, effort, refiner, stateH, counterexH, termination, dbgout );
    std::cout << "safety " << safety.first << std::endl;
    std::cout << "counterexample " << std::endl;
    for( auto& n : safety.second )
    {
	auto nval = rtree.nodeValue( n );
	if( nval )
	    std::cout << nval.value().get().getEnclosure();
	else
	    std::cout << "[outside]";
	std::cout << std::endl;
	
	// use after merging no-tree branch
	// auto pgval = graph::value( rtree.leafMapping(), n );
	// if( pgval->isInside() )
	//     std::cout << static_cast< InsideGraphValue< Ariadne::ExactBoxType >& >( *pgval ) << std::endl;
	// else
	//     std::cout << static_cast< OutsideGraphValue& >( *pgval ) << std::endl;
    }
    
    // trying graphical output
    auto fig = visualize( rtree, initialSet, effort );
    fig.write( "test.png" );
    
    return 0;
}
