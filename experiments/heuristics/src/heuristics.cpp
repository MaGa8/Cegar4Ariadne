#include "heuristics.hpp"

#include "guide.hpp" // others should also be removed from refinementTree, cegar and heuristics.hpp for cleanliness

#include "expression/space.hpp"
#include "expression/expression.hpp"
#include "function/function.hpp"

#include <array>
#include <sstream>

std::ostream& operator <<( std::ostream& os, const Metrics& mets )
{
    auto printAll = [&os] (auto& met) {
			os << "," << met.first << ",";
			std::for_each( met.begin(), met.end(), [&os] (auto& x) {os << x << ",";} );
			os << std::endl;
		    };
    auto metTuple = std::make_tuple( mets.mTotal, mets.mSearch, mets.mCheck, mets.mRefine, mets.mMisc, mets.mIterations, mets.mAverage, mets.mResult );
    applyTuple( printAll, metTuple );
    // os << mets.first << ",";
    // printAll( mets.mTotal ); os << std::endl;
    // os << std::endl << ",";
    // printAll( mets.mSearch ); os << std::endl;
    // os << std::endl << ",";
    // printAll( mets.mCheck ); os << std::endl;
    // os << std::endl << ",";
    // printAll( mets.mRefine ); os << std::endl;
    // os << std::endl << ",";
    // printAll( mets.mMisc ); os << std::endl;
    // os << std::endl << ",";
    // printAll( mets.mIterations ); os << std::endl;
    // os << std::endl << ",";
    // printAll( mets.mAverage ); os << std::endl;
    return os;
}

Ariadne::RealConstant make_constant( const std::string& name, const double& val )
{
    return Ariadne::RealConstant( name, Ariadne::Real( val ) );
}

System< Ariadne::ExactBoxType > logisticMap( const double& rx, const double& ry
					     , const double& wsafe, const double& winitial
					     , const double& delta, const Ariadne::Effort& effort )
{
    Ariadne::RealConstant rxc = make_constant( "rx", rx ), ryc = make_constant( "ry", ry );

    Ariadne::BoundedConstraintSet initial( { {0, winitial}, {0, winitial} } )
	, safe( { {-delta, wsafe}, {-delta, wsafe} } );

    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::EffectiveVectorFunction logistic = Ariadne::make_function( {x, y}, { rxc*x*(1 - x), ryc*y*(1 - y)} );

    std::stringstream ss;
    ss << "logistic_map_rx" << rx << "_ry" << ry << "_s" << wsafe << "_i" << winitial;

    return System< Ariadne::ExactBoxType >( initial, safe, logistic, effort, ss.str() );
}

System< Ariadne::ExactBoxType > henonMap( const double& a, const double& b
					  , const double& xsafeLower, const double& xsafeUpper
					  , const double& ysafeLower, const double& ysafeUpper
					  , const double& xinitialLower, const double& xinitialUpper
					  , const double& yinitialLower, const double& yinitialUpper
					  , const Ariadne::Effort& effort )
{
    Ariadne::RealConstant ac = make_constant( "a", a ), bc = make_constant( "b", b );

    Ariadne::BoundedConstraintSet initial( Ariadne::RealBox( { {xinitialLower, xinitialUpper}, {yinitialLower, yinitialUpper} } ) )
	, safe( Ariadne::RealBox( { {xsafeLower, xsafeUpper}, {ysafeLower, ysafeUpper} } ) );

    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::EffectiveVectorFunction logistic = Ariadne::make_function( {x, y}, {1 - ac*x*x + y, bc*x} );

    std::stringstream ss;
    ss << "henon_map_a" << a << "_b" << b;
    return System< Ariadne::ExactBoxType >( initial, safe, logistic, effort, ss.str() );
}

System< Ariadne::ExactBoxType > bogdanovMap( const double& epsilon, const double& k, const double& mu
					      , const double& wsafe, const double winitial, const Ariadne::Effort& effort )
{
    Ariadne::BoundedConstraintSet initial( { {-winitial, winitial}, {-winitial, winitial} } )
	, safe( { {-wsafe, wsafe}, {-wsafe, wsafe} } );

    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::RealConstant epsc = make_constant( "e", epsilon )
	, kc = make_constant( "k", k )
	, muc = make_constant( "mu", mu );
    Ariadne::EffectiveVectorFunction logistic = Ariadne::make_function( {x, y}, {x + y, y + epsc*y + kc*x*(x - 1) + muc*x*y} );
    
    std::stringstream ss;
    ss << "bogdanov_map_eps" << epsilon << "_k" << k << "_mu" << mu << "_s" << wsafe << "_i" << winitial;
    return System< Ariadne::ExactBoxType >( initial, safe, logistic, effort, ss.str() );
}

System< Ariadne::ExactBoxType > tinkerbell( const double& a, const double& b, const double& c, const double& d
					    , const double& xsafeLower, const double& xsafeUpper
					    , const double& ysafeLower, const double& ysafeUpper
					    , const double& xinitialLower, const double& xinitialUpper
					    , const double& yinitialLower, const double& yinitialUpper
					    , const Ariadne::Effort effort )
{
    Ariadne::BoundedConstraintSet initial( Ariadne::RealBox( { {xinitialLower, xinitialUpper}, {yinitialLower, yinitialUpper} } ) )
	, safe( Ariadne::RealBox( { {xsafeLower, xsafeUpper}, {ysafeLower, ysafeUpper} } ) ); 

    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::RealConstant ac = make_constant( "a", a )
	, bc = make_constant( "b", b )
	, cc = make_constant( "c", c )
	, dc = make_constant( "d", d );
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( {x,y}, {x*x - y*y + ac*x + bc*y, 2*x*y + cc*x + dc*y} );

    std::stringstream ss;
    ss << "tinkerbell_a" << a << "_b" << b << "_c" << c << "_d" << d;

    return System< Ariadne::ExactBoxType >( initial, safe, f, effort, ss.str() );
}

System< Ariadne::ExactBoxType > duffing( const double& a, const double& b
					  , const double& xsafeLower, const double& xsafeUpper
					  , const double& ysafeLower, const double& ysafeUpper
					  , const double& xinitialLower, const double& xinitialUpper
					  , const double& yinitialLower, const double& yinitialUpper
					  , const Ariadne::Effort effort )
{
    Ariadne::BoundedConstraintSet initial( Ariadne::RealBox( { {xinitialLower, xinitialUpper}, {yinitialLower, yinitialUpper} } ) )
	, safe( Ariadne::RealBox( { {xsafeLower, xsafeUpper}, {ysafeLower, ysafeUpper} } ) ); 

    Ariadne::RealVariable x( "x" ), y( "y" );
    Ariadne::RealConstant ac = make_constant( "a", a )
	, bc = make_constant( "b", b );
    Ariadne::EffectiveVectorFunction f = Ariadne::make_function( {x,y}, {y, -bc*x + ac*y - y*y*y} );

    std::stringstream ss;
    ss << "duffing_a" << a << "_b" << b;
    return System< Ariadne::ExactBoxType >( initial, safe, f, effort, ss.str() );
}



int main()
{
    typedef Ariadne::ExactBoxType EncT;
    Ariadne::Effort effort( 25 );
    std::vector< System< EncT > > systems;

    // std::vector< System< EncT > > systems = { logistic( 0.5, 1, 1, 0.4, 0.1, 10 ) };
    // std::vector< ExperimentConfiguration > configs = { ExperimentConfiguration
    System< EncT > logisticFalse = logisticMap( 0.98, 0.98, 0.35, 0.4, 0.1, effort );
    System< EncT > logisticEasy = logisticMap( 0.98, 0.98, 1.1, 0.4, 0.1, effort );
    System< EncT > logisticHard = logisticMap( 1.00, 0.99, 3, 0.4, 0.1, effort );
    System< EncT > logisticCrazy = logisticMap( 1.1, 0.99, 3, 0.4, 0.1, effort );

    System< EncT > henonFalse = henonMap( 1.4, 0.3
					  , -2, 2, -2, 2
					  , 0, 0.4, 0, 0.34, effort );

    System< EncT > henonEasy = henonMap( 1.375, 0.3
					 , -2, 2, -2, 1
					 , 0, 0.6, 0, 0.3, effort );
    System< EncT > henonHard = henonMap( 1.4, 0.3
					 , -3, 3, -3, 3
					 , 0, 0.25, 0, 0.25, effort );
    System< EncT > henonCrazy = henonMap( 1.425, 0.3
					  , -3, 3, -3, 3
					  , 0, 0.25, 0, 0.25, effort );
    
    System< EncT > tinkerbell1 = tinkerbell( 0.9, -0.6013, 2, 0.5
    					     , -2, 2, -2, 2
					     , -0.8, -0.6, -0.7, -0.5
    					     , effort );
    System< EncT > tinkerbell2 = tinkerbell( 0.3, 0.6, 2, 0.27
    					     , -0.5, 0.5, -0.5, 1.5
    					     , 0, 0.15, 0, 0.7, effort );

    System< EncT > duffing0 = duffing( 2.75, 0.2
    				       , -2, 2, -2, 2
    				       , -0.25, 1, -0.5, 0.8, effort );

    systems = {
	       // logisticFalse
	       // , logisticEasy
	       // , logisticHard
	       //, logisticCrazy
	        henonFalse
	       // , henonEasy
	       // , henonHard
	       //, henonCrazy
	       // , tinkerbell1
	       //, tinkerbell2
	       // , duffing0
    };
    
    const uint noSysConfigs = 4;

    // for( uint i = 0; i <= noSysConfigs; ++i )
    // {
    // 	// logistic map: values between rx in [0.9, 1.5], ry in [0.9, 1.1], positively correlated
    // 	systems.push_back( logisticMap( 0.9 + i * 0.6 / noSysConfigs, 0.9 + i * 0.2 / noSysConfigs
    // 					, 1.1, 0.5, 0.1, effort ) );
    // 	// henon map: a in [1.2, 1.5], b in [0.2, 0.4] positively correlated
    // 	systems.push_back( henonMap( 1.2 + i*0.3/noSysConfigs, 0.2 + i*0.2/noSysConfigs
    // 				     , 2.5, 0.25, effort ) );
    // 	// bogdanov map: eps in [0.0000005, 0.05], k in [0.000001, 0.1], mu in [0.0000001, 0.01]
    // 	systems.push_back( bogdanovMap( std::pow( 5, -1 - i), std::pow( 1, -i ), std::pow( 1, -1 - i )
    // 					, 1.25, 0.25, effort ) );
    // }

    // refiners
    LargestSideRefiner largestSide;
    // state heuristics
    RandomStateValue randomS;
    StateVolume volumeS;
    StateSideLength lengthS;
    StateVolumeDifference differenceS( 0.33 );
    // guides
    GreatestState greatestC;
    SumStates sumC;
    ExpSumStates expSumC;
    // termination
    LimitedTime< std::chrono::minutes > limitTime( std::chrono::minutes( 15 ) );

    const uint reps = 1;
    auto refiners = std::make_tuple( largestSide
				     );
    auto stateHs = std::make_tuple(
				   // randomS
				   volumeS
				   , lengthS
				   , differenceS
				     );
    auto cexHs = std::make_tuple( greatestC
				   , sumC
				   , expSumC
				   );

    uint noRefiners = tupleSize( refiners )
	, noStateHs = tupleSize( stateHs )
	, noCexHs = tupleSize( cexHs );
    
    std::cout << ",,"; // skip first two cols
    applyTuple( [&] (auto ref) { printHead( std::cout, ref, noStateHs * noCexHs ); }, refiners );
    std::cout << std::endl << ",,";
    for( uint i = 0; i < noRefiners; ++i )
	applyTuple( [&] (auto sh) { printHead( std::cout, sh, noCexHs ); }, stateHs );
    std::cout << std::endl << ",,";
    for( uint i = 0; i < noRefiners * noStateHs; ++i )
	applyTuple( [&] (auto ch) { printHead( std::cout, ch, 1 ); }, cexHs );
    std::cout << std::endl;
    
    for( auto& system : systems )
    {
	std::cout << system.name();
	Metrics mets;

	applyTuple( [&] (auto ref) {
		applyTuple( [&] (auto sh) {
			applyTuple( [&] (auto ch) {
					performWithConfig( system, reps, ref, sh, ch, limitTime, mets );
			    }
			    , cexHs );
		    }
		    , stateHs );
	    }
	    , refiners );
	std::cout << mets;
    }

}
