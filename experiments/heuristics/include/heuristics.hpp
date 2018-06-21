#ifndef HEURISTICS_HPP
#define HEURISTICS_HPP

#include "cegar.hpp"
#include "visualization.hpp"

#include "expression/constant.hpp"

#include <tuple>

template< typename E > class System;

/*! 
  \param R type of refinement
  \param L type of locator
*/
template< typename R, typename L, typename G, typename TermT >
struct ExperimentConfiguration
{
    ExperimentConfiguration( const uint& repetitions
			     , const R& refinement, const L& locator, const G& guide, const TermT& term )
	: mRepetitions( repetitions )
	, mRefinement( refinement )
	, mLocator( locator )
	, mGuide( guide )
	, mTerm( term )
    {}

    template< typename E, typename ... ObserversT >
    Ariadne::ValidatedKleenean run( const System< E >& sys, ObserversT& ... observers )
    {
	Ariadne::ValidatedKleenean safetyProof;
	for( uint i = 0; i < mRepetitions; ++i )
	{
	    RefinementTree< E > rtree = sys.refinementTree();
	    auto result = cegar( rtree, sys.initialSet(), rtree.effort(), mRefinement, mLocator, mGuide, mTerm, observers ... );
	    // auto mapVrange = graph::vertices( rtree.leafMapping() );
	    // std::cout << "nodes " << std::distance( mapVrange.first, mapVrange.second ) << " / " << rtree.tree().size()
		      // << "; system is safe: " << result.first << std::endl;
	    safetyProof = result.first;
	}
	return safetyProof;
    }
    
    const uint mRepetitions;
    const R mRefinement;
    const L mLocator;
    const G mGuide;
    const TermT mTerm;
};

template< typename E >
class System
{
  public:
    System( const Ariadne::BoundedConstraintSet& initial, Ariadne::BoundedConstraintSet& safe
	    , const Ariadne::EffectiveVectorFunction dynamics
	    , const Ariadne::Effort& effort, const std::string& name )
	: mInitial( initial ), mSafe( safe ), mDynamics( dynamics ), mEffort( effort ), mName( name )
    {}


    RefinementTree< E > refinementTree() const
    {
	return RefinementTree< E >( mSafe, mDynamics, mEffort );
    }

    const Ariadne::BoundedConstraintSet& initialSet() const
    {
	return mInitial;
    }

    std::string name() const {return mName; }

  private:
    Ariadne::BoundedConstraintSet mInitial, mSafe;
    Ariadne::EffectiveVectorFunction mDynamics;
    Ariadne::Effort mEffort;
    std::string mName;
};

template< typename T >
struct Metric : public std::pair< std::string, std::vector< T > >
{
    Metric( const std::string& name )
	: std::pair< std::string, std::vector< T > >( name, std::vector< T >() )
    {}

    void push_back( const T& val ) { this->second.push_back( val ); }

    typename std::vector< T >::iterator begin() { return this->second.begin(); }

    typename std::vector< T >::const_iterator begin() const { return this->second.begin(); }

    typename std::vector< T >::iterator end() { return this->second.end(); }

    typename std::vector< T >::const_iterator end() const { return this->second.end(); }
};

struct Metrics
{
    Metrics()
	: mTotal( "t_total" ), mSearch( "t_search" ), mCheck( "t_check" ), mRefine( "t_refine" ), mMisc( "t_misc" )
	, mIterations( "iterations" ), mAverage( "averages" ), mResult( "result" )
    {}
    
    Metric< unsigned long > mTotal, mSearch, mCheck, mRefine, mMisc;
    Metric< unsigned int > mIterations;
    Metric< double > mAverage;
    Metric< Ariadne::ValidatedKleenean > mResult;
};

std::ostream& operator <<( std::ostream& os, const Metrics& mets );

Ariadne::RealConstant make_constant( const std::string& name, const double& val );

// System< Ariadne::ExactBoxType > linearFlow( const Ariadne::RealConstant& iniw, const Ariadne::REalConstant& inih
					    // , const Ariadne::RealConstant& );

System< Ariadne::ExactBoxType > logisticMap( const double& rx, const double& ry
					  , const double& wsafe, const double& winitial
					  , const double& delta, const Ariadne::Effort& effort );

System< Ariadne::ExactBoxType > henonMap( const double& a, const double& b
				       , const double& wsafe, const double& winitial, const Ariadne::Effort& effort );

System< Ariadne::ExactBoxType > bogdanovMap( const double& epsilon, const double& k, const double& mu
					      , const double& wsafe, const double winitial, const Ariadne::Effort& effort );


template< typename CallableT, typename TupleT, size_t ... S >
decltype( auto ) applyTupleImpl( CallableT f, TupleT t, std::index_sequence< S ... > )
{
    return (f( std::get< S >( t ) ), ... );
}

template< typename CallableT, typename TupleT >
decltype( auto ) applyTuple( CallableT f, TupleT t )
{
    return applyTupleImpl( f, t, std::make_index_sequence< std::tuple_size< TupleT >::value >() );
}

template< typename Tpl >
constexpr size_t tupleSize( const Tpl& t )
{
    return std::tuple_size< Tpl >::value;
}

template< typename H >
void printHead( std::ostream& o, const H& heuristic, const uint padding )
{
    o << heuristic.name();
    for( uint i = 0; i < padding; ++i )
	o << ",";
}

template< typename E, typename R, typename L, typename G, typename TermT, typename ... Ms >
void performWithConfig( const System< E >& system, const uint& repetitions
			, const R& refinement, const L& locator, const G& guide, const TermT& term, Metrics& mets )
{
    ExperimentConfiguration< R, L, G, TermT > config( repetitions, refinement, locator, guide, term );

    CegarTimer< std::chrono::milliseconds > timer;
    IterationCounter iterations;
    CounterexampleLengthAverage cexLength;
    // DebugOutput debugOut;

    try
    {
	auto safety = config.run( system, timer, iterations, cexLength );

	mets.mTotal.push_back( timer.total() );
	mets.mSearch.push_back( timer.search() );
	mets.mCheck.push_back( timer.check() );
	mets.mRefine.push_back( timer.refine() );
	mets.mMisc.push_back( timer.other() );
	mets.mIterations.push_back( iterations.iterations() );
	mets.mAverage.push_back( cexLength.average() );
	mets.mResult.push_back( safety );

	// std::cout << system.name() << " + " << refinement.name() << " + " << locator.name() << " + " << guide.name() << std::endl
	// 	  << "total " << timer.total()
	// 	  << ", searching " << timer.search()
	// 	  << ", checking " << timer.check()
	// 	  << ", refining " << timer.refine()
	// 	  << ", misc " << timer.other()
	// 	  << ", iterations " << iterations.iterations()
	// 	  << std::endl;
	
	// maybe do a visualization as well?
    }
    catch( const std::exception& e )
    {
	std::cout << "failed because " << e.what() << std::endl;
    }
    // std::cout << std::endl;
}

#endif

