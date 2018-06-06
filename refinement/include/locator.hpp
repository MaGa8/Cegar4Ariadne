#ifndef LOCATOR_HPP
#define LOCATOR_HPP

#include "refinementTree.hpp"

#include <vector>
#include <functional>
#include <algorithm>
#include <random>
#include <cmath>
#include <limits>

/* \brief a state heuristic supports 
   template< typename Rtree > double operator ()( const Rtree&, const typename Rtree::NodeT )
   mapping states to real values and
   std::string name() const 
   returing a describing name
*/

class RandomStateValue
{
  public:
    RandomStateValue()
	: mRandom( std::random_device()() )
	, mDist( std::numeric_limits< double >::min(), std::numeric_limits< double >::max() )
    {}

    template< typename Rtree, typename IterT >
    double operator ()( const Rtree& rtree, const IterT& cexBegin, const IterT& cexEnd, const IterT& istate )
    {
	return mDist( mRandom );
    }

    std::string name() const
    {
	return "random_state_value";
    }
    
  private:
    std::default_random_engine mRandom;
    std::uniform_real_distribution< double > mDist;
};

struct StateVolume
{
    template< typename Rtree, typename IterT >
    double operator ()( const Rtree& rtree, const IterT& cexBegin, const IterT& cexEnd, const IterT& istate )
    {
	auto nval = rtree.nodeValue( *istate );
	if( nval )
	    return nval.value().get().getEnclosure().measure().get_d();
	return 0;
    }

    std::string name() const
    {
	return "state_volume";
    }
};

struct StateSideLength
{
    template< typename Rtree, typename IterT >
    double operator ()( const Rtree& rtree, const IterT& cexBegin, const IterT& cexEnd, const IterT& istate )
    {
	auto nval = rtree.nodeValue( *istate );
	if( nval )
	    return nval.value().get().getEnclosure().radius().get_d();
	return 0;
    }

    std::string name() const
    {
	return "state_side_length";
    }
};

struct StateVolumeDifference
{
    //! \param endpointFactor factor to multiply volume with if no next state is available
    StateVolumeDifference( const double& endpointFactor ) : mEndpointFactor( endpointFactor ) {}
    
    template< typename Rtree, typename IterT >
    double operator ()( const Rtree& rtree, const IterT& cexBegin, const IterT& cexEnd, const IterT& istate ) const
    {
	IterT iNextState = istate; ++iNextState;
	auto sval = rtree.nodeValue( *istate );
	if( iNextState != cexEnd && sval )
	{
	    auto nextVal = rtree.nodeValue( *iNextState );
	    if( nextVal )
		return std::abs( (sval.value().get().getEnclosure().measure() - nextVal.value().get().getEnclosure().measure() ).get_d() );
	}

	if( sval ) // next state either not there or outside
	    return sval.value().get().getEnclosure().measure().get_d() * mEndpointFactor;
	return 0.0; // is outside
    }

    std::string name() const { return "state_volume_difference"; }

    const double mEndpointFactor;
};

/*
  a locator should implement 
  NodeRefVec operator ()( const RefinementTree< EnclosureT >& rtree, IteratorT ibegin, const IteratorT& iend )
  for NodeRefVec = std::vector< std::reference_wrapper< typename Rtree::NodeT > >

  constexpr string name() const;
 */

template< typename Rtree >
using NodeRefVec = std::vector< std::reference_wrapper< typename Rtree::NodeT > >;

//! \class represents policy of refining all nodes in a counterexample
struct CompleteCounterexample
{
    template< typename Rtree, typename IterT >
    NodeRefVec< Rtree > operator ()( const Rtree& rtree, IterT ibegin, const IterT& iend )
    {
	NodeRefVec< Rtree > ns;
	ns.reserve( std::distance( ibegin, iend ) );
	for( ; ibegin != iend; ++ibegin )
	    ns.push_back( std::ref( *ibegin ) );
	return ns;
    }

    std::string name() const {return "refine_all_states"; }
};

//! \brief selects random state to refine
class RandomStates
{
  public:
    RandomStates() : mRandom( std::random_device()() ) {}

    RandomStates( const RandomStates& orig ) = default;
    
    template< typename Rtree, typename IterT >
    NodeRefVec< Rtree > operator ()( const Rtree& rtree, IterT ibegin, const IterT& iend )
    {
	std::advance( ibegin, mDist( mRandom ) % std::distance( ibegin, iend ) );
	return { *ibegin };
    }

    std::string name() const {return "refine_random_states"; }
  private:
    std::uniform_int_distribution<> mDist;
    std::default_random_engine mRandom;
};

//! \brief selects the box with the largest volume
struct LargestBox
{
    template< typename IntervalT >
    using Rtree = RefinementTree< Ariadne::Box< IntervalT > >;

    template< typename IntervalT, typename IterT >
    NodeRefVec< Rtree< IntervalT > > operator ()( const Rtree< IntervalT >& rtree, IterT ibegin, const IterT& iend )
    {
	auto imax = std::max_element( ibegin, iend
				      , [&rtree] (const typename Rtree< IntervalT >::NodeT& n1
						  , const typename Rtree< IntervalT >::NodeT& n2) {
					  auto val1 = rtree.nodeValue( n1 ), val2 = rtree.nodeValue( n2 );
					  if( !val1 )
					      return true;
					  if( !val2 )
					      return false;
					  return definitely( val1.value().get().getEnclosure().measure() < val2.value().get().getEnclosure().measure() );
					  } );
	return { *imax };
    }
    
    std::string name() const {return "refine_largest_box"; }
};

template< size_t NSamples >
class MaximumEntropy
{
  public:
    template< typename IntervalT >
    using Rtree = RefinementTree< Ariadne::Box< IntervalT > >;

    MaximumEntropy()
	: mRandom( std::default_random_engine( std::random_device()() ) )
    {}

    template< typename IntervalT >
    typename Rtree< IntervalT >::EnclosureT randomPointBox( const Ariadne::Box< IntervalT >& bx )
    {
	auto pt = bx.lower_bounds() + mDist( mRandom ) * bx.upper_bounds();

	Ariadne::Array< Ariadne::ExactIntervalType > intervals( pt.size() );
	std::transform( pt.array().begin(), pt.array().end(), intervals.begin()
			, [] (const Ariadne::Approximation< Ariadne::FloatDP >& x) {
			    return Ariadne::ExactIntervalType( cast_exact( x ), cast_exact( x ) ); } );
	return typename Rtree< IntervalT >::EnclosureT( Ariadne::Vector( intervals ) );
    }

    template< typename IntervalT >
    double safeUnsafeScore( const Rtree< IntervalT >& rtree, const Ariadne::Box< IntervalT >& bx )
    {
	uint safe = 0, unsafe = 0;
	for( uint cSample = 0; cSample < NSamples; ++cSample )
	{
	    auto sample = randomPointBox( bx );
	    if( definitely( rtree.constraints().covers( Ariadne::Box< IntervalT >( sample ) ) ) )
		++safe;
	    else
		++unsafe;
	}
	double safeRatio = safe / static_cast< const double& >( NSamples )
	    , unsafeRatio = unsafe / static_cast< const double& >( NSamples );
	return -safeRatio*std::log( safeRatio ) - unsafeRatio*std::log( unsafeRatio );
    }

    template< typename IntervalT, typename IterT >
    NodeRefVec< Rtree< IntervalT > > operator ()( const Rtree< IntervalT >& rtree, IterT ibegin, const IterT& iend )
    {
	std::vector< double > scores( std::distance( ibegin, iend ) );
	std::transform( ibegin, iend, scores.begin()
			, [&rtree, this] (const typename Rtree< IntervalT >::NodeT& n) {
			    auto nval = rtree.nodeValue( n );
			    if( !nval )
				return 0.0;
			    return safeUnsafeScore( rtree, nval.value().get().getEnclosure() ); } );
	std::advance( ibegin, std::distance( scores.begin(), std::max_element( scores.begin(), scores.end() ) ) );
	return { *ibegin };
    }

    std::string name() const {return "refine_most_heterogenous_state"; }

  private:
    std::uniform_real_distribution<> mDist;
    std::default_random_engine mRandom;
};

#endif
