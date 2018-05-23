#ifndef LOCATOR_HPP
#define LOCATOR_HPP

#include "refinementTree.hpp"

#include <vector>
#include <functional>
#include <algorithm>
#include <random>
#include <cmath>

#ifndef LOC_DEFS
#define LOC_DEFS typedef typename BaseT::EnclosureT EnclosureT; typedef typename BaseT::IteratorT IteratorT; typedef typename BaseT::NodeRefVec NodeRefVec
#endif

/*!
  \interface for locator function objects which shall be used to determine the state to refine
  \param EnclosureT type of enclosure to refine
  cannot add method to interface because cannot combine virtual methods and template methods
*/
template< typename E, typename I >
struct ILocator
{
    typedef E EnclosureT;
    typedef I IteratorT;
    typedef std::vector< std::reference_wrapper< typename RefinementTree< EnclosureT >::NodeT > > NodeRefVec;

    virtual NodeRefVec operator ()( const RefinementTree< EnclosureT >& rtree, IteratorT ibegin, const IteratorT& iend ) = 0;
};

//! \class represents policy of refining all nodes in a counterexample
template< typename E, typename I >
struct CompleteCounterexample : public ILocator< E, I >
{
    typedef ILocator< E, I > BaseT;
    LOC_DEFS;
    
    NodeRefVec operator ()( const RefinementTree< EnclosureT >& rtree, IteratorT ibegin, const IteratorT& iend )
    {
	NodeRefVec ns;
	ns.reserve( std::distance( ibegin, iend ) );
	for( ; ibegin != iend; ++ibegin )
	    ns.push_back( std::ref( *ibegin ) );
	return ns;
    }
};

//! \brief selects random state to refine
template< typename E, typename I >
class RandomStates : public ILocator< E, I >
{
  public:
    typedef ILocator< E, I > BaseT;
    LOC_DEFS;
    
    RandomStates() : mRandom( std::random_device()() ) {}

    RandomStates( const RandomStates& orig ) = default;
    
    NodeRefVec operator ()( const RefinementTree< EnclosureT >& rtree, IteratorT ibegin, const IteratorT& iend )
    {
	std::advance( ibegin, mDist( mRandom ) % std::distance( ibegin, iend ) );
	return { *ibegin };
    }
  private:
    std::uniform_int_distribution<> mDist;
    std::default_random_engine mRandom;
};

//! \brief selects the box with the largest volume
template< typename IntervalT, typename IterT >
struct LargestBox : public ILocator< Ariadne::Box< IntervalT >, IterT >
{
    typedef ILocator< Ariadne::Box< IntervalT >, IterT > BaseT;
    LOC_DEFS;
    
    NodeRefVec operator ()( const RefinementTree< EnclosureT >& rtree, IteratorT ibegin, const IteratorT& iend ) const
    {
	auto imax = std::max_element( ibegin, iend
				      , [&rtree] (const typename RefinementTree< EnclosureT >::NodeT& n1
						  , const typename RefinementTree< EnclosureT >::NodeT& n2) {
					  auto val1 = rtree.nodeValue( n1 ), val2 = rtree.nodeValue( n2 );
					  if( !n1 )
					      return true;
					  if( !n2 )
					      return false;
					  return val1.value().get().getEnclosure().measure() < val2.value().get().getEnclosure().measure();
					  } );
	return { *imax };
    }
};

template< typename IntervalT, typename IterT >
class MaximumEntropy : public ILocator< Ariadne::Box< IntervalT >, IterT >
{
  public:
    typedef ILocator< Ariadne::Box< IntervalT >, IterT > BaseT;
    LOC_DEFS;

    MaximumEntropy( const uint& noSamples )
	: mNoSamples( noSamples )
	, mRandom( std::default_random_engine( std::random_device()() ) )
    {}

    auto randomPoint( const EnclosureT& bx ) -> decltype( bx.lower_bounds() + bx.upper_bounds() )
    {
	return bx.lower_bounds() + mDist( mRandom ) * bx.upper_bounds();
    }

    double safeUnsafeScore( const RefinementTree< EnclosureT >& rtree, const EnclosureT& bx )
    {
	uint safe = 0, unsafe = 0;
	for( uint cSample = 0; cSample < mNoSamples; ++cSample )
	{
	    auto sample = randomPoint( bx );
	    if( definitely( rtree.constraints().covers( Ariadne::Box< IntervalT >( sample ) ) ) )
		++safe;
	    else
		++unsafe;
	}
	double safeRatio = safe / static_cast< const double& >( mNoSamples )
	    , unsafeRatio = unsafe / static_cast< const double& >( mNoSamples );
	return -safeRatio*std::log( safeRatio ) - unsafeRatio*std::log( unsafeRatio );
    }
    
    NodeRefVec operator ()( const RefinementTree< EnclosureT >& rtree, IteratorT& ibegin, const IteratorT& iend )
    {
	std::vector< double > scores( std::distance( ibegin, iend ) );
	std::transform( ibegin, iend, scores.begin()
			, [&rtree, this] (const typename RefinementTree< EnclosureT >::NodeT& n) {
			    return safeUnsafeScore( rtree, n ); } );
	return { *std::advance( ibegin, std::distance( ibegin
						       , std::max_element( scores.begin(), scores.end() ) ) ) };
    }

  private:
    uint mNoSamples;
    std::uniform_real_distribution<> mDist;
    std::default_random_engine mRandom;
};

#undef LOC_DEFS

#endif
