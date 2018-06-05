#ifndef COUNTEREXAMPLE_STORE_HPP
#define COUNTEREXAMPLE_STORE_HPP

#include "refinementTree.hpp"

#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

template< typename E >
using CounterexampleT = std::vector< typename RefinementTree< E >::NodeT >;

template< typename E >
struct ScoredCounterexample
{
    typedef std::vector< std::pair< typename RefinementTree< E >::NodeT, double > > ScorePathT;

    template< typename IterT, typename SH >
    static ScorePathT assignValues( const RefinementTree< E >& rtree, const IterT& cexBegin, const IterT& cexEnd, SH& sh )
    {
	ScorePathT svs;
	svs.reserve( std::distance( cexBegin, cexEnd ) );

	for( auto istate = cexBegin; istate != cexEnd; ++istate )
	    svs.push_back( std::make_pair( *istate, sh( rtree, cexBegin, cexEnd, istate ) ) ); // insertion intended

	return svs;
    }
    
    template< typename IterT, typename SH, typename CH >
    ScoredCounterexample( const RefinementTree< E >& rtree, const IterT& cexBegin, const IterT& cexEnd, SH& sh, CH& ch )
	: mStates( assignValues( rtree, cexBegin, cexEnd, sh ) )
	, mTotal( std::accumulate( mStates.begin(), mStates.end(), 0.0
				   , [&] (const double& val, auto& stateValue) { return ch( val, stateValue.second ); } ) )
    {}

    bool operator <( const ScoredCounterexample< E >& other ) const
    {
	return this->mTotal > other.mTotal;
    }
    
    ScorePathT mStates;
    double mTotal;
};

/*!
  \class maintains a store of counterexample for efficient retrieval of most highly scoring counterexample
  \param SH state heuristic: assign score to state, higher scores meaning more suitable for refinement
  \param CH counterexample heuristic: maps states in counterexample and their score to single score
*/
template< typename E, typename SH, typename CH >
class CounterexampleStore
{
  public:
    CounterexampleStore( const SH& stateH, const CH& counterexampleH )
	: mCSet()
	, mStateH( stateH )
	, mCounterexampleH( counterexampleH )
    {}

    std::pair< CounterexampleT< E >, typename RefinementTree< E >::NodeT > obtain()
    {
	if( !hasCounterexample() )
	    throw std::runtime_error( "counterexample store does not hold counterexamples but is asked for one" );
	
	CounterexampleT< E > ret;
	ret.reserve( mCSet.begin()->mStates.size() );
	std::transform( mCSet.begin()->mStates.begin(), mCSet.begin()->mStates.end(), std::back_inserter( ret )
			, [] (auto& stateVal) {return stateVal.first;} );
	auto maxStateVal = std::max_element( mCSet.begin()->mStates.begin(), mCSet.begin()->mStates.end()
					     , [] (auto& sv1, auto& sv2) {return sv1.second > sv2.second;} );

	mCSet.erase( mCSet.begin() );

	return std::make_pair( ret, maxStateVal->first );
    }
    
    bool hasCounterexample() const
    {
	return !mCSet.empty();
    }

    //! \todo find condition on which to terminate (find counterexample with value as high as last one?)
    bool terminateSearch() const
    {
	return false; // never terminate
    }

    //! \todo abandon in favor of clean
    void startSearch()
    {
	mCSet.clear();
    }

    template< typename IterT >
    void found( const RefinementTree< E >& rtree, const IterT& beginCex, const IterT& endCex )
    {
	mCSet.emplace( rtree, beginCex, endCex, mStateH, mCounterexampleH );
    }

    //! \todo remove from interface once definitely not needed
    void outOfCounterexamples() {} // don't care now

    void clear()
    {
	mCSet.clear();
    }

  private:
    typedef std::set< ScoredCounterexample< E > > CSetT;

    CSetT mCSet;
    SH mStateH;
    CH mCounterexampleH;
};

#endif
