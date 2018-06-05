#ifndef COUNTEREXAMPLE_STORE_HPP
#define COUNTEREXAMPLE_STORE_HPP

#include "refinementTree.hpp"

#include <vector>
#include <map>
#include <iterator>

template< typename E >
using CounterexampleT = std::vector< RefinementTree< E >::NodeT >;

template< typename E >
struct ScoredCounterexample
{
    typedef std::map< RefinementTree< E >::NodeT, double, RefinementTree< E >::NodeComparator > StateValueMapT;

    template< typename IterT, typename SH >
    StateValueMapT assignValues( const RefinementTree& rtree, const IterT& cexBegin, const IterT& cexEnd, SH& sh )
    {
	StateValueMapT svmap( RefinementTree< E >::NodeComparator( rtree ) );
	std::transform( cexBegin, cexEnd, std::inserter( svmap ), svmap.end()
			, [] (auto& s) {return std::make_pair( s, sh( s ) ); } );
	return svmap;
    }
    
    template< typename IterT, typename SH, typename CH >
    ScoredCounterexample( const RefinementTree< E >& rtree, const IterT& cexBegin, const IterT& cexEnd, SH& sh, CH& ch )
	: mStates( assignValues( rtree, cexBegin, cexEnd, sh ) )
	, mTotal( std::accumulate( mStates.begin(), mStates.end(),
				   [&] (const double& val, auto& stateValue) { return ch( val, stateValue.second ); } ) )
    {}

    std::operator <( const ScoredCounterexample< E >& other )
    {
	return this->mTotal < other.mTotal;
    }
    
    const StateValueMapT mStates;
    const double mTotal;
};

/*!
  \class maintains a store of counterexample for efficient retrieval of most highly scoring counterexample
  \param SH state heuristic: assign score to state, higher scores meaning more suitable for refinement
  \param CH counterexample heuristic: maps states in counterexample and their score to single score
*/
template< typename SH, typename CH, typename E >
class CounterexampleStore
{
  public:
    CounterexampleStore( const SH& stateH, const CH& counterexampleH )
	: mCSet()
	, mStateH( stateH )
	, mCounterexampleH( counterexampleH )
    {}

    CounterexampleT< E > obtain()
    {
	CounterexampleT< E > ret;
	if( hasCounterexample() )
	{
	    ret.reserve( mCSet.begin()->mStates.size() );
	    std::transform( mCSet.begin()->mStates.begin(), mCSet.begin()->mStates.end(), std::back_inserter( ret )
			    , [] (auto& stateVal) {return stateVal.first;} );
	    mCSet.erase( mCSet.begin() );
	}
	return ret;
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
	mCset.emplace( rtree, beginCex, endCex, mStateH, mCounterexampleH );
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
