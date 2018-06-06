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
	{
	    svs.push_back( std::make_pair( *istate, sh( rtree, cexBegin, cexEnd, istate ) ) );
	}
	

	return svs;
    }
    
    template< typename IterT, typename SH, typename CH >
    ScoredCounterexample( const RefinementTree< E >& rtree, const IterT& cexBegin, const IterT& cexEnd, SH& sh, CH& ch )
	: mStates( assignValues( rtree, cexBegin, cexEnd, sh ) )
	, mTotal( std::accumulate( mStates.begin(), mStates.end(), 0.0
				   , [&] (const double& val, auto& stateValue) { return ch( val, stateValue.second ); } ) )
    {  }

    bool operator <( const ScoredCounterexample< E >& other ) const
    {
	return this->mTotal < other.mTotal;
    }
    
    ScorePathT mStates;
    double mTotal;
};

template< typename E >
struct std::greater< ScoredCounterexample< E > >
{
    bool operator() (const ScoredCounterexample< E >& x1, const ScoredCounterexample< E >& x2 ) const
    {
	return !(x1 < x2);
    }
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

	std::pair< CounterexampleT< E >, typename RefinementTree< E >::NodeT > ret;
	ret.first.reserve( mCSet.begin()->mStates.size() );
	std::transform( mCSet.begin()->mStates.begin(), mCSet.begin()->mStates.end(), std::back_inserter( ret.first )
			, [] (auto& stateVal) {return stateVal.first;} );

	auto ispick = std::max_element( mCSet.begin()->mStates.begin(), mCSet.begin()->mStates.end()
					, [] (auto& sv1, auto& sv2) {return sv1.second < sv2.second;} );
	ret.second = ispick->first;

	mCSet.erase( mCSet.begin() );

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
	mCSet.emplace( rtree, beginCex, endCex, mStateH, mCounterexampleH );
    }

    //! \brief signals that n is being invalidated, thus no counterexample containing n should be handed out anymore
    void invalidate( const RefinementTree< E >& rtree, const typename RefinementTree< E >::NodeT& n )
    {
	typename RefinementTree< E >::NodeComparator ncomp( rtree );

	auto icset = mCSet.begin();
	while( icset != mCSet.end() )
	{
	    if( std::any_of( icset->mStates.begin(), icset->mStates.end(), [&] (auto& stateVal) {
				   return rtree.equal( n, stateVal.first ); } ) )
		icset = mCSet.erase( icset );
	    else
		++icset;
	}
    }

    //! \todo remove from interface once definitely not needed
    void outOfCounterexamples() {} // don't care for now

    void clear()
    {
	mCSet.clear();
    }

  private:
    typedef std::set< ScoredCounterexample< E >, std::greater< ScoredCounterexample< E > > > CSetT;

    CSetT mCSet;
    SH mStateH;
    CH mCounterexampleH;
};

#endif
