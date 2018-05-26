#ifndef GUIDE_HPP
#define GUIDE_HPP

#include "refinementTree.hpp"

#include <vector>
#include <cmath>

template< typename E >
using CounterexampleT = std::vector< typename RefinementTree< E >::NodeT >;

/* 
   guide class should provide
   
   template< typename IterT > void found( IterT beginCounterex, const IterT& endCounterex )
   bool terminateSearch()
   void outOfCounterexamples()
   bool hasCounterexample()
   CounterexampleT< E > obtain()
*/

/*!
  \class for each counterexample found: randomly decide whether or not to replace the one previously stored
  keep with probability p^k where 0 < p < 1 and k is the number of counterexamples found so far

 */
template< typename E >
class KeepRandomCounterexamples
{
  public:
    KeepRandomCounterexamples( double p )
	: mP( p )
	, mDist( 0.0, 1.0 )
	, mRandom( std::random_device()() )
	, mSeen( 0 )
    {}
    
    template< typename IterT >
    void found( const IterT& beginCounterex, const IterT& endCounterex )
    {
	// if search was terminated and new search begins
	if( mSeen == 0 )
	    mCounterexample.clear();
	if( mDist( mRandom ) < std::pow( mP, mSeen ) )
	    mCounterexample = CounterexampleT< E >( beginCounterex, endCounterex );
	++mSeen;
    }

    void outOfCounterexamples()
    {
	mSeen = 0;
    }

    bool terminateSearch() { return mSeen > 500; }

    bool hasCounterexample() { return !mCounterexample.empty(); }

    CounterexampleT< E > obtain() // move constructor is invoked, right?
    {
	auto tmp = mCounterexample;
	mCounterexample.clear();
	return tmp;
    }

  private:
    double mP;
    CounterexampleT< E > mCounterexample;
    std::uniform_real_distribution<> mDist;
    std::default_random_engine mRandom;
    uint mSeen;
};


#endif
