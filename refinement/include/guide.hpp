#ifndef GUIDE_HPP
#define GUIDE_HPP

#include "refinementTree.hpp"

#include <vector>
#include <cmath>

template< typename E >
using CounterexampleT = std::vector< typename RefinementTree< E >::NodeT >;

/* 
   guide class should provide
   
   void startSearch
   template< typename IterT > void found( IterT beginCounterex, const IterT& endCounterex )
   bool terminateSearch()
   void outOfCounterexamples()
   bool hasCounterexample()
   CounterexampleT< E > obtain()
*/
template< typename E >
struct Guide
{
    //! \brief called before the search is started
    void startSearch() {}

    //! \brief called whenever a counterexample is found
    template< typename IterT >
    void found( IterT beginCounterex, const IterT& encCounterex ) {}

    //! \return true if the search should be terminated
    bool terminateSearch() {}

    //! \brief called if all states have been visited
    void outOfCounterexamples() {}

    //! \return true if guide has at least one counterexample
    bool hasCounterexample() {}

    //! \return some counterexample
    CounterexampleT< E > obtain() {}

    std::string name() const {return "do_nothing"; }
};

/*!
  \class for each counterexample found: randomly decide whether or not to replace the one previously stored
  replace with probability p
  chance of i th counterexample to remain: p^(n - i)
  to balance: stop search with probability p at each step
 */
template< typename E >
class KeepRandomCounterexamples : public Guide< E >
{
  public:
    KeepRandomCounterexamples( double p )
	: mP( p )
	, mDist( 0.0, 1.0 )
	, mRandom( std::random_device()() )
	, mTerminate( false )
    {}

    void startSearch()
    {
	mTerminate = false;
	mCounterexample.clear();
    }
    
    template< typename IterT >
    void found( const IterT& beginCounterex, const IterT& endCounterex )
    {
	// if search was terminated and new search begins
	if( mCounterexample.empty() || mDist( mRandom ) < mP )
	    mCounterexample = CounterexampleT< E >( beginCounterex, endCounterex );

	if( mDist( mRandom ) < mP )
	    mTerminate = true;
    }

    bool terminateSearch() { return mTerminate; }

    bool hasCounterexample() { return !mCounterexample.empty(); }

    CounterexampleT< E > obtain() // move constructor is invoked, right?
    {
	auto tmp = mCounterexample;
	mCounterexample.clear();
	return tmp;
    }

    std::string name() const {return "process_random_counterexample";}

  private:
    double mP;
    CounterexampleT< E > mCounterexample;
    std::uniform_real_distribution<> mDist;
    std::default_random_engine mRandom;
    bool mTerminate;
};


#endif
