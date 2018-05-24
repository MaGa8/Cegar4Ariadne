#ifndef CEGAR_OBSERVER_HPP
#define CEGAR_OBSERVER_HPP

#include "refinementTree.hpp"

#include <chrono>
#include <vector>
#include <algorithm>
#include <iterator>
#include <functional>
#include <ostream>

/*!
  \class blueprint for observers passed to cegar function
  \note do not use as base class because member functions are not virtual
*/
// should template methods also be defined in header file rather than source file?
struct CegarObserver
{
    //! \brief called immediatly after the set of initial nodes has been determined, before the start of the loop
    template< typename Rtree >
    void initialized( const Rtree& rtree ) {}

    //! \brief first statement called in loop
    template< typename Rtree >
    void startIteration( const Rtree& rtree ) {}

    //! \brief immediatly called before the refinement tree is searched for counterexamples
    template< typename Rtree, typename IterT >
    void searchCounterexample( const Rtree& rtree, IterT& iAbstractionsBegin, IterT& iAbstractionsEnd ) {}
    
    //! \brief called immediatly after a counterexample has been found
    template< typename Rtree, typename IterT >
    void foundCounterexample( const Rtree& rtree, IterT& iCounterexBegin, IterT& iCounterexEnd ) {}

    //! \brief called immediatly before counterexample is checked
    template< typename Rtree, typename IterT >
    void checkSpurious( const Rtree& rtree, IterT& iCounterexBegin, IterT& iCounterexEnd ) {}

    //! \brief called immediatly after counterexample has been checked
    template< typename Rtree, typename IterT >
    void spurious( const Rtree& rtree, IterT& iCounterexBegin, IterT& iCounterexEnd, const Ariadne::ValidatedUpperKleenean& spurious ) {}

    //! \brief called immediatly before abstraction is refined
    template< typename Rtree >
    void startRefinement( const Rtree& rtree, const typename Rtree::NodeT& toRefine ) {}
    
    //! \brief called immediatly after an abstraction has been refined
    template< typename Rtree, typename IterT >
    void refined( const Rtree& rtree, const IterT& iRefinedBegin, const IterT& iRefinedEnd ) {}

    //! \brief last statement before return
    template< typename Rtree >
    void finished( const Rtree& rtree, const Ariadne::ValidatedKleenean safe ) {}
};

template< typename ObserverT, typename Rtree >
void callInitialized( ObserverT& obs, const Rtree& rtree )
{
    obs.initialized( rtree );
}


//! \class logs how much time is spent on
//! \param D duration type
template< typename D >
class CegarTimer : public CegarObserver
{
  public:
    unsigned long total() const { return mTotal; }
    unsigned long search() const { return mSearch; }
    unsigned long check() const { return mSpurious; }
    unsigned long refine() const { return mRefine; }
    unsigned long other() const { return mOther; }

    CegarTimer()
	: mTotal( 0 )
	, mSearch( 0 )
	, mSpurious( 0 )
	, mRefine( 0 )
	, mOther( 0 )
    {}

    //! \brief called immediatly after the set of initial nodes has been determined, before the start of the loop
    template< typename Rtree >
    void initialized( const Rtree& rtree )
    {
	mBeginLoop = std::chrono::high_resolution_clock::now();
	mLastMethodCall = mBeginLoop;
    }

    //! \brief immediatly called before the refinement tree is searched for counterexamples
    template< typename Rtree, typename IterT >
    void searchCounterexample( const Rtree& rtree, IterT& iAbstractionsBegin, IterT& iAbstractionsEnd )
    {
	measurement( mLastMethodCall, mOther );
    }
    
    //! \brief called immediatly after a counterexample has been found
    template< typename Rtree, typename IterT >
    void foundCounterexample( const Rtree& rtree, IterT& iCounterexBegin, IterT& iCounterexEnd )
    {
	measurement( mLastMethodCall, mSearch );
    }

    //! \brief called immediatly before counterexample is checked
    template< typename Rtree, typename IterT >
    void checkSpurious( const Rtree& rtree, IterT& iCounterexBegin, IterT& iCounterexEnd )
    {
	measurement( mLastMethodCall, mOther );
    }

    //! \brief called immediatly after counterexample has been checked
    template< typename Rtree, typename IterT >
    void spurious( const Rtree& rtree, IterT& iCounterexBegin, IterT& iCounterexEnd, const Ariadne::ValidatedUpperKleenean& spurious )
    {
	measurement( mLastMethodCall, mSpurious );
    }

    //! \brief called immediatly before abstraction is refined
    template< typename Rtree >
    void startRefinement( const Rtree& rtree, const typename Rtree::NodeT& toRefine )
    {
	measurement( mLastMethodCall, mOther );
    }
    
    //! \brief called immediatly after an abstraction has been refined
    template< typename Rtree, typename IterT >
    void refined( const Rtree& rtree, const IterT& iRefinedBegin, const IterT& iRefinedEnd )
    {
	measurement( mLastMethodCall, mRefine );
    }

    //! \brief last statement called in loop
    template< typename Rtree >
    void finished( const Rtree& rtree, const Ariadne::ValidatedKleenean& safe )
    {
	measurement( mBeginLoop, mTotal );
    }

  private:
    void measurement( std::chrono::high_resolution_clock::time_point& tpast, D& duration )
    {
	auto now = std::chrono::high_resolution_clock::now();
	duration += std::chrono::duration_cast< D >( now - tpast ).count();
	tpast = now;
    }
    
    std::chrono::high_resolution_clock::time_point mBeginLoop, mLastMethodCall;
    unsigned long mTotal, mSearch, mSpurious, mRefine, mOther;
};

class IterationCounter : public CegarObserver
{
  public:
    IterationCounter()
	: mIterations( 0 )
    {}

    uint iterations() const { return mIterations; }
    
    template< typename Rtree, typename IterT >
    void searchCounterexample( const Rtree& rtree, IterT& iAbstractionsBegin, IterT& iAbstractionsEnd )
    {
	++mIterations;
    }
  private:
    uint mIterations;
};


#endif
