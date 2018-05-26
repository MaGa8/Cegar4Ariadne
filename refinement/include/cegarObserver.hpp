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
    void searchCounterexample( const Rtree& rtree, IterT iAbstractionsBegin, const IterT& iAbstractionsEnd ) {}
    
    //! \brief called immediatly after a counterexample has been found
    template< typename Rtree, typename IterT >
    void foundCounterexample( const Rtree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd ) {}

    //! \brief called immediatly before counterexample is checked
    template< typename Rtree, typename IterT >
    void checkSpurious( const Rtree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd ) {}

    //! \brief called immediatly after counterexample has been checked
    template< typename Rtree, typename IterT >
    void spurious( const Rtree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd, const Ariadne::ValidatedUpperKleenean& spurious ) {}

    //! \brief called immediatly before abstraction is refined
    template< typename Rtree >
    void startRefinement( const Rtree& rtree, const typename Rtree::NodeT& toRefine ) {}
    
    //! \brief called immediatly after an abstraction has been refined
    template< typename Rtree, typename IterT >
    void refined( const Rtree& rtree, IterT iRefinedBegin, const IterT& iRefinedEnd ) {}

    //! \brief last statement before return
    template< typename Rtree >
    void finished( const Rtree& rtree, const Ariadne::ValidatedKleenean safe ) {}
};

template< typename ObserverT, typename Rtree >
void callInitialized( ObserverT& obs, const Rtree& rtree )
{
    obs.initialized( rtree );
}

template< typename ObserverT, typename Rtree >
void callStartIteration( ObserverT& obs, const Rtree& rtree )
{
    obs.startIteration( rtree );
}

template< typename ObserverT, typename Rtree, typename IterT >
void callSearchCounterexample( ObserverT& obs, const Rtree& rtree, const IterT& iAbstractionsBegin, const IterT& iAbstractionsEnd )
{
    obs.searchCounterexample( rtree, iAbstractionsBegin, iAbstractionsEnd );
}
    
//! \brief called immediatly after a counterexample has been found
template< typename ObserverT, typename Rtree, typename IterT >
void callFoundCounterexample( ObserverT& obs, const Rtree& rtree, const IterT& iCounterexBegin, const IterT& iCounterexEnd )
{
    obs.foundCounterexample( rtree, iCounterexBegin, iCounterexEnd );
}

//! \brief called immediatly before counterexample is checked
template< typename ObserverT, typename Rtree, typename IterT >
void callCheckSpurious( ObserverT& obs, const Rtree& rtree, const IterT& iCounterexBegin, const IterT& iCounterexEnd )
{
    obs.checkSpurious( rtree, iCounterexBegin, iCounterexEnd );
}

//! \brief called immediatly after counterexample has been checked
template< typename ObserverT, typename Rtree, typename IterT >
void callSpurious( ObserverT& obs, const Rtree& rtree, const IterT& iCounterexBegin, const IterT& iCounterexEnd, const Ariadne::ValidatedUpperKleenean& spurious )
{
    obs.spurious( rtree, iCounterexBegin, iCounterexEnd, spurious );
}

//! \brief called immediatly before abstraction is refined
template< typename ObserverT, typename Rtree >
void callStartRefinement( ObserverT& obs, const Rtree& rtree, const typename Rtree::NodeT& toRefine )
{
    obs.startRefinement( rtree, toRefine );
}
    
//! \brief called immediatly after an abstraction has been refined
template< typename ObserverT, typename Rtree, typename IterT >
void callRefined( ObserverT& obs, const Rtree& rtree, const IterT& iRefinedBegin, const IterT& iRefinedEnd )
{
    obs.refined( rtree, iRefinedBegin, iRefinedEnd );
}

//! \brief last statement before return
template< typename ObserverT, typename Rtree >
void callFinished( ObserverT& obs, const Rtree& rtree, const Ariadne::ValidatedKleenean safe )
{
    obs.finished( rtree, safe );
}

//! \class logs how much time is spent on
//! \param D duration type
template< typename D >
class CegarTimer : public CegarObserver
{
  public:
    typedef std::chrono::high_resolution_clock ClockT;
    
    unsigned long total() const { return std::chrono::duration_cast< D >( mTotal ).count(); }
    unsigned long search() const { return std::chrono::duration_cast< D >( mSearch ).count(); }
    unsigned long check() const { return std::chrono::duration_cast< D >( mSpurious ).count(); }
    unsigned long refine() const { return std::chrono::duration_cast< D >( mRefine ).count(); }
    unsigned long other() const { return std::chrono::duration_cast< D >( mOther ).count(); }

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
    void searchCounterexample( const Rtree& rtree, IterT iAbstractionsBegin, const IterT& iAbstractionsEnd )
    {
	measurement( mLastMethodCall, mOther );
    }
    
    //! \brief called immediatly after a counterexample has been found
    template< typename Rtree, typename IterT >
    void foundCounterexample( const Rtree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd )
    {
	measurement( mLastMethodCall, mSearch );
    }

    //! \brief called immediatly before counterexample is checked
    template< typename Rtree, typename IterT >
    void checkSpurious( const Rtree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd )
    {
	measurement( mLastMethodCall, mOther );
    }

    //! \brief called immediatly after counterexample has been checked
    template< typename Rtree, typename IterT >
    void spurious( const Rtree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd, const Ariadne::ValidatedUpperKleenean& spurious )
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
    void refined( const Rtree& rtree, IterT iRefinedBegin, const IterT& iRefinedEnd )
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
    void measurement( typename ClockT::time_point& tpast, typename ClockT::duration& duration )
    {
	auto now = std::chrono::high_resolution_clock::now();
	duration += now - tpast;
	tpast = now;
    }
    
    ClockT::time_point mBeginLoop, mLastMethodCall;
    ClockT::duration mTotal, mSearch, mSpurious, mRefine, mOther;
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

struct DebugOutput
{
    //! \brief called immediatly after the set of initial nodes has been determined, before the start of the loop
    template< typename Rtree >
    void initialized( const Rtree& rtree )
    {
	std::cout << "cegar initialized" << std::endl;
    }

    //! \brief first statement called in loop
    template< typename Rtree >
    void startIteration( const Rtree& rtree )
    {
	std::cout << "start iteration: tree of size " << rtree.tree().size() << std::endl;
    }

    //! \brief immediatly called before the refinement tree is searched for counterexamples
    template< typename Rtree, typename IterT >
    void searchCounterexample( const Rtree& rtree, IterT iAbstractionsBegin, const IterT& iAbstractionsEnd )
    {
	std::cout << "start searching for counterexample" << std::endl;
    }
    
    //! \brief called immediatly after a counterexample has been found
    template< typename Rtree, typename IterT >
    void foundCounterexample( const Rtree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd )
    {
	std::cout << "found counterexample" << std::endl;
    }

    //! \brief called immediatly before counterexample is checked
    template< typename Rtree, typename IterT >
    void checkSpurious( const Rtree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd )
    {
	std::cout << "check if spurious" << std::endl;
    }

    //! \brief called immediatly after counterexample has been checked
    template< typename Rtree, typename IterT >
    void spurious( const Rtree& rtree, IterT iCounterexBegin, const IterT& iCounterexEnd, const Ariadne::ValidatedUpperKleenean& spurious )
    {
	std::cout << "determined if spurious " << std::endl;
    }

    //! \brief called immediatly before abstraction is refined
    template< typename Rtree >
    void startRefinement( const Rtree& rtree, const typename Rtree::NodeT& toRefine )
    {
	std::cout << "start refinement " << std::endl;
    }
    
    //! \brief called immediatly after an abstraction has been refined
    template< typename Rtree, typename IterT >
    void refined( const Rtree& rtree, IterT iRefinedBegin, const IterT& iRefinedEnd )
    {
	std::cout << "performed refinement " << std::endl;
    }

    //! \brief last statement before return
    template< typename Rtree >
    void finished( const Rtree& rtree, const Ariadne::ValidatedKleenean safe )
    {
	std::cout << "cegar terminated " << std::endl;
    }
};


#endif
