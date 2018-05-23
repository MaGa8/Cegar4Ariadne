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
  \interface for observers passed to cegar function
*/
template< typename Rtree, typename InitialIterT, typename CounterexampleIterT, typename RefinedIterT >
struct CegarObserver
{
    //! \brief called immediatly after the set of initial nodes has been determined, before the start of the loop
    virtual void initialized( const Rtree& rtree, const InitialIterT& iInitialNodesBegin, const InitialIterT& iInitialNodesEnd ) {}

    virtual void searchCounterexample();
    
    //! \brief called immediatly after a counterexample has been found
    virtual void foundCounterexample( const CounterexampleIterT& iCounterexampleBegin, const CounterexampleIterT& iCounterexampleEnd ) {}

    virtual void checkSpurious() {}
    
    virtual void spurious( const Ariadne::ValidatedUpperKleenean& spurious ) {}
    
    virtual void startRefinement( const typename Rtree::NodeT& toRefine ) {}
    
    //! \brief called immediatly after an abstraction has been refined
    virtual void refined( const Rtree& rtree, const RefinedIterT& iRefinedBegin, const RefinedIterT& iRefinedEnd ) {}

    //! \brief last statement called in loop
    virtual void finished() {}
};

// //! \class logs how much time is spent on 
// template< typename Rtree, typename IterT >
// class CegarTimer : public CegarObserver
// {
//   public:
//     typedef std::vector< typename std::chrono::high_resolution_clock::time_point > TimePointListT;
//     template< typename D >
//     using DurationListT = std::vector< D >;

//     template< typename D >
//     DurationList< D > toDuration( const TimePointListT& tps )
//     {
// 	DurationListT durations;
// 	durations.reserve( tps.size() );
// 	for( uint i = 0; i < tps.size() - 1; ++i )
// 	    durations.push_back( std::duration_cast< D >( tps[ i + 1 ] - tps[ i ] ) );
// 	return durations;
//     }

//     const TimePointListT& iterations() const { return mIterations; }

//     const TimePointListT& counterexamples() const { return mCounterexamples; }

//     const TimePointListT& refinements() const { return mRefinements; }

//     const TimePointListT& finishes() const { return mFinishes; }

//     void newIteration() { measurement( mIterations ); }

//     void counterexample( const IterT& iCounterexampleBegin, const IterT& iCounterexampleEnd ) { measurement( mCounterexamples ); }

//     virtual void spurious( const Ariadne::ValidatedUpperKleenean& spurious ) {}

//     void refinement( const Rtree& rtree ) { measurement( mRefinements ); }

//     void finishedIteration() { measurements( mFinishes ); }
    
//   private:
//     void measurement( const TimePointListT& l )
//     {
// 	l.push_back( std::chrono::high_resolution_clock::now() );
//     }
    
//     TimePointListT mIterations, mCounterexamples, mRefinements, mFinishes;  
// };


#endif
