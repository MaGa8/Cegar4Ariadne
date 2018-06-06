#ifndef CEGAR_HPP
#define CEGAR_HPP

#include "refinementTree.hpp"
#include "refinement.hpp"
#include "locator.hpp"
#include "cegarObserver.hpp"
#include "termination.hpp"
#include "counterexampleStore.hpp"

#include "geometry/geometry.hpp"
#include "geometry/set.hpp"
#include "geometry/function_set.hpp"
#include "geometry/set_interface.hpp"

#include <stack>
#include <map>
#include <functional>

template< typename E >
using CounterexampleT = std::vector< typename RefinementTree< E >::NodeT >;

template< typename E >
using NodeSet = std::set< typename RefinementTree< E >::NodeT, typename RefinementTree< E >::NodeComparator >;

template< typename E >
using VisitMap = std::map< typename RefinementTree< E >::NodeT, bool, typename RefinementTree< E >::NodeComparator >;

/*!
  runs DFS to find counterexample
  any path terminates in
  1) loop leading back to state along path
  2) state with violated safety conditions
  \param iImgBegin iterator to beginning of refinement tree nodes describing the image of the initial set, should dereference to RefinementTree< E >::NodeT
  \return vector of nodes terminated by a possibly unsafe node
  \todo add parameter to control ordering of branches in dfs exploration 
  \todo remember which nodes were already explored & safe: if encountered again, no need to check further as it leads to known result!
*/
template< typename E, typename NodeIterT, typename SH, typename CH >
void findCounterexample( const RefinementTree< E >& rtree
			 , NodeIterT iImgBegin, NodeIterT iImgEnd
			 , CounterexampleStore< E, SH, CH >& counterStore
			 , VisitMap< E >& visitMap
			 , const std::vector< typename RefinementTree< E >::NodeT >& path = {}
			 )
{
    while( iImgBegin != iImgEnd && !counterStore.terminateSearch() )
    {
	auto iVisited = visitMap.find( *iImgBegin );
	
	if( iVisited == visitMap.end() || !iVisited->second )
	{
	    if( iVisited == visitMap.end() )
		visitMap.emplace( *iImgBegin, true );
	    else
		iVisited->second = true;

	    CounterexampleT< E > copyPath( path.begin(), path.end() );
	    copyPath.push_back( *iImgBegin );
	    // counterexample found (could not happen if node was visited before)
	    if( !definitely( rtree.isSafe( *iImgBegin ) ) )
		counterStore.found( rtree, copyPath.begin(), copyPath.end() );
	    else
	    {
		auto posts =  rtree.postimage( *iImgBegin );
		findCounterexample( rtree, posts.begin(), posts.end(), counterStore, visitMap, copyPath );
	    }
	}
	++iImgBegin;
    }
    counterStore.outOfCounterexamples();
}

template< typename E, typename NodeIterT, typename SH, typename CH >
void findCounterexample( RefinementTree< E >& rtree
			 , NodeIterT iAbstractionsBegin, NodeIterT iAbstractionsEnd
			 , CounterexampleStore< E, SH, CH >& counterStore
			 )
{
    VisitMap< E > visitMap( {}, typename RefinementTree< E >::NodeComparator( rtree ) );
    counterStore.startSearch();
    findCounterexample( rtree, iAbstractionsBegin, iAbstractionsEnd, counterStore, visitMap );
}

/*! 
  \param ibegin iterator over sequence of refinement tree nodes
  \return iterator to node pt lies in, according to eval
*/
template< typename E, typename NumberT, typename NodeIterT >
NodeIterT findContaining( const RefinementTree< E >& rtree, const Ariadne::Point< NumberT > pt
			  , NodeIterT ibegin, const NodeIterT& iend
			  , const std::function< bool( const Ariadne::ValidatedKleenean& ) >& eval )
{
    // center is not contained in initial image
    return std::find_if( ibegin, iend
			 ,[&] ( const typename RefinementTree< E >::NodeT& n )
			 {
			     auto val = rtree.nodeValue( n );
			     if( !val )
				 return false;			  
			     return eval( val.value().get().getEnclosure().contains( pt ) );
			 } );
}

template< typename F >
Ariadne::ExactBoxType boundsPoint2Box( const Ariadne::Point< Ariadne::Bounds< F > >& pt )
{
    Ariadne::Array< Ariadne::ExactIntervalType > intervals( pt.dimension() );
    std::transform( pt.array().begin(), pt.array().end(), intervals.begin()
		    , [] (const Ariadne::Bounds< F >& i) {
			return Ariadne::ExactIntervalType( cast_exact( i.lower() ), cast_exact( i.upper() ) ); } );
    return Ariadne::Box( Ariadne::Vector( intervals ) );
}

// implement this using lower kleenean?
/*! 
  \param beginCounter and endCounter iterators to beginning and end of counterexample trajectory, should dereference to typename RefinementTree< E >::NodeT
  \param beginImage and endImage iterators to beginning and end of image of initial set obtained from refinement tree, should dereference to typename RefinementTree< E >::NodeT as well
  \return false if there definitely exists a point that is mapped to the terminal state of the counterexample, indeterminate otherwise, including if there does not possibly exist such a point 
  why upper kleenean?
  if return false, know for sure that counterexample is not spurious because a point exist with trajectory leading to unsafe state
  if return true center point did not map along trajectory
  \todo allow divergence from supposed counterexample, i.e. follow trajectory of center point until loop
*/
template< typename E, typename PathIterT >
Ariadne::ValidatedUpperKleenean isSpurious( const RefinementTree< E >& rtree
					    , PathIterT beginCounter, PathIterT endCounter
					    , const Ariadne::BoundedConstraintSet& initialSet
					    , const Ariadne::Effort& effort )
{
    typedef RefinementTree< E > Rtree;
    
    // determine: initial set and first state of counterexample intersect
    std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > oBeginCex = rtree.nodeValue( *beginCounter );

    if( !oBeginCex )
    {
	// INCLUDE AGAIN ONCE ISSUE WITH INSIDE IS SETTLED
	// Ariadne::ValidatedKleenean initialInsideSafe = inside( initialSet, rtree.constraints(), 0.00001 );
	// if( definitely( !initialInsideSafe ) ) // replace this arbitrary number by parameter
	//    return false;

	return true;
    }

    // find correct typedef to work with
    Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > currPoint = oBeginCex.value().get().getEnclosure().centre();

    // std::function< bool( const typename Rtree::NodeT& ) > contains2 = [&currPoint, &rtree] (auto& n) {
    // 	std::function< Ariadne::ValidatedLowerKleenean( const typename Rtree::EnclosureT&, const Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > >& ) > contains =
    // 	  [] (auto& enc, auto& pt ) { return enc.contains( pt ); };
    // 	return possibly( rtree.relates( n, currPoint, contains ) ); };
    // if( std::none_of( beginImage, endImage, contains2 ) )
    // 	return true;

    Ariadne::ExactBoxType pointAsExactBox = boundsPoint2Box( currPoint );
    if( possibly( !initialSet.covers( pointAsExactBox ) ) )
	return true;
    
    //map forward
    const typename Rtree::EnclosureT& rtEnc = tree::value( rtree.tree(), tree::root( rtree.tree() ) )->getEnclosure();
    for( PathIterT nextCounter = beginCounter + 1; nextCounter != endCounter; beginCounter = nextCounter++ )
    {
	auto oNext = rtree.nodeValue( *nextCounter );
	Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > mappedPoint = rtree.dynamics().evaluate( currPoint );
	Ariadne::ValidatedKleenean containsMapped;

	if( oNext )
	    containsMapped = oNext.value().get().getEnclosure().contains( mappedPoint );
	else
	    containsMapped = !rtEnc.contains( mappedPoint );
	if( definitely( !containsMapped ) )
	    return true; // should be indeterminate
	currPoint = mappedPoint;
    }
    return false;
}

// can only prove that there exists a true counterexample -> system is unsafe
/*
  find counterexample: 
  (1) if eventually possibly unsafe -> refine                              not( definitely( covers( safeSet, bx ) )
  (2) if eventually definitely unsafe  and  not spurious -> return         definitely( separate( safeSet, bx ) )      
  (2) is subcase of (1)
*/
/*!
  \param rtree refinement tree to work on
  \param initialBegin begin of range of set of boxes describing the initial state
  \param effort effort to use for calculations
  \param refinementStrat strategy to use for refining individual box
  \param maxNodes number of nodes in tree after which to stop iterations
  \return pair of kleenean describing safety and sequence of nodes that forms a trajectory starting from the initial set
*/
template< typename E, typename RefinementT, typename SH, typename CH, typename TermT, typename ... ObserversT >
std::pair< Ariadne::ValidatedKleenean, CounterexampleT< E > > cegar( RefinementTree< E >& rtree
								     , const Ariadne::BoundedConstraintSet& initialSet
								     , const Ariadne::Effort& effort
								     , RefinementT refinement
								     , const SH& stateH
								     , const CH& counterexampleH
								     , TermT termination
								     , ObserversT& ... observers )
{
    typedef RefinementTree< E > Rtree;

    std::function< Ariadne::ValidatedUpperKleenean( const typename Rtree::EnclosureT&, const Ariadne::BoundedConstraintSet& ) > interPred =
	[effort] (auto& enc, auto& cset) {return !(cset.separated( enc ).check( effort ) ); };
    
    NodeSet< E > initialImage = NodeSet< E >( typename RefinementTree< E >::NodeComparator( rtree ) );
    {
	auto img = rtree.intersection( initialSet, interPred );
	initialImage.insert( img.begin(), img.end() );
    }
    CounterexampleStore< E, SH, CH > counters( stateH, counterexampleH );

    (callInitialized(observers, rtree), ...);

    termination.start( rtree );
    bool terminate = false;
    
    while( !terminate )
    {
	(callStartIteration( observers, rtree ), ... );
	(callSearchCounterexample(observers, rtree, initialImage.begin(), initialImage.end() ), ... );
	
	findCounterexample( rtree, initialImage.begin(), initialImage.end(), counters );

	(callSearchTerminated( observers, rtree ), ... );
	
	if( !counters.hasCounterexample() )
	{
	    (callFinished( observers, rtree, Ariadne::ValidatedKleenean( true ) ), ... );
	    return std::make_pair( Ariadne::ValidatedKleenean( true ), std::vector< typename Rtree::NodeT >() );
	}

	while( counters.hasCounterexample() && !terminate )
	{
	    auto counterexample = counters.obtain();
	    (callProcessCounterexample( observers, rtree, counterexample.first.begin(), counterexample.first.end() ), ... );
	    
	    (callCheckSpurious( observers, rtree, counterexample.first.begin(), counterexample.first.end() ), ... );
	
	    Ariadne::ValidatedUpperKleenean spurious = isSpurious( rtree, counterexample.first.begin(), counterexample.first.end(), initialSet, effort );

	    (callSpurious( observers, rtree, counterexample.first.begin(), counterexample.first.end(), spurious ), ... );
	
	    if( definitely( !spurious ) &&
		definitely( !rtree.isSafe( counterexample.first.back() ) ) )
	    {
		(callFinished( observers, rtree, Ariadne::ValidatedKleenean( false ) ), ... );
		return std::make_pair( Ariadne::ValidatedKleenean( false ), counterexample.first );
	    }

	    // \todo get rid of vector
	    std::vector< std::reference_wrapper< typename Rtree::NodeT > > nodesToRefine = { counterexample.second };

	    for( const typename Rtree::NodeT& refine : nodesToRefine )
	    {
		if( rtree.nodeValue( refine ) )
		{
		    const typename Rtree::RefinementT::NodeT& treeNodeRef =
			static_cast< const InsideGraphValue< typename Rtree::RefinementT::NodeT >& >( *graph::value( rtree.leafMapping(), refine ) ).treeNode();
		    auto iRefined = initialImage.find( refine );

		    (callStartRefinement( observers, rtree, refine ), ... );

		    counters.invalidate( rtree, refine );
		    rtree.refine( refine, refinement );

		    // find a way to prevent this statement from executing if no observers are passed
		    auto refinedNodes = rtree.leaves( treeNodeRef ); 
		    (callRefined( observers, rtree, refinedNodes.begin(), refinedNodes.end() ), ... );
	 
		    if( iRefined != initialImage.end() )
		    {
			initialImage.erase( iRefined );
			auto refinedInitials = rtree.intersection( treeNodeRef, initialSet, interPred );
			initialImage.insert( refinedInitials.begin(), refinedInitials.end() );
		    }
		}
	    }
	    terminate = termination( rtree ); // check in each inner loop for quick response
	}
    }
    (callFinished( observers, rtree, Ariadne::ValidatedKleenean( Ariadne::indeterminate ) ), ... );
    return make_pair( Ariadne::ValidatedKleenean( Ariadne::indeterminate ), std::vector< typename Rtree::NodeT >() );
}

#endif
