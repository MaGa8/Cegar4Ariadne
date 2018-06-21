#ifndef CEGAR_HPP
#define CEGAR_HPP

#include "refinementTree.hpp"
#include "refinement.hpp"
#include "locator.hpp"
#include "cegarObserver.hpp"
#include "termination.hpp"
#include "counterexampleStore.hpp"
#include "graphValuePrinter.hpp"

#include "geometry/geometry.hpp"
#include "geometry/set.hpp"
#include "geometry/function_set.hpp"
#include "geometry/set_interface.hpp"

#include <stack>
#include <map>
#include <unordered_set>
#include <list>
#include <functional>

#include <omp.h>

template< typename E >
using CounterexampleT = std::vector< typename RefinementTree< E >::NodeT >;

template< typename E >
using NodeSet = std::unordered_set< typename RefinementTree< E >::NodeT, NodeHash< E >, NodeEqual< E > >;

template< typename E >
using VisitMap = std::map< typename RefinementTree< E >::NodeT, bool, typename RefinementTree< E >::NodeComparator >;


/*!
  runs BFS to find counterexample
  any path terminates in
  1) loop leading back to state along path
  2) state with violated safety conditions
  \param iImgBegin iterator to beginning of refinement tree nodes describing the image of the initial set, should dereference to RefinementTree< E >::NodeT
  \return vector of nodes terminated by a possibly unsafe node
  \todo add parameter to control ordering of branches in dfs exploration 
  \todo remember which nodes were already explored & safe: if encountered again, no need to check further as it leads to known result!
*/
template< typename E, typename IterT, typename SH, typename CH >
void findCounterexample( const RefinementTree< E >& rtree
			 , const IterT& beginInitial, const IterT& endInitial
			 , CounterexampleStore< E, SH, CH >& cstore )
{
    std::vector< CounterexampleT< E > > paths, newPaths;
    NodeSet< E > visited( beginInitial, endInitial, graph::size( rtree.graph() ), NodeHash( rtree ), NodeEqual( rtree ) );
    std::transform( beginInitial, endInitial, std::back_inserter( paths )
		    , [] (auto& n) { return CounterexampleT< E >( {n} ); } );

    while( !paths.empty() && !cstore.terminateSearch() )
    {
#pragma omp parallel for
	for( uint np = 0; np < paths.size(); ++np )
	{
	    auto& boundaryNode = paths[ np ].back();

	    // std::cout << "exploreing counterexample of length " << ip->size() << std::endl;
	    // std::cout << "boundary " << GraphValuePrinter< E >( *graph::value( rtree.graph(), boundaryNode ), false, true, true, true ) << std::endl;
	    
	    if( possibly( !rtree.isSafe( boundaryNode ) ) )
	    {
#pragma omp critical
		cstore.found( rtree, paths[ np ].begin(), paths[ np ].end() );
	    }
	    else if( possibly( !rtree.isTransSafe( boundaryNode ) ) )
	    {
		auto img = rtree.postimage( boundaryNode );
		for( auto iImg = img.begin(); iImg != img.end(); ++iImg ) // try not to copy last path expanded
		{
		    if( visited.find( *iImg ) == visited.end() )
		    {
			CounterexampleT< E > copy( paths[ np ].begin(), paths[ np ].end() );
			copy.push_back( *iImg );
#pragma omp critical
			newPaths.push_back( copy );
		    }
#pragma omp critical		
		    visited.insert( *iImg );
		}
	    }
	}
	paths = std::move( newPaths );
	newPaths.clear();
    }
    cstore.outOfCounterexamples();
    // std::cout << "returning" << std::endl;
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

template< typename F >
Ariadne::ExactBoxType exactPoint2Box( const Ariadne::ExactPoint& pt )
{
    Ariadne::Array< Ariadne::ExactIntervalType > intervals( pt.dimension() );
    std::transform( pt.array().begin(), pt.array().end(), intervals.begin()
		    , [] (const F& x) {return Ariadne::ExactIntervalType( cast_exact( x ), cast_exact( x ) );} );
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
    auto oBeginCex = rtree.nodeValue( *beginCounter );

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
    const typename Rtree::EnclosureT& rtEnc = rtree.initialEnclosure();
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

/*!
  \brief test whether pt begins a trajectory inside sn that begins within the initial set and eventually leafes the safe set.
  \return true if pt certainly maps to unsafe state. false if abstract path corresponding to the trajectory of pt contains a loop before reaching an unsafe state
*/
template< typename E >
Ariadne::ValidatedUpperKleenean safeTrajectory( const RefinementTree< E >& rtree
						, const typename RefinementTree< E >::NodeT& sn
						, const Ariadne::ValidatedPoint& pt
						, const Ariadne::BoundedConstraintSet& initialSet
						, const Ariadne::Effort& effort )
{
    typedef RefinementTree< E > R;

    if( definitely( initialSet.separated( boundsPoint2Box( pt ) ) ) )
	return true;

    // return validated kleenean whether pt is in n
    auto ptInNode = [&rtree] (const Ariadne::ValidatedPoint& pt, const typename R::NodeT& n ) {
			auto nval = rtree.nodeValue( n );
			if( nval )
			    return nval.value().get().getEnclosure().contains( pt );
			return !rtree.initialEnclosure().contains( pt );
		    };
    
    typename R::NodeComparator ncomp( rtree );
    std::set< typename R::NodeT, typename R::NodeComparator > visited( ncomp );
    auto cs = sn;
    Ariadne::ValidatedPoint ptm = pt;
    while( visited.find( cs ) == visited.end() &&
    	   possibly( rtree.isSafe( cs ) ) &&
    	   possibly( ptInNode( ptm, cs ) ) ) // last clause: catch pt not in sn
    {
    	visited.insert( cs );
    	ptm = rtree.dynamics().evaluate( ptm );
    	auto img = rtree.postimage( cs ); 
    	auto iIn = std::find_if( img.begin(), img.end(), [&ptInNode, &ptm] (auto& n) { return possibly( ptInNode( ptm, n ) ); } );

    	if( iIn == img.end() )
    	    throw std::logic_error( "mapped point needs to map to one abstract state in abstract image" );
    	cs = *iIn;
    }

    return possibly( rtree.isSafe( cs ) );
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
    
    NodeSet< E > initialImage = NodeSet< E >( 0, NodeHash( rtree ), NodeEqual( rtree ) );
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
	
	    auto cexBeginVal = rtree.nodeValue( counterexample.first.front() );
	    Ariadne::ValidatedUpperKleenean safe = false;
	    if( cexBeginVal )
	    {
		Ariadne::Point< Ariadne::ValidatedNumericType > cen = cexBeginVal.value().get().getEnclosure().centre();
		safe = safeTrajectory( rtree, counterexample.first.front()
					 , cexBeginVal.value().get().getEnclosure().centre()
					 , initialSet, effort );
	    }

	    if( definitely( !safe ) )
	    {
		(callFinished( observers, rtree, Ariadne::ValidatedKleenean( false ) ), ... );
		return std::make_pair( Ariadne::ValidatedKleenean( false ), counterexample.first );
	    }

	    // \todo get rid of vector
	    std::vector< std::reference_wrapper< typename Rtree::NodeT > > nodesToRefine = { counterexample.second };

	    for( const typename Rtree::NodeT& refine : nodesToRefine )
	    {
		if( rtree.nodeValue( refine ) && possibly( rtree.isSafe( refine ) ) )
		{
		    auto iRefined = initialImage.find( refine );

		    (callStartRefinement( observers, rtree, refine ), ... );

		    counters.invalidate( rtree, refine );
		    auto refinedNodes = rtree.refine( refine, refinement );
		    
		    (callRefined( observers, rtree, refinedNodes.begin(), refinedNodes.end() ), ... );
	 
		    if( iRefined != initialImage.end() )
		    {
			initialImage.erase( iRefined );
			for( auto& nrefd : refinedNodes )
			{
			    if( possibly( rtree.overlapsConstraints( initialSet, nrefd ) ) )
				initialImage.insert( nrefd );
			}
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
