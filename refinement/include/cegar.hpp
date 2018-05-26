#ifndef CEGAR_HPP
#define CEGAR_HPP

#include "refinementTree.hpp"
#include "refinement.hpp"
#include "locator.hpp"
#include "cegarObserver.hpp"

#include <stack>
#include <functional>

template< typename E >
class NodeComparator
{
  public:
    NodeComparator( const RefinementTree< E >& rtree ) : mRtree( rtree ) {}
		       
    bool operator ()( const typename RefinementTree< E >::NodeT& n1
		      , const typename RefinementTree< E >::NodeT& n2 ) const
    {
	std::optional< std::reference_wrapper< const InteriorTreeValue< E > > > otval1 = mRtree.get().nodeValue( n1 )
	    , otval2 = mRtree.get().nodeValue( n2 );
	// always unsafe node is always equal to always unsafe node
	if( !otval1 && !otval2 )
	    return true;
	if( !otval1 || !otval2 )
	    return false;
	return otval1.value().get().id() < otval2.value().get().id();
    }

  private:
    std::reference_wrapper< const RefinementTree< E > > mRtree;
};

template< typename E >
using CounterexampleT = std::vector< typename RefinementTree< E >::NodeT >;

template< typename E >
using NodeSet = std::set< typename RefinementTree< E >::NodeT, NodeComparator< typename RefinementTree< E >::EnclosureT > >;

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
template< typename E, typename NodeIterT, typename GuideT >
void findCounterexample( RefinementTree< E >& rtree
					 , NodeIterT iImgBegin, NodeIterT iImgEnd
					 , GuideT& guide
					 , const std::vector< typename RefinementTree< E >::NodeT >& path = {} )
{
    typedef RefinementTree< E > Rtree;

    // std::stack< std::pair< NodeIterT, NodeIterT > > iterStack;
    // while( !iterStack.empty() && !guide.terminateSearch()

    while( iImgBegin != iImgEnd && !guide.terminateSearch() )
    {
	auto iLoop = std::find_if( path.begin(), path.end()
				   , std::bind( &RefinementTree< E >::equal, &rtree, *iImgBegin, std::placeholders::_1 ) );
	if( iLoop == path.end() )
	{
	    CounterexampleT< E > copyPath( path.begin(), path.end() );
	    copyPath.push_back( *iImgBegin );
	    // counterexample found (could not happen if node was visited before)
	    if( !definitely( rtree.isSafe( *iImgBegin ) ) )
	    {
		guide.found( copyPath.begin(), copyPath.end() );
		// return copyPath;
	    }
	    else
	    {
		auto posts =  rtree.postimage( *iImgBegin );
		findCounterexample( rtree, posts.begin(), posts.end(), guide, copyPath );
	    }
	}
	++iImgBegin;
    }
    guide.outOfCounterexamples();
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
					    , const Ariadne::ConstraintSet& initialSet
					    , const Ariadne::Effort& effort )
{
    typedef RefinementTree< E > Rtree;
    
    // determine: initial set and first state of counterexample intersect
    std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > oBeginCex = rtree.nodeValue( *beginCounter );

    if( !oBeginCex )
    {
	std::function< Ariadne::ValidatedLowerKleenean( const typename Rtree::EnclosureT&, const Ariadne::ConstraintSet& ) > covers =
	    [&effort] (auto& nodeEnc, auto& cset ) {return cset.covers( nodeEnc ).check( effort ); };

	if( definitely( rtree.relates( rtree.outside(), initialSet, covers ) ) ) // passed in outside: tests whether initial box is not covered
	    return false;
	// if( std::any_of( beginImage, endImage, [&rtree, &notCovers] (const typename Rtree::NodeT& n) {
		    // return possibly( rtree.relates( n, rtree.initialEnclosure(), notCovers ) ); } ) ) 
	    // return false;
	else
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
template< typename E, typename RefinementT, typename LocatorT, template< typename EE > typename GuideT, typename ... ObserversT >
std::pair< Ariadne::ValidatedKleenean, CounterexampleT< E > > cegar( RefinementTree< E >& rtree
								     , const Ariadne::ConstraintSet& initialSet
								     , const Ariadne::Effort& effort
								     , const RefinementT& refinement
								     , LocatorT locator
								     , GuideT< E > guide
								     , const uint maxNodes
								     , ObserversT& ... observers )
{
    typedef RefinementTree< E > Rtree;

    std::function< Ariadne::ValidatedLowerKleenean( const typename Rtree::EnclosureT&, const Ariadne::ConstraintSet& ) > interPred =
	[effort] (auto& enc, auto& cset) {return cset.overlaps( enc ).check( effort ); };
    
    NodeSet< E > initialImage = NodeSet< E >( NodeComparator( rtree ) );
    {
	auto img = rtree.intersection( initialSet, interPred );
	initialImage.insert( img.begin(), img.end() );
    }

    (callInitialized(observers, rtree), ...);
        
    while( rtree.tree().size() < maxNodes )
    {
	(callStartIteration( observers, rtree ), ... );
	(callSearchCounterexample(observers, rtree, initialImage.begin(), initialImage.end() ), ... );
	
	// look for counterexample
	guide.startSearch(); // move this to find counterexample once search is no longer recursive
	findCounterexample( rtree, initialImage.begin(), initialImage.end(), guide );

	// found possibly many counterexamples
	(callSearchTerminated( observers, rtree ), ... );
	
	if( !guide.hasCounterexample() )
	{
	    (callFinished( observers, rtree, Ariadne::ValidatedKleenean( true ) ), ... );
	    return std::make_pair( Ariadne::ValidatedKleenean( true ), std::vector< typename Rtree::NodeT >() );
	}

	while( guide.hasCounterexample() )
	{
	    auto counterexample = guide.obtain();
	    (callProcessCounterexample( observers, rtree, counterexample.begin(), counterexample.end() ), ... );
	    
	    (callCheckSpurious( observers, rtree, counterexample.begin(), counterexample.end() ), ... );
	
	    Ariadne::ValidatedUpperKleenean spurious = isSpurious( rtree, counterexample.begin(), counterexample.end(), initialSet, effort );

	    (callSpurious( observers, rtree, counterexample.begin(), counterexample.end(), spurious ), ... );
	
	    if( definitely( !spurious ) &&
		definitely( !rtree.isSafe( counterexample.back() ) ) )
	    {
		(callFinished( observers, rtree, Ariadne::ValidatedKleenean( false ) ), ... );
		return std::make_pair( Ariadne::ValidatedKleenean( false ), counterexample );
	    }

	    std::vector< std::reference_wrapper< typename Rtree::NodeT > > nodesToRefine = locator( rtree, counterexample.begin(), counterexample.end() );
	    for( const typename Rtree::NodeT& refine : nodesToRefine )
	    {
		if( rtree.nodeValue( refine ) )
		{
		    const typename Rtree::RefinementT::NodeT& treeNodeRef =
			static_cast< const InsideGraphValue< typename Rtree::RefinementT::NodeT >& >( *graph::value( rtree.leafMapping(), refine ) ).treeNode();
		    auto iRefined = initialImage.find( refine );

		    (callStartRefinement( observers, rtree, refine ), ... );
	
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
	}
    }
    (callFinished( observers, rtree, Ariadne::ValidatedKleenean( Ariadne::indeterminate ) ), ... );
    return make_pair( Ariadne::ValidatedKleenean( Ariadne::indeterminate ), std::vector< typename Rtree::NodeT >() );
}

#endif
