#ifndef CEGAR_HPP
#define CEGAR_HPP

#include "refinementTree.hpp"

#include <functional>

/*!
  runs DFS to find counterexample
  any path terminates in
  1) loop leading back to state along path
  2) state with violated safety conditions
  \param iImgBegin iterator to beginning of refinement tree nodes describing the image of the initial set, should dereference to RefinementTree< IntervalT >::NodeT
  \return vector of nodes terminated by a possibly unsafe node
  \todo add parameter to control ordering of branches in dfs exploration 
  \todo remember which nodes were already explored & safe: if encountered again, no need to check further as it leads to known result!
*/
template< typename IntervalT, typename NodeIterT >
std::vector< typename RefinementTree< IntervalT >::NodeT > findCounterexample( RefinementTree< IntervalT >& rtree
									       , NodeIterT iImgBegin, NodeIterT iImgEnd
									       , const std::vector< typename RefinementTree< IntervalT >::NodeT >& path = {} )
{
    typedef RefinementTree< IntervalT > Rtree;
    for( ; iImgBegin != iImgEnd; ++iImgBegin )
    {
	auto iLoop = std::find_if( path.begin(), path.end()
				   , std::bind( &RefinementTree< IntervalT >::equal, &rtree, *iImgBegin, std::placeholders::_1 ) );
	if( iLoop == path.end() )
	{
	    std::vector< typename RefinementTree< IntervalT >::NodeT > copyPath( path.begin(), path.end() );
	    copyPath.push_back( *iImgBegin );
	    // counterexample found (could not happen if node was visited before)
	    if( !definitely( rtree.isSafe( *iImgBegin ) ) )
		return copyPath;

	    // recurse & return
	    auto posts =  rtree.postimage( *iImgBegin );
	    std::vector< typename RefinementTree< IntervalT >::NodeT > cex = findCounterexample( rtree, posts.begin(), posts.end(), copyPath );
	    if( !cex.empty() )
		return cex;
	}
    }
    return {};
}

// implement this using lower kleenean?
/*! 
  \param beginCounter and endCounter iterators to beginning and end of counterexample trajectory, should dereference to typename RefinementTree< IntervalT >::NodeT
  \param beginImage and endImage iterators to beginning and end of image of initial set obtained from refinement tree, should dereference to typename RefinementTree< IntervalT >::NodeT as well
  \return false if there definitely exists a point that is mapped to the terminal state of the counterexample, indeterminate otherwise, including if there does not possibly exist such a point 
  why upper kleenean?
  if return false, know for sure that counterexample is not spurious because a point exist with trajectory leading to unsafe state
  if return true center point did not map along trajectory
  \todo allow divergence from supposed counterexample, i.e. follow trajectory of center point until loop
*/
template< typename IntervalT, typename PathIterT, typename ImageIterT >
Ariadne::ValidatedUpperKleenean isSpurious( const RefinementTree< IntervalT >& rtree
					    , PathIterT beginCounter, PathIterT endCounter
					    , ImageIterT beginImage, ImageIterT endImage
					    , const Ariadne::Effort& effort )
{
    typedef RefinementTree< IntervalT > Rtree;
    
    std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > oBeginCex = rtree.nodeValue( *beginCounter );
    // always unsafe node is terminal state because it is unsafe
    // need this check, because first node in counterexample might be always unsafe
    if( !oBeginCex )
	return false;
    // one or less nodes in counterexample
    if( beginCounter + 1 >= endCounter )
	return false;
    // center is not contained in initial image
    if( std::none_of( beginImage, endImage,
		      [&rtree, &oBeginCex, &effort] ( const typename Rtree::NodeT& imgNode )
		      {
			  const typename Rtree::EnclosureT& startCexEnc = oBeginCex.value().get().getEnclosure();
			  std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > oImg = rtree.nodeValue( imgNode );
			  if( !oImg )
			      return false;			  
			  Ariadne::Kleenean contained = possibly( oImg.value().get().getEnclosure().contains( startCexEnc.centre() ) );
			  return possibly( contained.check( effort ) );
		      } ) )
	return true;
	
    Ariadne::Box< IntervalT > currBox = oBeginCex.value().get().getEnclosure();
    Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > currPoint = currBox.centre();
    const typename Rtree::EnclosureT& rtEnc = tree::value( rtree.tree(), tree::root( rtree.tree() ) )->getEnclosure();
    for( ; beginCounter != endCounter - 1; ++beginCounter )
    {
	std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > oNext = rtree.nodeValue( *(beginCounter + 1) );
	Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > mappedPoint = rtree.dynamics().evaluate( currPoint );
	Ariadne::ValidatedKleenean containsMapped;

	if( oNext )
	{
	    Ariadne::Box< IntervalT > nextBox = oNext.value().get().getEnclosure();
	    containsMapped = nextBox.contains( mappedPoint );
	}
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
template< typename IntervalT >
std::pair< Ariadne::ValidatedKleenean
	   , std::vector< typename RefinementTree< IntervalT >::NodeT > > cegar( RefinementTree< IntervalT >& rtree
										 , const typename RefinementTree< IntervalT >::EnclosureT& initialSet
										 , const Ariadne::Effort& effort
										 , const IRefinementStrategy< IntervalT >& refinementStrat
										 , const uint maxNodes )
{
    class NodeComparator;
    typedef RefinementTree< IntervalT > Rtree;
    typedef std::set< typename Rtree::NodeT, NodeComparator > NodeSet;

    class NodeComparator
    {
      public:
	NodeComparator( const Rtree& rtree ) : mRtree( rtree ) {}
		       
	// NodeComparator( const NodeComparator& orig ) : mRtree( orig.mRtree ) {}

	// NodeComparator& operator =( const NodeComparator& orig ) { mRtree = orig.mRtree; return *this; }

	bool operator ()( const typename Rtree::NodeT& n1
			  , const typename Rtree::NodeT& n2 ) const
	{
	    std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > otval1 = mRtree.get().nodeValue( n1 )
		, otval2 = mRtree.get().nodeValue( n2 );
	    // always unsafe node is always equal to always unsafe node
	    if( !otval1 && !otval2 )
		return true;
	    if( !otval1 || !otval2 )
		return false;
	    return otval1.value().get().id() < otval2.value().get().id();
	}

      private:
	std::reference_wrapper< const Rtree > mRtree;
    };

    // \todo test wether repopulating image anew every iteration
    // or identifying and removing refined nodes is more efficient
    NodeSet initialImage = NodeSet( NodeComparator( rtree ) );
    {
	auto img = rtree.image( initialSet );
	initialImage.insert( img.begin(), img.end() );
    }
    
    while( rtree.tree().size() < maxNodes )
    {
	// look for counterexample
	auto counterexample = findCounterexample( rtree, initialImage.begin(), initialImage.end() );
	if( counterexample.empty() )
	{
	    return std::make_pair( true, std::vector< typename Rtree::NodeT >() );
	}

	bool definitelyNotSpurious = definitely( !isSpurious( rtree
							      , counterexample.begin(), counterexample.end()
							      , initialImage.begin(), initialImage.end()
							      , effort ) )
	    , definitelyUnsafe = definitely( !rtree.isSafe( counterexample.back() ) );
	
	if( definitelyNotSpurious && definitelyUnsafe )
	{
	    return std::make_pair( false, counterexample );
	}

	for( uint i = 0; i < counterexample.size(); ++i )
	{
	    std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > oref = rtree.nodeValue( counterexample[ i ] );
	    if( oref )
	    {
		const IGraphValue& graphValRef = *graph::value( rtree.leafMapping(), counterexample[ i ] ); // dont use this after refinement
		const typename Rtree::RefinementT::NodeT& treeNodeRef = static_cast< const InsideGraphValue< typename Rtree::RefinementT::NodeT >& >( graphValRef ).treeNode();
		rtree.refine( counterexample[ i ], refinementStrat );
		// iterate over components of initial set and add their images starting from refined node for best performance
		for( const typename Rtree::NodeT& initial : initialImage )
		{
		    std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > treeValInitial = rtree.nodeValue( initial );
		    if( treeValInitial )
		    {
			auto refinedImg = rtree.image( initialSet, treeNodeRef );
			initialImage.insert( refinedImg.begin(), refinedImg.end() );
		    }
		}
		initialImage.erase( counterexample[ i ] );
	    }
	}
    }
    return make_pair( Ariadne::ValidatedKleenean( Ariadne::indeterminate ), std::vector< typename Rtree::NodeT >() );
}

#endif
