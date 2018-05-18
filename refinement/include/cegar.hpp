#ifndef CEGAR_HPP
#define CEGAR_HPP

#include "refinementTree.hpp"

/*!
  runs DFS to find counterexample
  any path terminates in
  1) loop leading back to state along path
  2) state with violated safety conditions
  \param iImgBegin iterator to beginning of refinement tree nodes describing the image of the initial set, should dereference to RefinementTree< IntervalT >::NodeT
  \return vector of nodes terminated by a possibly unsafe node
  \todo add parameter to control ordering of branches in dfs exploration 
  \todo speed up: use references in vector storing path passed to recursive calls -> don't need to copy nodes
  \todo use const refinement tree once constness issues in tree+graph are fixed
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
	// std::cout << "looking for counterex at " << rtree.nodeValue( *iImgBegin ).getEnclosure()
	// << " which is " << rtree.nodeValue( *iImgBegin ).isSafe() << " safe " << std::endl;

	// counterexample found
	// image should never contain always-unsafe-node --> that assumption is wrong!
	std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > oImg = rtree.nodeValue( *iImgBegin );
	if( !definitely( rtree.isSafe( *iImgBegin ) ) )
	{
	    std::vector< typename RefinementTree< IntervalT >::NodeT > copyPath( path.begin(), path.end() );
	    copyPath.push_back( *iImgBegin );
	    // std::cout << "counterexample of length " << copyPath.size() << std::endl;
	    return copyPath;
	}

	// look for loop										       
	typename std::vector< typename RefinementTree< IntervalT >::NodeT >::const_iterator iBeginLoop =
	    std::find_if( path.begin(), path.end()
			  , [&rtree, &iImgBegin] (const typename RefinementTree< IntervalT >::NodeT& n) {
			      std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > on = rtree.nodeValue( n );
			      if( !on )          // image is valid, otherwise would have returned 
				  return false;
			      return on.value().get() == rtree.nodeValue( *iImgBegin ).value().get(); } );
	// no loop found
	if( iBeginLoop == path.end() )
	{
	    // std::cout << "no loop, recurse!" << std::endl;
	    // recurse & return
	    std::vector< typename RefinementTree< IntervalT >::NodeT > copyPath( path.begin(), path.end() );
	    copyPath.push_back( *iImgBegin );
	    auto posts =  rtree.postimage( *iImgBegin );
	    std::vector< typename RefinementTree< IntervalT >::NodeT > cex = findCounterexample( rtree, posts.begin(), posts.end(), copyPath );
	    if( !cex.empty() )
		return cex;
	}
    }
    // no counterexample found
    // std::cout << "nothing found, go home now" << std::endl;
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
	// if( !oNext )     // always unsafe node is terminal, yes, but I'm here to check that it's reachable
	// return false;
	Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > mappedPoint = rtree.dynamics().evaluate( currPoint );
	Ariadne::ValidatedKleenean containsMapped;

	// std::cout << currPoint << "  to  " << mappedPoint << std::endl;

	if( oNext )
	{
	    Ariadne::Box< IntervalT > nextBox = oNext.value().get().getEnclosure();
	    containsMapped = nextBox.contains( mappedPoint );
	}
	else
	    containsMapped = !rtEnc.contains( mappedPoint );
	// this is nonsense! should map centre of current box
	if( definitely( !containsMapped ) )
	{
	    // std::cout << mappedPoint << " is inside " << oNext << "?" << std::endl;
	    return true; // should be indeterminate
	}
	
	currPoint = mappedPoint;
    }
    return false;
    
    
    // // intersect image with first element of counterex trajectory -> new image -> recurse
    // std::vector< typename RefinementTree< IntervalT >::NodeT > mappedImage;
    // for( ; beginImage != endImage; ++beginImage )
    // {
    // 	auto imagePosts = rtree.postimage( *beginImage );
    // 	Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > mappedCenter = rtree.dynamics().evaluate( rtree.nodeValue( *beginImage ).getEnclosure().centre() );
    // 	auto ipost = imagePosts.begin();
    // 	while( ipost != imagePosts.end() &&
    // 	       !possibly( rtree.nodeValue( *ipost ).getEnclosure().contains( mappedCenter ) ) )
    // 	    ++ipost;
	
    // 	// replaced in while loop above
    // 	// Ariadne::ValidatedLowerKleenean doesMap = trgBox.contains( mappedCentre );

    // 	// certainly cannot trace any given set -> certainly spurious
    // 	if( ipost == imagePosts.end() )
    // 	    return true;

    // 	Ariadne::ValidatedUpperKleenean continueSpurious = true; // should be indeterminate
    // 	if( beginCounter != endCounter )
    // 	    continueSpurious = isSpurious( rtree, beginCounter + 1, endCounter, ipost, ipost + 1, effort ); // ipost is not end, so can do +1

    // 	DO NOT HANDLE OTHER CASE: CURRENT BOX IS COMPLETELY UNSAFE -> this should be the last one

    // 	// can surely trace some set -> certainly not spurious
    // 	// READ THIS! \todo should add condition: unsafe node needs to be separate from safe set
    // 	if( definitely( !continueSpurious ) )
    // 	    return false;
    // }
    // // neither found instance to prove nor could disprove that the counterexample was spurious
    // return true; // should be indeterminate
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
	std::cout << "new iteration, number of nodes " << rtree.tree().size() << "/" << maxNodes << std::endl;

	// look for counterexample
	auto counterexample = findCounterexample( rtree, initialImage.begin(), initialImage.end() );
	if( counterexample.empty() )
	{
	    std::cout << "no counterexample found" << std::endl;
	    return std::make_pair( true, std::vector< typename Rtree::NodeT >() );
	}

	std::cout << "found counterexample " << std::endl;
	for( auto trajNode : counterexample )
	{
	    std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > otrajValue = rtree.nodeValue( trajNode );
	    if( otrajValue )
		std::cout << otrajValue.value().get().getEnclosure();
	    else
		std::cout << "[unsafe]";
	    std::cout << "  ->  ";
	}
	std::cout << std::endl;

	bool definitelyNotSpurious = definitely( !isSpurious( rtree
							      , counterexample.begin(), counterexample.end()
							      , initialImage.begin(), initialImage.end()
							      , effort ) )
	    , definitelyUnsafe = definitely( !rtree.isSafe( counterexample.back() ) );
	
	if( definitelyNotSpurious && definitelyUnsafe )
	{
	    std::cout << "non spurious counterexample found" << std::endl;
	    return std::make_pair( false, counterexample );
	}
	else if( !definitelyNotSpurious )
	    std::cout << "it's spurious" << std::endl;
	else
	    std::cout << "terminal node is not completely unsafe" << std::endl;

	// want to refine last state as well
	std::cout << "collecting image in refined tree" << std::endl;
	for( uint i = 0; i < counterexample.size(); ++i )
	{
	    std::optional< std::reference_wrapper< const InteriorTreeValue< typename Rtree::EnclosureT > > > oref = rtree.nodeValue( counterexample[ i ] );
	    std::cout << "refining box ";
	    if( oref )
		std::cout << oref.value().get().getEnclosure();
	    else
		std::cout << "[unsafe]";
	    std::cout << std::endl;

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
			// use shortcut later
			auto refinedImg = rtree.image( initialSet, treeNodeRef );
			// surely works
			// auto refinedImg = rtree.image( itvInitial.value().get().getEnclosure() );
			initialImage.insert( refinedImg.begin(), refinedImg.end() );
		    }
		}
		initialImage.erase( counterexample[ i ] );
	    }
	}

	// \todo test against this code
	// replaced by code in for loop over counterexs above
	// NodeSet newInitialImage = NodeSet( NodeComparator( rtree ) );
	// for( const typename Rtree::NodeT& prevImage : initialImage )
	// {
	//     auto imageNow = rtree.image( rtree.nodeValue( prevImage ).getEnclosure() );
	    
	//     newInitialImage.insert( imageNow.begin(), imageNow.end() );
	// }
	// initialImage = std::move( newInitialImage );

	std::cout << "cegar iteration done " << std::endl;
    }
    return make_pair( Ariadne::ValidatedKleenean( Ariadne::indeterminate ), std::vector< typename Rtree::NodeT >() );
}

#endif
