#ifndef REFINEMENT_TREE_HPP
#define REFINEMENT_TREE_HPP

#include "adjacencyDiGraph.hpp"
#include "linkedFixedBranchTree.hpp"
#include "refinementStrategy.hpp"

#include "geometry/box.hpp"
#include "geometry/function_set.hpp"
#include "numeric/logical.hpp"
#include "function/constraint.hpp"

#include <exception>
#include <set>
#include <functional>
#include <algorithm>
#include <limits>
#include <memory>

template< typename I >
class RefinementTree
{
  public:
    class InteriorTreeValue;
    struct IGraphValue;
    class InsideGraphValue;
    struct OutsideGraphValue;
    
    typedef I IntervalT;
    typedef Ariadne::Box< IntervalT > EnclosureT;
    // use shared_ptr because data structures require copy constructibility
    // refinement tree: stores pointers to values that vary depending on whether they are stored in leafs or interior nodes
    typedef tree::LinkedFixedBranchTree< std::shared_ptr< InteriorTreeValue >, 2 > RefinementT;
    // mapping graph: stores pointers to values that are either an "always-unsafe-node" or regular node storing a tree node
    typedef graph::AdjacencyDiGraph< std::shared_ptr< IGraphValue >, graph::VecMap, graph::InVec, graph::InVec > MappingT;
    typedef typename MappingT::VertexT NodeT;

    //! \class value to store in refinement tree
    class InteriorTreeValue
    {
      public:
	InteriorTreeValue()
	    : mEnclosure( Ariadne::Box< IntervalT >::zero( 0 ) )
	    , mId( -1 )                     // should give some weird value that's easy to spot
	    , mSafe( false )
	{}
	
	InteriorTreeValue( const unsigned long& id, const EnclosureT& e, const Ariadne::ValidatedLowerKleenean& safe )
	    : mId( id )
	    , mEnclosure( e )
	    , mSafe( safe )
	{}

	InteriorTreeValue( const InteriorTreeValue& orig ) = default;

	InteriorTreeValue& operator =( const InteriorTreeValue& orig ) = default;

	//! \return unique identifier of this value
	const unsigned long& id() const { return mId; }

	//! \return box stored
	const EnclosureT& getEnclosure() const { return mEnclosure; }

	//! \return true if the box stored is completely covered by all constraints, indeterminate if it lies within the initial abstraction but is not covered completely by all constraints and false if it lies outside of the initial abstraction
	Ariadne::ValidatedLowerKleenean isSafe() const { return mSafe; }

	bool operator ==( const InteriorTreeValue& tv ) const
	{
	    return this->mId == tv.mId;
	}

	virtual bool isLeaf() const { return false; }
	
      private:
	unsigned long mId;
	EnclosureT mEnclosure;
	Ariadne::ValidatedLowerKleenean mSafe;
    };

    //! \class value to store in leaves of refinement tree, does map to graph
    class LeafTreeValue : public InteriorTreeValue
    {
      public:
	LeafTreeValue() {}

	LeafTreeValue( const unsigned long& id, const EnclosureT& e, const Ariadne::ValidatedLowerKleenean& safe )
	    : InteriorTreeValue( id, e, safe )
	{}

	LeafTreeValue( const LeafTreeValue& orig ) = default;

	LeafTreeValue& operator =( const LeafTreeValue& orig ) = default;

	const typename MappingT::VertexT& graphNode() const { return mGraphNode; }
	
	bool isLeaf() const { return true; }

	void setGraphNode( const typename MappingT::VertexT& v ) { mGraphNode = v; }
      private:
	typename MappingT::VertexT mGraphNode;
    };

    //! \class abstract base class for value to store in graph
    struct IGraphValue
    {
	virtual bool isInside() const = 0;
    };

    //! \class value to store in graph of region inside first initial abstraction
    class InsideGraphValue : public IGraphValue
    {
      public:
	InsideGraphValue( typename RefinementT::NodeT& treeNode )
	    : mTreeNode( treeNode )
	{}

	InsideGraphValue( const InsideGraphValue& orig ) = default;

	InsideGraphValue& operator =( const InsideGraphValue& orig ) = default;

	const typename RefinementT::NodeT& treeNode() const { return mTreeNode; }

	bool isInside() const { return true; }
      private:
	typename RefinementT::NodeT mTreeNode;
    };

    //! \class value to store in graph of region outside first initial abstractio
    struct OutsideGraphValue : public IGraphValue
    {
	bool isInside() const { return false; }
    };

    //! \return box containing the entire state space representable by doubles
    static Ariadne::Box< IntervalT > universeBox( size_t dimension )
    {
	Ariadne::Array< IntervalT > is( dimension );
	std::generate( is.begin(), dimension, [] () { return IntervalT( std::numeric_limits< double >::min()
									, std::numeric_limits< double >::max() ); } );
	return Ariadne::Box< IntervalT >( Ariadne::Vector< IntervalT >( is ) );
    }

    //! \return perform nsplits on toSplit on dimension ndim
    //! \todo move to refiner class
    static std::vector< Ariadne::Box< IntervalT > > splitN( const Ariadne::Box< IntervalT >& toSplit, size_t ndim, size_t nsplits )
    {
	Ariadne::Array< IntervalT > intervals = toSplit.array();
	const typename IntervalT::LowerBoundType& l = toSplit.get( ndim ).lower();
	const typename IntervalT::UpperBoundType& u = toSplit.get( ndim ).upper();
	const typename IntervalT::WidthType modW = toSplit.get( ndim ).width() / nsplits;
	std::vector< Ariadne::Box< IntervalT > > splits;
	splits.reserve( nsplits );
	for( uint i = 0; i < nsplits; ++i )
	{
	    intervals[ ndim ] = IntervalT( l + i * modW, l + (i + 1) * modW );
	    splits.emplace_back( intervals );
	}
	return splits;
    }

    /*!
      \param initial abstraction using a single box
      \param constraints safety constraints to be satisfied
      \param dynamics function describing evolution of points
    */
    RefinementTree( const EnclosureT& initial
		    , const Ariadne::ConstraintSet& constraints
		    , const Ariadne::EffectiveVectorFunction dynamics
		    , const Ariadne::Effort effort
		    )
	: mConstraints( constraints )
	, mDynamics( dynamics )
	, mRefinements( typename RefinementT::ValueT( new LeafTreeValue( 0, initial, constraints.covers( initial ).check( effort ) )
						      //				   , universeBox( initial.dimension() )
						      //, constraints.covers( universeBox( initial.dimension() ) ).check( effort )
						      ) )
	, mEffort( effort )
    {
	// add always unsafe node
	auto iAddedUnsafe = graph::addVertex( mMappings, typename MappingT::ValueT( new OutsideGraphValue() ) );
	if( iAddedUnsafe == graph::vertices( mMappings ).second )
	    throw std::logic_error( "always unsafe node added but iterator to end returned" );
	mOutsideNode = *iAddedUnsafe;
	
	// set up root
	typename RefinementT::NodeT rt = root( mRefinements );
	NodeT initialNode = addToGraph( rt );
	const InteriorTreeValue& otv = nodeValue( initialNode ).value(); // should never fail
	// self loop root -> root
	if( possibly( isReachable( otv, initialNode ) ) )
	    addEdge( mMappings, initialNode, initialNode );
	// transition to unsafe: only can transition to unsafe, not from it
	if( possibly( isReachable( otv, mOutsideNode ) ) )
	    addEdge( mMappings, initialNode, mOutsideNode );
    }

    //! \return constraints determining the safe set
    const Ariadne::ConstraintSet& constraints() const
    {
	return mConstraints;
    }

    const Ariadne::EffectiveVectorFunction& dynamics() const
    {
	return mDynamics;
    }

    //! \return tree storing iterative refinements
    const RefinementT& tree() const
    {
	return mRefinements;
    }

    //! \return graph storing reachability between leaves
    const MappingT& leafMapping() const
    {
	return mMappings;
    }

    //! \return tree value stored at node v, storing the box and safety
    std::optional< std::reference_wrapper< const InteriorTreeValue > > nodeValue( const NodeT& v ) const
    {
	const typename MappingT::ValueT& gval = graph::value( mMappings, v );
	if( gval->isInside() )
	{
	    const typename RefinementT::ValueT& tval = tree::value( mRefinements, static_cast< InsideGraphValue& >( *gval ).treeNode() );
	    return std::make_optional( std::reference_wrapper< const InteriorTreeValue >( *tval ) );
	}
	else
	    return std::nullopt;
    }

    //! \return true if n is certain to be safe, false if it is certain to be unsafe, indeterminate otherwise
    Ariadne::ValidatedKleenean isSafe( const NodeT& n ) const
    {
	std::optional< std::reference_wrapper< const InteriorTreeValue > > on = nodeValue( n );
	if( on )
	{
	    if( definitely( on.value().get().isSafe() ) )
		return true;
	    else
		return Ariadne::indeterminate;
	}
	else
	    return false;
    }

    //! \return the always unsafe node used
    const NodeT& outside() const { return mOutsideNode; }

    //! \param from abstraction for which to find image in leaves of tree; needs to be of type that can be intersected with EnclosureT
    //! \return most refined boxes intersecting with from, including outside node
    template< typename EnclosureT2 >
    std::vector< NodeT > image( const EnclosureT2& from ) const
    {
	return image( from, root( mRefinements ) );
    }

    //! \param from abstraction for which to find image in leaves of tree; needs to be of type that can be intersected with EnclosureT
    //! \param subtreeRoot require leaves to be rooted at this node
    //! \return most refined boxes intersecting with from, including outside node
    template< typename EnclosureT2 >
    std::vector< NodeT > image( const EnclosureT2& from, const typename RefinementT::NodeT& subtreeRoot ) const
    {
	std::vector< NodeT > img = imageRecursive( from, subtreeRoot );
	// code non functional as of now, requires tree value to store ValidatedKleenean as safety flag
	auto iUnsafe = std::find_if( img.begin(), img.end(), [this] (const NodeT& n) { return definitely( !isSafe( n ) ); } );
	if( iUnsafe != img.end() )
	    img.push_back( outside() );
	return img;
    }

    // all node accessors: can be declared const, because tree or graph need to be accessed in order to change anything
    //! \return all leaves of the tree
    std::vector< NodeT > leaves() const 
    {
	return leaves( tree::root( mRefinements ) );
    }
    
    //! \return all leaves in subtree at v
    std::vector< NodeT > leaves( const NodeT& v ) const
    {
	const typename MappingT::ValueT& gval = graph::value( mMappings, v );
	if( gval->isInside() )
	    return leaves( static_cast< const InsideGraphValue& >( *gval ).treeNode() );
	else
	    return {};
    }
    
    //! \return all leaves in subtree at subRoot
    std::vector< NodeT > leaves( const typename RefinementT::NodeT& tn ) const
    {
    	std::vector< NodeT > ls;
	const InteriorTreeValue& tv = *tree::value( mRefinements, tn );
	
    	if( tree::isLeaf( mRefinements, tn ) )
    	{
	    //! \todo use static_cast later
	    ls.push_back( static_cast< const LeafTreeValue& >( tv ).graphNode() );
    	}
    	else
    	{
    	    typename tree::FixedBranchTreeTraits< RefinementT >::CRangeT cs = tree::children( mRefinements, tn );
    	    for( ; cs.first != cs.second; ++cs.first )
    	    {
    		std::vector< NodeT > rls = leaves( *cs.first );
    		ls.insert( ls.end(), rls.begin(), rls.end() );
    	    }
    	}
    	return ls;
    }

    //! \return all leaves in refinement tree mapping to from
    std::vector< NodeT > preimage( const NodeT& from ) const
    {
	std::vector< NodeT > preimg;
	typename graph::DiGraphTraits< MappingT >::InRangeT ins = graph::inEdges( mMappings, from );
	for( ; ins.first != ins.second; ++ins.first )
	    preimg.push_back( graph::source( mMappings, *ins.first ) );

	return preimg;
    }

    //! \return all leaves in refinement tree from maps to
    std::vector< NodeT > postimage( const NodeT& from ) const 
    {
	std::vector< NodeT > postimg;
	typename graph::DiGraphTraits< MappingT >::OutRangeT outs = graph::outEdges( mMappings, from );
	for( ; outs.first != outs.second; ++outs.first )
	    postimg.push_back( graph::target( mMappings, *outs.first ) );
	return postimg;
    }

    //! \return true if trg is deemed reachable from srcVal
    //! \note if always unsafe node is passed as second argument, it is reachable iff srcVal maps outside the initial abstraction
    // \todo eventually parametrize this
    // \todo does this always return a validated Kleenean? test with effective boxes at some point
    Ariadne::ValidatedUpperKleenean isReachable( const InteriorTreeValue& srcVal, const NodeT& trg ) const
    {
	const EnclosureT& srcBox = srcVal.getEnclosure();
	return isReachable( srcVal.getEnclosure(), trg );
    }
    
    //! \return true if there exists some point in src s.t. there exists a point in trg that can be reached
    Ariadne::ValidatedUpperKleenean isReachable( const EnclosureT& src, const NodeT& trg ) const
    {
	Ariadne::UpperBoxType ubMapped =  Ariadne::image( src, mDynamics );
	std::optional< std::reference_wrapper< const InteriorTreeValue > > trgVal = nodeValue( trg );
	if( trgVal )
	{
	    auto mapIntersection = Ariadne::intersection( src, trgVal.value().get().getEnclosure() );
	    Ariadne::ValidatedUpperKleenean doesInter = !mapIntersection.is_empty();
	    return doesInter;
	}
	else
	{
	    const EnclosureT& initialRefEnc = tree::value( tree(), tree::root( tree() ) )->getEnclosure();
	    Ariadne::ValidatedLowerKleenean intersectionEqual = Ariadne::intersection( initialRefEnc, ubMapped ) == src;
	    return !intersectionEqual;
	}
    }

    /*! 
      \param v leaf node in refinement tree
      \brief refines node v using r and updates
    */
    void refine( NodeT& v, const IRefinementStrategy< IntervalT >& r )
    {
	const typename MappingT::ValueT& gval = graph::value( mMappings, v );
	// interior node cannot be refined -> do nothing
	if( !gval->isInside() )
	    return;
	
	InsideGraphValue& inGval = static_cast< InsideGraphValue& >( *gval );
    	typename RefinementT::NodeT treev = inGval.treeNode();
	// interior node is already refined - THIS TEST IS WRONG! leaves are to be refined!
	// if( tree::isLeaf( mRefinements, treev ) )
	//     return;
	
    	// if v is leaf:
    	// obtain refinements of node v as list of appropriate length
    	EnclosureT obox = tree::value( mRefinements, treev )->getEnclosure();
	// make new tree values
    	std::vector< EnclosureT > refined = r.refine( obox );
    	std::array< std::shared_ptr< InteriorTreeValue >, RefinementT::N > tvals;
	for( uint i = 0; i < refined.size(); ++i )
	{
	    const EnclosureT& refd = refined[ i ];
	    tvals[ i ] = typename RefinementT::ValueT( new LeafTreeValue( tree().size() + i
									  , refined[ i ]
									  , constraints().covers( refined[ i ] ).check( mEffort ) ) );
	}
	tree::expand( mRefinements, treev, tvals );

    	// add expansions to the graph
    	std::pair< typename RefinementT::CIterT, typename RefinementT::CIterT > cs = tree::children( mRefinements, treev );
    	std::vector< NodeT > refinedNodes;
    	for( ; cs.first != cs.second; ++cs.first )
	{
	    typename RefinementT::NodeT refd = *cs.first;
	    refinedNodes.push_back( addToGraph( refd ) );
	}
	
    	// add edges
	//! \todo see whether it's worthwhile to use InsideGraphValues directly instead of relying on pre- and postimage to do the job
	std::vector< NodeT > pres( preimage( v ) );
	std::vector< NodeT > posts( postimage( v ) );
	// append all children of refined node, because center mapping refinements may introduce new links absent at parent level
	// simply adding v is fine, because pres and posts will be traced down to leaf
	pres.push_back( v );
	posts.push_back( v );
	
    	for( NodeT& refined : refinedNodes )
    	{
    	    // determine which elements of the preimage of v map to which refined component
    	    for( NodeT& pre : pres )
    	    {
    		for( NodeT& preLeaf : leaves( pre ) )
    		{
		    const InteriorTreeValue& preVal = nodeValue( preLeaf ).value(); // obtained as leaf, so should have value
    		    if( possibly( isReachable( preVal, refined ) ) )
			graph::addEdge( mMappings, preLeaf, refined );
    		}
    	    }

    	    // determine which elements of the postimage the components of v map to
	    const InteriorTreeValue& refinedVal = nodeValue( refined ).value(); // refined node is also in tree
    	    for( NodeT& post : posts )
    	    {
		// outside value has no leaves
		if( !graph::value( mMappings, post )->isInside() && possibly( isReachable( refinedVal, post ) ) )
		    graph::addEdge( mMappings, refined, post );
		
    		for( NodeT& postLeaf : leaves( post ) )
    		{
		    if( possibly( isReachable( refinedVal, postLeaf ) ) )
			graph::addEdge( mMappings, refined, postLeaf );
    		}
    	    }
    	}

	// <--- inefficient alternative below
	// to bring this back to life: map whole box and test image as intersection

	// auto lvs = leaves( tree::root( mRefinements ) );
	// for( NodeT& refined : refinedNodes )
	// {
	//     for( NodeT& lf : lvs )
	//     {
	// 	if( possibly( isReachable( lf, refined ) ) )
	// 	    graph::addEdge( mMappings, lf, refined );
	// 	if( possibly( isReachable( refined, lf ) ) )
	// 	    graph::addEdge( mMappings, refined, lf );
	//     }
	// }
	
	// unlink v from the graph
	removeFromGraph( v );
    }

  private:

    // helper function for recursive calls of image
    // abstract same code as leaves in DFS controller
    /*!
      \param from enclosure whose image to find
      \param to image boxes located in subtree rooted at to
      \return all nodes intersecting with from, but EXCLUDING always unsafe
    */
    template< typename EnclosureT2 >
    std::vector< NodeT > imageRecursive( const EnclosureT2& from, const typename RefinementT::NodeT& to ) const
    {
    	// vector because there can't be duplicates
    	std::vector< NodeT > parts;
	const InteriorTreeValue& tv = *tree::value( mRefinements, to );
    	const EnclosureT& boxTo = tv.getEnclosure();
    	auto inter = Ariadne::intersection( from, boxTo );
	// result is either Boolean for ExactInterval or LowerKleenean otherwise -> convert implicitly to latter if necessary
	//! \todo would like to use LowerKleenean only and check with specific effort
	Ariadne::ValidatedLowerKleenean isEmpty = inter.is_empty();

	// if there's no chance that to and from intersect: return empty
	if( definitely( isEmpty/*.check( mEffort )*/ ) )
    	    return parts;
    	else if( tree::isLeaf( mRefinements, to ) )
    	{
	    parts.push_back( static_cast< const LeafTreeValue& >( tv ).graphNode() );
    	}
    	else
    	{
    	    typename tree::FixedBranchTreeTraits< RefinementT >::CRangeT cs = tree::children( mRefinements, to );
    	    for( ; cs.first != cs.second; ++cs.first )
    	    {
    		std::vector< NodeT > rns = imageRecursive( from, *cs.first );
    		parts.insert( parts.end(), rns.begin(), rns.end() );
    	    }
    	}
    	return parts;
    }

    //! \brief add interior vertex to graph and ensure that tn holds the vertex
    typename MappingT::VertexT addToGraph( typename RefinementT::NodeT& tn )
    {
	auto pAddedVal = std::shared_ptr< IGraphValue >( new InsideGraphValue( tn ) );
	const typename MappingT::VertexT& vadded = *mMappings.addVertex( pAddedVal );
	static_cast< LeafTreeValue& >( *tree::value( mRefinements, tn ) ).setGraphNode( vadded );
	return vadded;
    }

    //! \brief removes v from the graph and transforms the tree node stored into an interior node
    void removeFromGraph( NodeT& n )
    {
	const typename MappingT::ValueT& gval = graph::value( mMappings, n );
	if( gval->isInside() )
	{
	    InsideGraphValue& inGval = static_cast< InsideGraphValue& >( *gval );
	    typename RefinementT::ValueT& tval = tree::value( mRefinements, inGval.treeNode() );
	    tval.reset( new InteriorTreeValue( *tval ) ); // copy & slice
	}
	graph::removeVertex( mMappings, n );
    }
    
    Ariadne::ConstraintSet mConstraints;
    Ariadne::EffectiveVectorFunction mDynamics;
    Ariadne::Effort mEffort;
    RefinementT mRefinements;
    MappingT mMappings;
    NodeT mOutsideNode;
};

/*!
  SHOULDN'T I USE REFINEMENT TREE NODES RATHER THAN BARE BONES TREE NODES? why?!

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
	std::optional< std::reference_wrapper< const typename Rtree::InteriorTreeValue > > oImg = rtree.nodeValue( *iImgBegin );
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
			      std::optional< std::reference_wrapper< const typename Rtree::InteriorTreeValue > > on = rtree.nodeValue( n );
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
    
    std::optional< std::reference_wrapper< const typename Rtree::InteriorTreeValue > > oBeginCex = rtree.nodeValue( *beginCounter );
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
			  std::optional< std::reference_wrapper< const typename Rtree::InteriorTreeValue > > oImg = rtree.nodeValue( imgNode );
			  if( !oImg )
			      return false;			  
			  Ariadne::Kleenean contained = possibly( oImg.value().get().getEnclosure().contains( startCexEnc.centre() ) );
			  return possibly( contained.check( effort ) );
		      } ) )
	return true;
	
    Ariadne::Box< IntervalT > currBox = oBeginCex.value().get().getEnclosure();
    for( ; beginCounter != endCounter - 1; ++beginCounter )
    {
	std::optional< std::reference_wrapper< const typename Rtree::InteriorTreeValue > > oNext = rtree.nodeValue( *(beginCounter + 1) );
	if( !oNext )     // always unsafe node is terminal
	    return false;
	Ariadne::Box< IntervalT > nextBox = oNext.value().get().getEnclosure();
	if( definitely( !nextBox.contains( currBox.centre() ) ) )
	    return true; // should be indeterminate
	currBox = nextBox;
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
  \return sequence of nodes that forms a trajectory starting from the initial set
 */
template< typename IntervalT >
std::vector< typename RefinementTree< IntervalT >::NodeT > cegar( RefinementTree< IntervalT >& rtree
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
	    std::optional< std::reference_wrapper< const typename Rtree::InteriorTreeValue > > otval1 = mRtree.get().nodeValue( n1 )
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
	    return std::vector< typename Rtree::NodeT >();
	}

	std::cout << "found counterexample " << std::endl;
	for( auto trajNode : counterexample )
	{
	    std::optional< std::reference_wrapper< const typename Rtree::InteriorTreeValue > > otrajValue = rtree.nodeValue( trajNode );
	    if( otrajValue )
		std::cout << otrajValue.value().get().getEnclosure();
	    else
		std::cout << "[unsafe]";
	    std::cout << "  ->  ";
	}
	std::cout << std::endl;

	std::optional< std::reference_wrapper< const typename Rtree::InteriorTreeValue > > otermVal = rtree.nodeValue( counterexample.back() );
	bool definitelyNotSpurious = definitely( !isSpurious( rtree
							      , counterexample.begin(), counterexample.end()
							      , initialImage.begin(), initialImage.end()
							      , effort ) )
	    , definitelyUnsafe = false;
	if( !otermVal )
	    definitelyUnsafe = true;
	else
	    definitely( rtree.constraints().separated( otermVal.value().get().getEnclosure() ) );
	
	if( definitelyNotSpurious && definitelyUnsafe )
	{
	    std::cout << "non spurious counterexample found" << std::endl;
	    return counterexample;
	}
	else if( !definitelyNotSpurious )
	    std::cout << "it's spurious" << std::endl;
	else
	    std::cout << "terminal node is not completely unsafe" << std::endl;

	// want to refine last state as well
	std::cout << "collecting image in refined tree" << std::endl;
	for( uint i = 0; i < counterexample.size(); ++i )
	{
	    std::optional< std::reference_wrapper< const typename Rtree::InteriorTreeValue > > oref = rtree.nodeValue( counterexample[ i ] );
	    std::cout << "refining box ";
	    if( oref )
		std::cout << oref.value().get().getEnclosure();
	    else
		std::cout << "[unsafe]";
	    std::cout << std::endl;

	    if( oref )
	    {
		const typename Rtree::IGraphValue& graphValRef = *graph::value( rtree.leafMapping(), counterexample[ i ] ); // dont use this after refinement
		const typename Rtree::RefinementT::NodeT& treeNodeRef = static_cast< const typename Rtree::InsideGraphValue& >( graphValRef ).treeNode();
		rtree.refine( counterexample[ i ], refinementStrat );
		// iterate over components of initial set and add their images starting from refined node for best performance
		for( const typename Rtree::NodeT& initial : initialImage )
		{
		    std::optional< std::reference_wrapper< const typename Rtree::InteriorTreeValue > > treeValInitial = rtree.nodeValue( initial );
		    if( treeValInitial )
		    {
			// use shortcut later
			auto refinedImg = rtree.image( treeValInitial.value().get().getEnclosure(), treeNodeRef );
			// surely works
			// auto refinedImg = rtree.image( itvInitial.value().get().getEnclosure() );
			initialImage.insert( refinedImg.begin(), refinedImg.end() );
		    }
		// }
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
    return {};
}

#endif
