#ifndef REFINEMENT_TREE_HPP
#define REFINEMENT_TREE_HPP

#include "adjacencyDiGraph.hpp"
#include "linkedFixedBranchTree.hpp"
#include "refinementStrategy.hpp"
#include "treeValue.hpp"
#include "graphValue.hpp"

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
    typedef I IntervalT;
    typedef Ariadne::Box< IntervalT > EnclosureT;
    // use shared_ptr because data structures require copy constructibility
    // refinement tree: stores pointers to values that vary depending on whether they are stored in leafs or interior nodes
    typedef tree::LinkedFixedBranchTree< std::shared_ptr< InteriorTreeValue< EnclosureT >  >, 2 > RefinementT;
    // mapping graph: stores pointers to values that are either an "always-unsafe-node" or regular node storing a tree node
    typedef graph::AdjacencyDiGraph< std::shared_ptr< IGraphValue >, graph::VecMap, graph::InVec, graph::InVec > MappingT;
    typedef typename MappingT::VertexT NodeT;

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
	, mEffort( effort )
	, mNodeIdCounter( 0 )
	, mRefinements( typename RefinementT::ValueT( makeLeaf( initial )
						      //				   , universeBox( initial.dimension() )
						      //, constraints.covers( universeBox( initial.dimension() ) ).check( effort )
						      ) )
    {
	// add always unsafe node
	auto iAddedUnsafe = graph::addVertex( mMappings, typename MappingT::ValueT( new OutsideGraphValue() ) );
	if( iAddedUnsafe == graph::vertices( mMappings ).second )
	    throw std::logic_error( "always unsafe node added but iterator to end returned" );
	mOutsideNode = *iAddedUnsafe;
	
	// set up root
	typename RefinementT::NodeT rt = root( mRefinements );
	NodeT initialNode = addToGraph( rt );
	const InteriorTreeValue< EnclosureT >& otv = nodeValue( initialNode ).value(); // should never fail
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
    std::optional< std::reference_wrapper< const InteriorTreeValue< EnclosureT > > > nodeValue( const NodeT& v ) const
    {
	const typename MappingT::ValueT& gval = graph::value( mMappings, v );
	if( gval->isInside() )
	{
	    const typename RefinementT::ValueT& tval = tree::value( mRefinements, static_cast< InsideGraphValue< typename RefinementT::NodeT >& >( *gval ).treeNode() );
	    return std::make_optional( std::reference_wrapper< const InteriorTreeValue< EnclosureT > >( *tval ) );
	}
	else
	    return std::nullopt;
    }

    //! \return true if n is certain to be safe, false if it is certain to be unsafe, indeterminate otherwise
    Ariadne::ValidatedKleenean isSafe( const NodeT& n ) const
    {
	std::optional< std::reference_wrapper< const InteriorTreeValue< EnclosureT > > > on = nodeValue( n );
	if( on )
	    return on.value().get().isSafe();
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
	// image recursive does not add outside
	if( definitely( !(Ariadne::intersection( from, tree::value( tree(), tree::root( tree() ) )->getEnclosure() ) == from ) ) )
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
	    return leaves( static_cast< const InsideGraphValue< typename RefinementT::NodeT >& >( *gval ).treeNode() );
	else
	    return {};
    }
    
    //! \return all leaves in subtree at subRoot
    std::vector< NodeT > leaves( const typename RefinementT::NodeT& tn ) const
    {
    	std::vector< NodeT > ls;
	const InteriorTreeValue< EnclosureT >& tv = *tree::value( mRefinements, tn );
	
    	if( tree::isLeaf( mRefinements, tn ) )
    	{
	    //! \todo use static_cast later
	    ls.push_back( static_cast< const LeafTreeValue< EnclosureT, typename MappingT::VertexT >& >( tv ).graphNode() );
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
    Ariadne::ValidatedUpperKleenean isReachable( const InteriorTreeValue< EnclosureT >& srcVal, const NodeT& trg ) const
    {
	const EnclosureT& srcBox = srcVal.getEnclosure();
	return isReachable( srcVal.getEnclosure(), trg );
    }
    
    //! \return true if there exists some point in src s.t. there exists a point in trg that can be reached
    Ariadne::ValidatedUpperKleenean isReachable( const EnclosureT& src, const NodeT& trg ) const
    {
	Ariadne::UpperBoxType ubMapped =  Ariadne::image( src, mDynamics );
	std::optional< std::reference_wrapper< const InteriorTreeValue< EnclosureT > > > trgVal = nodeValue( trg );
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
	
	InsideGraphValue< typename RefinementT::NodeT >& inGval = static_cast< InsideGraphValue< typename RefinementT::NodeT >& >( *gval );
    	typename RefinementT::NodeT treev = inGval.treeNode();
	// interior node is already refined - THIS TEST IS WRONG! leaves are to be refined!
	// if( tree::isLeaf( mRefinements, treev ) )
	//     return;
	
    	// if v is leaf:
    	// obtain refinements of node v as list of appropriate length
    	EnclosureT obox = tree::value( mRefinements, treev )->getEnclosure();
	// make new tree values
    	std::vector< EnclosureT > refined = r.refine( obox );
    	std::array< std::shared_ptr< InteriorTreeValue< EnclosureT > >, RefinementT::N > tvals;
	for( uint i = 0; i < refined.size(); ++i )
	{
	    const EnclosureT& refd = refined[ i ];
	    tvals[ i ] = typename RefinementT::ValueT( makeLeaf( refined[ i ] ) );
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
	//! \todo see whether it's worthwhile to use InsideGraphValue< typename RefinementT::NodeT >s directly instead of relying on pre- and postimage to do the job
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
		    const InteriorTreeValue< EnclosureT >& preVal = nodeValue( preLeaf ).value(); // obtained as leaf, so should have value
    		    if( possibly( isReachable( preVal, refined ) ) )
			graph::addEdge( mMappings, preLeaf, refined );
    		}
    	    }

    	    // determine which elements of the postimage the components of v map to
	    const InteriorTreeValue< EnclosureT >& refinedVal = nodeValue( refined ).value(); // refined node is also in tree
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
	const InteriorTreeValue< EnclosureT >& tv = *tree::value( mRefinements, to );
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
	    parts.push_back( static_cast< const LeafTreeValue< EnclosureT, typename MappingT::VertexT >& >( tv ).graphNode() );
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

    //! \return pointer to newly allocated leaf value ensuring that the safety flag is correctly initialized
    LeafTreeValue< EnclosureT, typename MappingT::VertexT >* makeLeaf( const EnclosureT& enc )
    {
	return new LeafTreeValue< EnclosureT, typename MappingT::VertexT >( mNodeIdCounter++, enc
				  , definitely( constraints().covers( enc ).check( mEffort ) )
				  ? Ariadne::ValidatedKleenean( true )
				  : (definitely( constraints().separated( enc ).check( mEffort ) )
				     ? Ariadne::ValidatedKleenean( false )
				     : Ariadne::indeterminate) );
    }
    
    //! \brief add interior vertex to graph and ensure that tn holds the vertex
    typename MappingT::VertexT addToGraph( typename RefinementT::NodeT& tn )
    {
	auto pAddedVal = std::shared_ptr< IGraphValue >( new InsideGraphValue< typename RefinementT::NodeT >( tn ) );
	const typename MappingT::VertexT& vadded = *mMappings.addVertex( pAddedVal );
	static_cast< LeafTreeValue< EnclosureT, typename MappingT::VertexT >& >( *tree::value( mRefinements, tn ) ).setGraphNode( vadded );
	return vadded;
    }

    //! \brief removes v from the graph and transforms the tree node stored into an interior node
    void removeFromGraph( NodeT& n )
    {
	const typename MappingT::ValueT& gval = graph::value( mMappings, n );
	if( gval->isInside() )
	{
	    InsideGraphValue< typename RefinementT::NodeT >& inGval = static_cast< InsideGraphValue< typename RefinementT::NodeT >& >( *gval );
	    typename RefinementT::ValueT& tval = tree::value( mRefinements, inGval.treeNode() );
	    tval.reset( new InteriorTreeValue< EnclosureT >( *tval ) ); // copy & slice
	}
	graph::removeVertex( mMappings, n );
    }
    
    Ariadne::ConstraintSet mConstraints;
    Ariadne::EffectiveVectorFunction mDynamics;
    Ariadne::Effort mEffort;
    unsigned long mNodeIdCounter;
    RefinementT mRefinements;
    MappingT mMappings;
    NodeT mOutsideNode;
};

#endif
