#ifndef REFINEMENT_TREE_HPP
#define REFINEMENT_TREE_HPP

#include "adjacencyDiGraph.hpp"
#include "linkedFixedBranchTree.hpp"
#include "refinement.hpp"
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

/*!
  \class stores iterative refinements in tree and maintains links between leaf nodes
  \param E type of enclosure stored
*/
template< typename E >
class RefinementTree
{
  public:
    typedef E EnclosureT;
    // use shared_ptr because data structures require copy constructibility
    // refinement tree: stores pointers to values that vary depending on whether they are stored in leafs or interior nodes
    typedef tree::LinkedFixedBranchTree< std::shared_ptr< InteriorTreeValue< EnclosureT >  >, 2 > RefinementT;
    // mapping graph: stores pointers to values that are either an "always-unsafe-node" or regular node storing a tree node
    typedef graph::AdjacencyDiGraph< std::shared_ptr< IGraphValue >, graph::VecMap, graph::InVec, graph::InVec > MappingT;
    typedef typename MappingT::VertexT NodeT;

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
		return false;
	    if( !otval1 )
		return false;
	    if( !otval2 )
		return true;

	    return otval1.value().get().id() < otval2.value().get().id();
	}

      private:
	std::reference_wrapper< const RefinementTree< E > > mRtree;
    };

    template< typename SpaceT >
    static std::function< Ariadne::ValidatedLowerKleenean( const EnclosureT&, const SpaceT& ) > mDummyInsidePredicate = [] (auto& enc, auto& s) {
	throw std::logic_error( "this predicate should have never been called" ); };

    template< typename SpaceT >
    static std::function< Ariadne::ValidatedUpperKleenean( const EnclosureT&, const SpaceT& ) > mDummyOutsidePredicate = [] (auto& enc, auto& s) {
	throw std::logic_error( "this predicate should have never been called" ); };

    static Ariadne::ExactBoxType upper2ExactBox( const Ariadne::UpperBoxType& ub )
    {
	Ariadne::Array< Ariadne::ExactIntervalType > ints( ub.dimension() );
	for( uint i = 0; i < ub.dimension(); ++i )
	    ints[ i ] = Ariadne::ExactIntervalType( cast_exact( ub[ i ].lower() ), cast_exact( ub[ i ].upper() ) );
	return Ariadne::ExactBoxType( ints );
    }
    
    /*!
      \param initial abstraction using a single box
      \param constraints safety constraints to be satisfied
      \param dynamics function describing evolution of points
    */
    RefinementTree( const Ariadne::BoundedConstraintSet& safeSet
		    , const Ariadne::EffectiveVectorFunction dynamics
		    , const Ariadne::Effort effort
		    )
	: mSafeSet( safeSet )
	, mDynamics( dynamics )
	, mEffort( effort )
	, mNodeIdCounter( 0 )
	, mRefinements( typename RefinementT::ValueT( makeLeaf( upper2ExactBox( safeSet.bounding_box() ) ) ) )
    {
	// add outside node
	auto iAddedOutside = graph::addVertex( mMappings, typename MappingT::ValueT( new OutsideGraphValue() ) );
	if( iAddedOutside == graph::vertices( mMappings ).second )
	    throw std::logic_error( "always unsafe node added but iterator to end returned" );
	mOutsideNode = *iAddedOutside;
	
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
    const Ariadne::BoundedConstraintSet& constraints() const
    {
	return mSafeSet;
    }

    const Ariadne::EffectiveVectorFunction& dynamics() const
    {
	return mDynamics;
    }

    const Ariadne::Effort effort() const { return mEffort; }

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
	    return Ariadne::ValidatedKleenean( false );
    }

    //! \return the always unsafe node used
    const NodeT& outside() const { return mOutsideNode; }

    const EnclosureT& initialEnclosure() const { return tree::value( tree(), tree::root( tree() ) )->getEnclosure(); }

    //! \param from abstraction for which to find image in leaves of tree; needs to be of type that can be intersected with EnclosureT
    //! \return most refined boxes intersecting with from, including outside node
    template< typename EnclosureT2 >
    std::vector< NodeT > intersection( const EnclosureT2& from ) const
    {
	return intersection( from, root( mRefinements ) );
    }

    //! \param from abstraction for which to find image in leaves of tree; needs to be of type that can be intersected with EnclosureT
    //! \param subtreeRoot require leaves to be rooted at this node
    //! \return most refined boxes intersecting with from, including outside node
    template< typename EnclosureT2 >
    std::vector< NodeT > intersection( const EnclosureT2& with, const typename RefinementT::NodeT& subtreeRoot ) const
    {
	std::function< Ariadne::ValidatedUpperKleenean( const EnclosureT&, const EnclosureT2& ) > inter =
	    [this] (const EnclosureT& n, const EnclosureT2& with ) {
	    Ariadne::ValidatedLowerKleenean emptyIntersection = Ariadne::intersection( n, with ).is_empty();
	    return !emptyIntersection; };
	
	std::vector< NodeT > abs = intersectionRecursive( subtreeRoot, with, inter );
	// image recursive does not add outside
	if( possibly( overlaps( with, outside() ) ) )
	    abs.push_back( outside() );
	return abs;
    }

    //! \brief generalization of other intersection overloads
    //! \param pred function overapproximating whether enclosure and space intersect
    template< typename ConstraintSetT >
    std::vector< NodeT > intersection( const typename RefinementT::NodeT& subRoot, const ConstraintSetT& s
				       , const std::function< Ariadne::ValidatedUpperKleenean( const EnclosureT&, const ConstraintSetT& ) >& pred )
    {
	std::vector< NodeT > abs = intersectionRecursive( subRoot, s, pred );
	// relates is applied to outside node, so dummy is never called
	if( possibly( overlapsConstraints( s, outside() ) ) )
	    abs.push_back( outside() );
	return abs;
    }

    template< typename ConstraintSetT >
    std::vector< NodeT > intersection( const ConstraintSetT& s
				       , const std::function< Ariadne::ValidatedUpperKleenean( const EnclosureT&, const ConstraintSetT& ) >& pred )
    {
	return intersection( tree::root( tree() ), s, pred );
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

    //! \return true if trg can be reached from src
    //! \note nothing can be reached from the outside node
    Ariadne::ValidatedUpperKleenean isReachable( const NodeT& src, const NodeT& trg )
    {
	auto optSrcVal = nodeValue( src );
	if( optSrcVal )
	    return isReachable( optSrcVal.value().get().getEnclosure(), trg );
	else
	    return false;
    }
    
    //! \return true if trg is deemed reachable from srcVal
    //! \note if always unsafe node is passed as second argument, it is reachable iff srcVal maps outside the initial abstraction
    // \todo eventually parametrize this
    // \todo does this always return a validated Kleenean? test with effective boxes at some point
    // \todo remove this one eventually
    Ariadne::ValidatedUpperKleenean isReachable( const InteriorTreeValue< EnclosureT >& srcVal, const NodeT& trg ) const
    {
	return isReachable( srcVal.getEnclosure(), trg );
    }
    
    //! \return true if there exists some point in src s.t. there exists a point in trg that can be reached
    //! \todo prepare for generalization of boxes
    Ariadne::ValidatedUpperKleenean isReachable( const EnclosureT& src, const NodeT& trg ) const
    {
	Ariadne::UpperBoxType ubMapped =  Ariadne::image( src, mDynamics );
	std::optional< std::reference_wrapper< const InteriorTreeValue< EnclosureT > > > trgVal = nodeValue( trg );
	if( trgVal )
	{
	    auto mapIntersection = Ariadne::intersection( ubMapped, trgVal.value().get().getEnclosure() );
	    Ariadne::ValidatedUpperKleenean doesInter = !mapIntersection.is_empty();
	    return doesInter;
	}
	else
	{
	    Ariadne::RealBox initialAbs( initialEnclosure() );
	    Ariadne::BoundedConstraintSet insideInitialAbs( initialAbs );
	    Ariadne::ExactBoxType mappedOverapprox = upper2ExactBox( ubMapped );
	    Ariadne::ValidatedLowerKleenean imageContainedInside = insideInitialAbs.covers( mappedOverapprox, mEffort );
	    // Ariadne::ValidatedUpperKleenean intersectionEqual = Ariadne::intersection( initialEnclosure(), ubMapped ) == ubMapped;
	    // return intersectionEqual;
	    return !imageContainedInside;
	}
    }

    //! \return true if n1 overlaps with n2
    Ariadne::ValidatedUpperKleenean overlaps( const NodeT& n1, const NodeT& n2 ) const
    {
    	auto v1 = nodeValue( n1 ), v2 = nodeValue( n2 );
    	if( v1 && v2 )
    	    return !Ariadne::intersection( v1.value().get().getEnclosure(), v2.value().get().getEnclosure() ).is_empty();
	else if( !v1 && !v2 )
	    return Ariadne::ValidatedUpperKleenean( true );          // outside always overlaps with itself
	else if( !v1 )
    	    return !(Ariadne::intersection( initialEnclosure(), v2.value().get().getEnclosure() ) == v2.value().get().getEnclosure() );
    	else
    	    return overlaps( n2, n1 );
    }

    template< typename E2 >
    Ariadne::ValidatedUpperKleenean overlaps( const E2& otherEnclosure, const NodeT& n ) const
    {
	auto nval = nodeValue( n );
    	if( nval )
    	    return !(Ariadne::intersection( otherEnclosure, nval.value().get().getEnclosure() ).is_empty() );
    	else
	    return !(Ariadne::intersection( initialEnclosure(), otherEnclosure ) == otherEnclosure );
    }

    //! \return true if n overlaps with constraint
    template< typename ConstraintSetT >
    Ariadne::ValidatedUpperKleenean overlapsConstraints( const ConstraintSetT& constraintSet, const NodeT& n ) const
    {
    	auto nval = nodeValue( n );
    	if( nval )
    	    return !(constraintSet.separated( nval.value().get().getEnclosure() ).check( mEffort ) );
    	else
    	{
    	    const EnclosureT initialAbs = initialEnclosure();
    	    if( std::all_of( constraintSet.constraints().begin(), constraintSet.constraints().end()                    
    	    		     , [&initialAbs] (auto& c) {
				 Ariadne::List< Ariadne::EffectiveConstraint > single = {c};
				 return definitely( Ariadne::ConstraintSet( single ).covers( initialAbs ) ); } ) )
    	    	return Ariadne::ValidatedUpperKleenean( false );
    	    return Ariadne::ValidatedUpperKleenean( true );
    	}
    }

    // //! \return true if n1 covers n2
    // Ariadne::ValidatedKleenean covers( const NodeT& n1, const NodeT& n2 )
    // {
    // 	auto v1 = nodeValue( n1 ), v2 = nodeValue( n2 );
    // 	if( v1 && v2 )
    // 	    return toKleenean( Ariadne::intersection( v1.value().get().getEnclosure(), v2.value().get().getEnclosure() ) ==
    // 			       v2.value().get().getEnclosure() );
    // 	else if( !v1 && !v2 )
    // 	    return Ariadne::ValidatedKleenean( true );
    // 	else if( !v1 )
    // 	    return toKleenean( !Ariadne::intersection( initialEnclosure(), v2.value().get().getEnclosure() ).is_empty() );
    // 	else
    // 	    return Ariadne::ValidatedKleenean( false );
    // }p

    // //! \return true if constraint covers n
    // template< typename ConstraintT >
    // Ariadne::ValidatedKleenean covers( const ConstraintT& constraint, const NodeT& n )
    // {
    // 	auto nval = nodeValue( n );
    // 	if( nval )
    // 	    return toKleenean( constraint.covers( nval.value().get().getEnclosure() ) );
    // 	Ariadne::ValidatedKleenean( false ); // good enough: cannot cover the outside fully
    // }

    // //! \return true if n1 is separated from n2
    // Ariadne::ValidatedKleenean separated( const NodeT& n1, const NodeT& n2 )
    // {
    // 	auto v1 = nodeValue( n1 ), v2 = nodeValue( n2 );
    // 	if( v1 && v2 )
    // 	    return !overlaps( n1, n2 );
    // 	else if( !v1 && !v2 )
    // 	    return false;  // outside cannot be separated
    // 	else if( !v1 )
    // 	    return toKleenean( !(Ariadne::intersection( initialEnclosure(), v2.value().get().getEnclosure() )
    // 				 == v2.value().get().getEnclosure() ) );
    // 	else
    // 	    return separated( n2, n1 ); // is associative
    // }

    // //! \return true if n is separated from constraint
    // template< typename ConstraintSetT >
    // Ariadne::ValidatedKleenean separated( const ConstraintSetT& cset, const NodeT& n )
    // {
    // 	auto nval = nodeValue( n );
    // 	if( nval )
    // 	    return toKleenean( cset.separated( nval.value().get().getEnclosure() ) );
    // }
    
    // /*!
    //   \param predIn predicate to verify when applied to n and s if n is an inside node
    //   \param predOut predicate to falisify when applied to the initial abstraction and s if n is an outside node
    //   \return true if either n is an inside node and predIn( n, s ) is definitely satisfied or if n is an outside node and predOut( initialAbstraction, s ) is definitely not satisfied; indeterminate otherwise
    // */
    // template< typename SpaceT >
    // Ariadne::ValidatedKleenean relates( const NodeT& n, const SpaceT& s
    // 					, const std::function< Ariadne::ValidatedLowerKleenean( const EnclosureT&, const SpaceT& ) >& predIn
    // 					, const std::function< Ariadne::ValidatedUpperKleenean( const EnclosureT&, const SpaceT& ) >& predOut) const
    // {
    // 	auto val = nodeValue( n );
    // 	if( val )
    // 	{
    // 	    if( definitely( predIn( val.value().get().getEnclosure(), s ) ) )
    // 		return true;
    // 	    else
    // 		return Ariadne::indeterminate;
    // 	}
    // 	else
    // 	{
    // 	    if( definitely( !predOut( initialEnclosure(), s ) ) )
    // 		return true;
    // 	    else
    // 		return Ariadne::indeterminate;
    // 	}
    // }
    
    //! \return true if nodes are equal based on their identifiers, false otherwise
    bool equal( const NodeT& n1, const NodeT& n2 ) const
    {
	auto val1 = nodeValue( n1 ), val2 = nodeValue( n2 );
	if( val1 && val2 )
	    return val1.value().get() == val2.value().get();
	else if( !val1 && !val2 )
	    return true;
	return false;
    }

    /*! 
      \param v leaf node in refinement tree
      \brief refines node v using r and updates
    */
    template< typename R >
    void refine( const NodeT& v, R& r )
    {
	const typename MappingT::ValueT& gval = graph::value( mMappings, v );
	// interior node cannot be refined -> do nothing
	if( !gval->isInside() )
	    return;
	
	InsideGraphValue< typename RefinementT::NodeT >& inGval = static_cast< InsideGraphValue< typename RefinementT::NodeT >& >( *gval );
    	typename RefinementT::NodeT treev = inGval.treeNode();
    	EnclosureT obox = tree::value( mRefinements, treev )->getEnclosure();
	// make new tree values
    	std::vector< EnclosureT > refined = r( obox );
    	std::array< std::shared_ptr< InteriorTreeValue< EnclosureT > >, RefinementT::N > tvals;
	for( uint i = 0; i < refined.size(); ++i )
	{
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

	for( NodeT& refined : refinedNodes )
	    refineEdges( v, refined );
	// unlink v from the graph after its connectivity is no longer needed
	removeFromGraph( v );
    }

  private:

    /*!
      \param from enclosure whose image to find
      \param to image boxes located in subtree rooted at to
      \return all nodes intersecting with from, but EXCLUDING always unsafe
    */
    template< typename E2 >
    std::vector< NodeT > intersectionRecursive( const typename RefinementT::NodeT& subRoot, const E2& with
						, const std::function< Ariadne::ValidatedUpperKleenean( const EnclosureT&, const E2& ) >& pred ) const
    {
    	// vector because there can't be duplicates
    	std::vector< NodeT > parts;
	const InteriorTreeValue< EnclosureT >& subRootVal = *tree::value( mRefinements, subRoot );
    	const EnclosureT& subRootBox = subRootVal.getEnclosure();

	// will always be false: pred yields validated lower kleenean
	// so negation is validated upper kleenean which cannot be affirmed
	if( definitely( !pred( subRootBox, with ) ) )
    	    return parts;
    	else if( tree::isLeaf( mRefinements, subRoot ) )
	    parts.push_back( static_cast< const LeafTreeValue< EnclosureT, typename MappingT::VertexT >& >( subRootVal ).graphNode() );
    	else
    	{
    	    typename tree::FixedBranchTreeTraits< RefinementT >::CRangeT cs = tree::children( mRefinements, subRoot );
    	    for( ; cs.first != cs.second; ++cs.first )
    	    {
    		std::vector< NodeT > rns = intersectionRecursive( *cs.first, with, pred );
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

    //! \note adapts edges of parent node after refinement
    void refineEdges( const NodeT& parent, const NodeT& child )
    {
	// add edges
	std::vector< NodeT > pres( preimage( parent ) );
	std::vector< NodeT > posts( postimage( parent ) );
	// simply adding v is fine, because pres and posts will be traced down to leaf
	// but parent will be in pre and postimage if it has a self loop, otherwise irrelevant - I think...
	// pres.push_back( v );
	// posts.push_back( v );
	
	for( NodeT& pre : pres )
	    connectAllLeaves( pre, child );

	for( NodeT& post : posts )
	    connectAllLeaves( child, post ); 
    }

    //! \brief adds links from all leaves in subtree at src to all leaves in subtree at trg that are reachable
    void connectAllLeaves( const NodeT& src, const NodeT& trg )
    {
	// auto srcVal = nodeValue( src ), trgVal = nodeValue( trg );

	std::vector< NodeT > srcLeavesOrOutside = leaves( src )
	    , trgLeavesOrOutside = leaves( trg );
	// outside has no leaves
	if( !graph::value( leafMapping(), src )->isInside() )
	    srcLeavesOrOutside.push_back( src );
	if( !graph::value( leafMapping(), trg )->isInside() )
	    trgLeavesOrOutside.push_back( trg );

	for( NodeT& srcLeaf : srcLeavesOrOutside )
	{
	    for( NodeT& trgLeaf : trgLeavesOrOutside )
	    {
		if( possibly( isReachable( srcLeaf, trgLeaf ) ) )
		    graph::addEdge( mMappings, srcLeaf, trgLeaf );
	    }
	}
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
    void removeFromGraph( const NodeT& n )
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
    
    Ariadne::BoundedConstraintSet mSafeSet;
    Ariadne::EffectiveVectorFunction mDynamics;
    Ariadne::Effort mEffort;
    unsigned long mNodeIdCounter;
    RefinementT mRefinements;
    MappingT mMappings;
    NodeT mOutsideNode;
};

#endif
