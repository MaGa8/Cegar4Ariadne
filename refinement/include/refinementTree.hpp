#ifndef REFINEMENT_TREE_HPP
#define REFINEMENT_TREE_HPP

#include "adjacencyDiGraph.hpp"
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
	    auto otval1 = mRtree.get().nodeValue( n1 )
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
	, mInitialEnclosure( upper2ExactBox( safeSet.bounding_box() ) )
	  // \todo rename make leaf to addState
    {
	// set up root
	NodeT initialNode = addState( mInitialEnclosure );

	// add outside node
	auto iAddedOutside = graph::addVertex( mMapping, typename MappingT::ValueT( new OutsideGraphValue() ) );
	if( iAddedOutside == graph::vertices( mMapping ).second )
	    throw std::logic_error( "always unsafe node added but iterator to end returned" );
	mOutsideNode = *iAddedOutside;

	// self loop root -> root
	if( possibly( isReachable( initialNode, initialNode ) ) )
	    addEdge( mMapping, initialNode, initialNode );
	// transition to unsafe: only can transition to unsafe, not from it
	if( possibly( isReachable( initialNode, mOutsideNode ) ) )
	    addEdge( mMapping, initialNode, mOutsideNode );
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

    //! \return graph storing reachability between leaves
    const MappingT& graph() const
    {
	return mMapping;
    }

    //! \return tree value stored at node v, storing the box and safety
    std::optional< std::reference_wrapper< const InsideGraphValue< EnclosureT > > > nodeValue( const NodeT& v ) const
    {
	const typename MappingT::ValueT& gval = graph::value( mMapping, v );
	if( gval->isInside() )
	    return std::make_optional( std::ref( static_cast< const InsideGraphValue< E >& >( *gval ) ) );
	else
	    return std::nullopt;
    }

    //! \return true if n is certain to be safe, false if it is certain to be unsafe, indeterminate otherwise
    Ariadne::ValidatedKleenean isSafe( const NodeT& n ) const
    {
	auto nval = nodeValue( n );
	if( nval )
	    return nval.value().get().isSafe();
	else
	    return Ariadne::ValidatedKleenean( false );
    }

    //! \return the always unsafe node used
    const NodeT& outside() const { return mOutsideNode; }

    const EnclosureT& initialEnclosure() const { return mInitialEnclosure; }

    //! \param from abstraction for which to find image in leaves of tree; needs to be of type that can be intersected with EnclosureT
    //! \return most refined boxes intersecting with from, including outside node
    template< typename EnclosureT2 >
    std::vector< NodeT > intersection( const EnclosureT2& from ) const
    {
	std::vector< NodeT > inters;
	auto vrange = graph::vertices( graph() );
	std::copy_if( vrange.first, vrange.second, std::back_inserter( inters ),
		      [this, &from] (auto& gn) {
			  auto gnval = nodeValue( gn );
			  if( gnval )
			      return possibly( !Ariadne::intersection( from, gnval.value().get().getEnclosure() ).is_empty() );
			  return possibly( !(Ariadne::intersection( from, initialEnclosure() ) == from) );
		      } );

	return inters;
    }

    //! \brief generalization of other intersection overloads
    //! \param pred function overapproximating whether enclosure and space intersect
    std::vector< NodeT > intersection( const Ariadne::BoundedConstraintSet& s
				       , const std::function< Ariadne::ValidatedUpperKleenean( const EnclosureT&, const Ariadne::BoundedConstraintSet& ) >& pred ) const
    {
	double largeM = 100 * initialEnclosure().radius().get_d();
	
	std::vector< NodeT > inters;
	auto vrange = graph::vertices( graph() );
	std::copy_if( vrange.first, vrange.second, std::back_inserter( inters )
		      , [this, &s, &largeM] (auto& gn) {
			    return possibly( overlapsConstraints( s, gn ) ); } );
	return inters;
    }

    //! \return all leaves in refinement tree mapping to from
    std::vector< NodeT > preimage( const NodeT& to ) const
    {
	std::vector< NodeT > preimg;
	typename graph::DiGraphTraits< MappingT >::InRangeT ins = graph::inEdges( mMapping, to );

	for( ; ins.first != ins.second; ++ins.first )
	{
	    const typename MappingT::EdgeT edge = *ins.first;
	    auto edgeSource = graph::source( mMapping, edge );
	    preimg.push_back( graph::source( mMapping, *ins.first ) );
	}
	return preimg;
    }

    //! \return all leaves in refinement tree from maps to
    std::vector< NodeT > postimage( const NodeT& from ) const 
    {
	std::vector< NodeT > postimg;
	typename graph::DiGraphTraits< MappingT >::OutRangeT outs = graph::outEdges( mMapping, from );
	for( ; outs.first != outs.second; ++outs.first )
	    postimg.push_back( graph::target( mMapping, *outs.first ) );
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
    
    //! \return true if there exists some point in src s.t. there exists a point in trg that can be reached
    //! \todo prepare for generalization of boxes
    Ariadne::ValidatedUpperKleenean isReachable( const EnclosureT& src, const NodeT& trg ) const
    {
	Ariadne::UpperBoxType ubMapped =  Ariadne::image( src, mDynamics );
	std::optional< std::reference_wrapper< const InsideGraphValue< EnclosureT > > > trgVal = nodeValue( trg );
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

    //! \note: don't use this method! breaks in case n == outsideNode()
    //! \return true if n overlaps with constraint
    Ariadne::ValidatedUpperKleenean overlapsConstraints( const Ariadne::BoundedConstraintSet& constraintSet, const NodeT& n ) const
    {
    	auto nval = nodeValue( n );
    	if( nval )
    	    return !(constraintSet.separated( nval.value().get().getEnclosure() ).check( mEffort ) );
    	else
    	{
	    // safe set determines initial enclosure
	    return !( constraints().covers( upper2ExactBox( constraintSet.bounding_box() ) ).check( effort() ) );
	    
	    // attempted generalized solution
	    // EnclosureT initial = initialEnclosure();
	    // size_t ndim = initial.dimension();
	    // // construct complement( initialEnclosure ) set
	    // Ariadne::List< Ariadne::EffectiveConstraint > complementConstraints;
	    // for( uint i = 0; i < ndim; ++i )
	    // {
	    // 	Ariadne::EffectiveScalarFunction xi = Ariadne::EffectiveScalarFunction::coordinate( Ariadne::EuclideanDomain( ndim ), i );
	    // 	Ariadne::Real l( initial[ i ].lower().get_d() )
	    // 	    , u( initial[ i ].upper().get_d() );          // is there a better way?

	    // 	complementConstraints.append( xi >= l );
	    // 	complementConstraints.append( xi <= u );
	    // }
	    // Ariadne::RealBox pseudoUniverse = { {-m, m}, {-m, m} };
	    // Ariadne::BoundedConstraintSet initialComplement( pseudoUniverse, complementConstraints )
	    // 	, universe( pseudoUniverse, Ariadne::List< Ariadne::EffectiveConstraint >() );
	    // if constraintSet intersection initialComplement overlaps, then constraintSet overlaps outside
	    // return Ariadne::overlap( Ariadne::intersection( initialComplement, constraintSet ), universe );
    	}
    }
    
    //! \return true nodes are equal based on their identifiers, false otherwise
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
      \return vector of refined nodes
    */
    template< typename R >
    std::vector< NodeT > refine( NodeT v, R& r ) // may not reference v, because modification of mMapping may cause reallocation of vertex container
    {
	auto vval = nodeValue( v );
	std::vector< NodeT > refinedStates;
	if( !vval )
	    return refinedStates;

	// make new tree values
    	std::vector< EnclosureT > refinedEnclosures = r( vval.value().get().getEnclosure() );
	refinedStates.reserve( refinedEnclosures.size() );
	// map to outside node directly
	for( auto& refEnc : refinedEnclosures )
	    refinedStates.push_back( addState( refEnc ) );

	refineEdges( v, refinedStates.begin(), refinedStates.end() );
    	
	// unlink v from the graph after its connectivity is no longer needed
	graph::removeVertex( mMapping, v );
	return refinedStates;
    }

  private:

    //! \return pointer to newly allocated leaf value ensuring that the safety flag is correctly initialized
    const NodeT& addState( const EnclosureT& enc )
    {
	Ariadne::ValidatedKleenean safety = definitely( constraints().covers( enc ).check( mEffort ) )
	    ? Ariadne::ValidatedKleenean( true )
	    : (definitely( constraints().separated( enc ).check( mEffort ) )
	       ? Ariadne::ValidatedKleenean( false )
	       : Ariadne::indeterminate);
	auto pvalue = std::shared_ptr< IGraphValue >( new InsideGraphValue< EnclosureT >( mNodeIdCounter++, enc, safety  ) );
	auto iadded = graph::addVertex( mMapping, pvalue );
	return *iadded;
    }

    //! \note adapts edges of parent node after refinement
    template< typename IterT >
    void refineEdges( const NodeT& parent, const IterT& beginChildren, const IterT& endChildren )
    {
	// add edges
	std::vector< NodeT > pres( preimage( parent ) );
	auto iParentInPres = std::find_if( pres.begin(), pres.end()
					   , std::bind( &RefinementTree< E >::equal, &*this, parent, std::placeholders::_1 ) );
	if( iParentInPres != pres.end() )
	    pres.erase( iParentInPres );
	
	std::vector< NodeT > posts( postimage( parent ) );
	auto iParentInPosts = std::find_if( posts.begin(), posts.end()
					    , std::bind( &RefinementTree< E >::equal, &*this, parent, std::placeholders::_1 ) );
	if( iParentInPosts != posts.end() )
	    posts.erase( iParentInPosts );

	for( auto ichild = beginChildren; ichild != endChildren; ++ichild )
	{
	    // connect parent's preimage minus self-loop
	    for( auto& pre : pres )
	    {
		if( possibly( isReachable( pre, *ichild ) ) )
		    graph::addEdge( mMapping, pre, *ichild );
	    }
	    // connect parent's post image minus self-loop
	    for( auto& post : posts )
	    {
		if( possibly( isReachable( *ichild, post ) ) )
		    graph::addEdge( mMapping, *ichild, post );
	    }
	    // connect all children from-to and to-from
	    for( auto iChildUp = ichild; iChildUp != endChildren; ++iChildUp )
	    {
		if( possibly( isReachable( *ichild, *iChildUp ) ) )
		    graph::addEdge( mMapping, *ichild, *iChildUp );
		if( ichild != iChildUp && possibly( isReachable( *iChildUp, *ichild ) ) ) // don't make a 2nd self-loop
		    graph::addEdge( mMapping, *iChildUp, *ichild );
	    }
	}
    }
    
    Ariadne::BoundedConstraintSet mSafeSet;
    Ariadne::EffectiveVectorFunction mDynamics;
    Ariadne::Effort mEffort;
    unsigned long mNodeIdCounter;
    MappingT mMapping;
    E mInitialEnclosure;
    NodeT mOutsideNode;
};

#endif

