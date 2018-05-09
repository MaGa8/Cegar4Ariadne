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

template< typename I >
class RefinementTree
{
  public:

    typedef I IntervalT;
    typedef Ariadne::Box< IntervalT > EnclosureT;

    //! \class value to store in refinement tree holding a box and safety variable
    class TreeValue
    {
      public:
	TreeValue()
	    : mEnclosure( Ariadne::Box< IntervalT >::zero( 0 ) )
	    , mSafe( false )
	{}
	
	TreeValue( const EnclosureT& e, Ariadne::ValidatedLowerKleenean safe )
	    : mEnclosure( e )
	    , mSafe( safe )
	{}

	TreeValue( const TreeValue& orig )
	    : mEnclosure( orig.mEnclosure )
	    , mSafe( orig.mSafe )
	{}

	TreeValue& operator =( const TreeValue& orig )
	{
	    this->mEnclosure = orig.mEnclosure;
	    this->mSafe = orig.mSafe;
	    return *this;
	}

	const EnclosureT& getEnclosure() const { return mEnclosure; }

	Ariadne::ValidatedLowerKleenean isSafe() const { return mSafe; }

	bool operator ==( const TreeValue& tv )
	{
	    Ariadne::LowerKleenean areEq = this->mEnclosure == tv.mEnclosure;
	    return definitely( areEq );
	}
      private:
	EnclosureT mEnclosure;
	Ariadne::ValidatedLowerKleenean mSafe;
    };

    typedef tree::LinkedFixedBranchTree< TreeValue, 2 > RefinementT;
    typedef graph::AdjacencyDiGraph< typename RefinementT::NodeT, graph::VecMap, graph::InVec, graph::InVec > MappingT;
    typedef typename MappingT::VertexT NodeT;

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
	, mRefinements( TreeValue( initial, constraints.overlaps( initial ).check( effort ) ) )
	, mEffort( effort )
    {
	NodeT initialNode = *addVertex( mMappings, root( mRefinements ) );
	if( possibly( isReachable( initialNode, initialNode ) ) )
	    addEdge( mMappings, initialNode, initialNode );
    }

    const Ariadne::ConstraintSet& getConstraints()
    {
	return mConstraints;
    }

    const Ariadne::EffectiveVectorFunction& getDynamics()
    {
	return mDynamics;
    }

    const RefinementT& getTree() const
    {
	return mRefinements;
    }

    const MappingT& getLeafMapping() const
    {
	return mMappings;
    }

    // TODO: implement const image function at some point
    //! \param from abstraction for which to find image in leaves of tree; needs to be of type that can be intersected with EnclosureT
    //! \return most refined boxes intersecting with from
    template< typename EnclosureT2 >
    std::vector< NodeT > image( const EnclosureT2& from )
    {
	std::vector< NodeT > constImage = imageRecursive( from, root( mRefinements ) );
	return constImage;
    }

    // helper function for recursive calls of image
    // abstract same code as leaves in DFS controller
    template< typename EnclosureT2 >
    std::vector< NodeT > imageRecursive( const EnclosureT2& from, const typename RefinementT::NodeT& to )
    {
    	// vector because there can't be duplicates
    	std::vector< NodeT > parts;
    	const EnclosureT& boxTo = tree::value( mRefinements, to ).getEnclosure();

    	auto inter = Ariadne::intersection( from, boxTo );

	// result is either Boolean for ExactInterval or LowerKleenean otherwise -> convert implicitly to latter if necessary

	//! \todo would like to use LowerKleenean only
	Ariadne::ValidatedLowerKleenean isEmpty = inter.is_empty();


	// if there's no chance that to and from intersect: return empty
	if( definitely( isEmpty/*.check( mEffort )*/ ) )
    	    return parts;
    	else if( tree::isLeaf( mRefinements, to ) )
    	{
    	    typename MappingT::VIterT iGraphNode = graph::findVertex( mMappings, to );
    	    if( iGraphNode == vertices( mMappings ).second )
    		throw std::logic_error( "no graph node could be found for a tree node" );
    	    parts.push_back( *iGraphNode );
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

    // all node accessors: can be declared const, because tree or graph need to be accessed in order to change anything
    //! \return all leaves of the tree
    std::vector< NodeT > leaves()
    {
	return leaves( tree::root( mRefinements ) );
    }
    
    //! \return all leaves in subtree at v
    std::vector< NodeT > leaves( const NodeT& v )
    {
	return leaves( graph::value( mMappings, v ) );
    }
    
    //! \return all leaves in subtree at subRoot
    std::vector< NodeT > leaves( const typename RefinementT::NodeT& treev ) 
    {
    	std::vector< NodeT > ls;
	
    	if( tree::isLeaf( mRefinements, treev ) )
    	{
    	    typename MappingT::VIterT iv = graph::findVertex( mMappings, treev );
    	    if( iv == vertices( mMappings ).second )
    		throw std::logic_error( "no matching graph vertex found for tree vertex" );
    	    else
    		ls.push_back( *iv );
    	}
    	else
    	{
    	    typename tree::FixedBranchTreeTraits< RefinementT >::CRangeT cs = tree::children( mRefinements, treev );
    	    for( ; cs.first != cs.second; ++cs.first )
    	    {
    		std::vector< NodeT > rls = leaves( *cs.first );
    		ls.insert( ls.end(), rls.begin(), rls.end() );
    	    }
    	}
    	return ls;
    }

    //! \return all leaves in refinement tree mapping to from
    std::vector< NodeT > preimage( const NodeT& from ) 
    {
	std::vector< NodeT > preimg;
	typename graph::DiGraphTraits< MappingT >::InRangeT ins = graph::inEdges( mMappings, from );
	for( ; ins.first != ins.second; ++ins.first )
	    preimg.push_back( graph::source( mMappings, *ins.first ) );

	return preimg;
    }

    //! \return all leaves in refinement tree from maps to
    std::vector< NodeT > postimage( const NodeT& from )
    {
	std::vector< NodeT > postimg;
	typename graph::DiGraphTraits< MappingT >::OutRangeT outs = graph::outEdges( mMappings, from );
	for( ; outs.first != outs.second; ++outs.first )
	    postimg.push_back( graph::target( mMappings, *outs.first ) );
	return postimg;
    }

    const TreeValue& getNodeValue( const NodeT& v ) const
    {
	return tree::value( mRefinements, graph::value( mMappings, v ) );
    }

    //! \return true if v is safe with certainty, false otherwise
    bool isSafe( const NodeT& v ) const
    {
    	return definitely( getNodeValue( v ).isSafe() );
    }

    Ariadne::ValidatedLowerKleenean nodesEqual( const NodeT& n1, const NodeT& n2 ) const
    {
	Ariadne::LowerKleenean areEqual = tree::value( mRefinements, graph::value( mMappings, n1 ) ).getEnclosure() ==
	    tree::value( mRefinements, graph::value( mMappings, n2 ) ).getEnclosure();
	return areEqual.check( mEffort );
    }

    //! \return true if trg is deemed reachable from src
    // \todo eventually parametrize this
    // \todo does this always return a validated Kleenean? test with effective boxes at some point
    Ariadne::ValidatedKleenean isReachable( const NodeT& src, const NodeT& trg ) const
    {
	const EnclosureT& srcBox = tree::value( mRefinements, graph::value( mMappings, src ) ).getEnclosure()
	    , trgBox = tree::value( mRefinements, graph::value( mMappings, trg ) ).getEnclosure();

	Ariadne::Point< Ariadne::Bounds< Ariadne::FloatDP > > mappedCentre = mDynamics.evaluate( srcBox.centre() );
    	return trgBox.contains( mappedCentre );
    }

    /*! 
      \param v leaf node in refinement tree
      \brief refines node v using r and updates
    */
    void refine( NodeT& v, const IRefinementStrategy< IntervalT >& r )
    {
    	typename RefinementT::NodeT treev = graph::value( mMappings, v );
    	// handle case: if v is not a leaf
	
    	// if v is leaf:
    	// obtain refinements of node v as list of appropriate length
    	EnclosureT obox = tree::value( mRefinements, treev ).getEnclosure();

    	std::vector< EnclosureT > refined = r.refine( obox );
    	// std::vector< TreeValue > tvals;
    	std::array< TreeValue, RefinementT::N > tvals;

	std::transform( refined.begin(), refined.end(), tvals.begin()
    			, [this] (const EnclosureT& e) { return TreeValue( e, mConstraints.overlaps( e ).check( mEffort ) ); } );

	tree::expand( mRefinements, treev, tvals );
    	// add expansions to the graph
    	std::pair< typename RefinementT::CIterT, typename RefinementT::CIterT > cs = tree::children( mRefinements, treev );
    	std::vector< NodeT > refinedNodes;
    	for( ; cs.first != cs.second; ++cs.first )
    	{
	    typename MappingT::VIterT ivadd = graph::addVertex( mMappings, *cs.first );
    	    refinedNodes.push_back( *ivadd );
    	}

    	// add edges
	// NO! NONE OF THIS WORKS! EDGES MAY APPEAR "OUT OF NOTHING" DUE TO THE WAY BOXES ARE SAID TO MAP INTO EACH OTHER! --->
	// std::vector< NodeT > pres( preimage( v ) );
	// std::vector< NodeT > posts( postimage( v ) );
	// // append all children of refined node, because center mapping refinements may introduce new links absent at parent level
	// // simply adding v is fine, because pres and posts will be traced down to leaf
	// pres.push_back( v );
	// posts.push_back( v );
	
    	// for( NodeT& refined : refinedNodes )
    	// {
    	//     // determine which elements of the preimage of v map to which refined component
    	//     for( NodeT& pre : pres )
    	//     {
    	// 	for( NodeT& preLeaf : leaves( pre ) )
    	// 	{
    	// 	    if( possibly( isReachable( preLeaf, refined ) ) )
	// 		graph::addEdge( mMappings, preLeaf, refined );
    	// 	}
    	//     }
    	//     // determine which elements of the postimage the components of v map to
    	//     for( NodeT& post : posts )
    	//     {
    	// 	for( NodeT& postLeaf : leaves( post ) )
    	// 	{
	// 	    if( possibly( isReachable( refined, postLeaf ) ) )
	// 		graph::addEdge( mMappings, refined, postLeaf );
    	// 	}
    	//     }
    	// }
    	// unlink v from the graph
	// <--- inefficient alternative below
	// to bring this back to life: map whole box and test image as intersection

	auto lvs = leaves( tree::root( mRefinements ) );
	for( NodeT& refined : refinedNodes )
	{
	    for( NodeT& lf : lvs )
	    {
		if( possibly( isReachable( lf, refined ) ) )
		    graph::addEdge( mMappings, lf, refined );
		if( possibly( isReachable( refined, lf ) ) )
		    graph::addEdge( mMappings, refined, lf );
	    }
	}
	
	graph::removeVertex( mMappings, v );
    }

  private:
    Ariadne::ConstraintSet mConstraints;
    Ariadne::EffectiveVectorFunction mDynamics;
    Ariadne::Effort mEffort;
    RefinementT mRefinements;
    MappingT mMappings;
};

/*!
  SHOULDN'T I USE REFINEMENT TREE NODES RATHER THAN BARE BONES TREE NODES?

  runs DFS to find counterexample
  any path terminates in
  1) loop leading back to state along path
  2) state with violated safety conditions
  \param initialBegin iterator to beginning of refinement tree nodes describing the image of the initial set, should dereference to box
  \return vector of nodes terminated by an unsafe node in REVERSE ORDER (terminal state first)
  \todo add parameter to control ordering of branches in dfs exploration 
  \todo speed up: use references in vector storing path passed to recursive calls -> don't need to copy nodes
  \todo use const refinement tree once constness issues in tree+graph are fixed
*/
template< typename IntervalT, typename NodeIterT >
std::vector< typename RefinementTree< IntervalT >::NodeT > findCounterexample( RefinementTree< IntervalT >& rtree
									       , NodeIterT iImgBegin, NodeIterT iImgEnd
									       , const std::vector< typename RefinementTree< IntervalT >::NodeT >& path = {} )
{
    for( ; iImgBegin != iImgEnd; ++iImgBegin )
    {
        // look for loop										       
	typename std::vector< typename RefinementTree< IntervalT >::NodeT >::const_iterator iBeginLoop =
	    std::find_if( path.begin(), path.end()
			  , [&rtree, &iImgBegin] (const typename RefinementTree< IntervalT >::NodeT& n) {return definitely( rtree.nodesEqual( *iImgBegin, n ) ); } );
	// loop found
	if( iBeginLoop != path.end() )
	{
	    return {};
	}

	std::vector< typename RefinementTree< IntervalT >::NodeT > copyPath( path.begin(), path.end() );
	copyPath.push_back( *iImgBegin );
	// counterexample found
	if( !Ariadne::definitely( rtree.isSafe( *iImgBegin ) ) )
	{
	    std::cout << "counterexample of length " << copyPath.size() << std::endl;
	    return copyPath;
	}

	// recurse & return
	auto posts =  rtree.postimage( *iImgBegin );
	std::vector< typename RefinementTree< IntervalT >::NodeT > cex = findCounterexample( rtree, posts.begin(), posts.end(), copyPath );
	if( !cex.empty() )
	    return cex;
    }
    // no counterexample found
    return {};
}

// implement this using lower kleenean?
/*! 
  \todo enforce constness of refinement tree
  \param beginCounter and endCounter iterators to beginning and end of counterexample trajectory, should dereference to typename RefinementTree< IntervalT >::NodeT
  \param beginImage and endImage iterators to beginning and end of image obtained from refinement tree
  \return false if there definitely exists a point that is mapped to the terminal state of the counterexample, true otherwise, including if there does not possibly exist such a point 
 */
template< typename IntervalT, typename PathIterT, typename ImageIterT >
Ariadne::ValidatedLowerKleenean isSpurious( RefinementTree< IntervalT >& rtree
					    , PathIterT beginCounter, PathIterT endCounter
					    , ImageIterT beginImage, ImageIterT endImage
					    , const Ariadne::Effort& effort )
{
    // intersect image with first element of counterex trajectory -> new image -> recurse
    std::vector< typename RefinementTree< IntervalT >::NodeT > intersectionImage;
    for( ; beginImage != endImage; ++beginImage )
    {
	Ariadne::UpperBoxType ubMapped =  Ariadne::image( rtree.getNodeValue( *beginImage ).getEnclosure(), rtree.getDynamics() );
	auto mapIntersection = Ariadne::intersection( ubMapped, rtree.getNodeValue( *beginCounter ).getEnclosure() );
	// certainly cannot trace any given set -> certainly spurious
	Ariadne::ValidatedLowerKleenean doesInter = mapIntersection.is_empty();
	if( !possibly( doesInter ) )
	    return true;

	Ariadne::ValidatedLowerKleenean doesContinue = true;
	if( beginCounter != endCounter )
	{
	    auto imageMapped = rtree.image( ubMapped );
	    doesContinue = isSpurious( rtree, beginCounter + 1, endCounter, imageMapped.begin(), imageMapped.end(), effort );
	}

	Ariadne::ValidatedLowerKleenean doesInterAndContinue = doesInter && doesContinue;
	// can surely trace some set -> certainly not spurious
	if( definitely( doesInterAndContinue ) )
	    return false;
    }
    // neither found instance to prove nor could disprove that the counterexample was spurious
    return true;
}

// can only prove that there exists a true counterexample -> system is unsafe
template< typename IntervalT, typename ImageIterT >
std::vector< typename RefinementTree< IntervalT >::NodeT > cegar( RefinementTree< IntervalT >& rtree
								  , ImageIterT imageBegin, ImageIterT imageEnd
								  , const Ariadne::Effort& effort
								  , const IRefinementStrategy< IntervalT >& refinementStrat
								  , const uint maxNodes )
{
    while( rtree.getTree().size() < maxNodes )
    {
	std::cout << "new iteration, number of nodes " << maxNodes << std::endl;
	// look for counterexample
	auto counterexample = findCounterexample( rtree, imageBegin, imageEnd );
	if( counterexample.empty() )
	    return {};
	if( !isSpurious( rtree, counterexample.begin(), counterexample.end(), imageBegin, imageEnd, effort ) )
	    return counterexample;

	// don't want to refine last state
	for( uint i = 0; i < counterexample.size() - 1; ++i )
	    rtree.refine( counterexample[ i ], refinementStrat );
    }
    return {};
}

#endif
