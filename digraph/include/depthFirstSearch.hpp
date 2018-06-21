#ifndef DEPTH_FIRST_SEARCH
#define DEPTH_FIRST_SEARCH

#include <unordered_set>
#include <functional>

namespace graph
{

    //! \brief do-nothing blueprint for depth first control classes
    //! \note const parameters may also be mutable if needed
    struct BidirectionalDFTControl
    {
	//! \brief called when tranversal is begun with initial node v
	template< typename DiGraphT >
	void init( const DiGraphT& g, const typename DiGraphT::VertexT& v ) {}

	//! \brief called when (in/out) edge is explored
	template< typename DiGraphT >
	void explore( const DiGraphT& g, const typename DiGraphT::EdgeT& e ) {}

	//! \brief called when node is visited after corresponding edge has been explored
	template< typename DiGraphT, typename VisitSetT >
	void visit( const DiGraphT& g, const typename DiGraphT::VertexT& v, const VisitSetT& visited ) {}

	//! \return true if traversal is to be terminated
	template< typename DiGraphT, typename VisitSetT >
	bool terminate( const DiGraphT& g, const typename DiGraphT::VertexT& v, const VisitSetT& visited ) { return false; } 

	//! \return true if traversal should backtrack on v before visiting it
	template< typename DiGraphT, typename VisitSetT >
	bool isBacktrack( const DiGraphT& g, const typename DiGraphT::VertexT& v, const VisitSetT& visited ) { return false; } 

	//! \brief called when backtracking on v after exploring corresponding edge
	template< typename DiGraphT, typename VisitSetT >
	void backtrack( const DiGraphT& g, const typename DiGraphT::VertexT& v, const VisitSetT& visited ) {} 

	//! \brief called when returning from exploring e
	template< typename DiGraphT >
	void returned( const DiGraphT& g, const typename DiGraphT::EdgeT& e ) {}

	//! \brief called when visit of v is done
	template< typename DiGraphT, typename VisitSetT >	
	void leave( const DiGraphT& g, const typename DiGraphT::VertexT& v, const VisitSetT& visited ) {} 

	//! \brief called when traversal finishes
	template< typename DiGraphT >
	void finish( const DiGraphT& g, const typename DiGraphT::VertexT& v ) {}
    };

    //! \brief base for control classes that support only forward depth first traversal
    struct ForwardDFTControl : public BidirectionalDFTControl {};

    //! \brief base for control classes that support only backward depth first traversal
    struct BackwardDFTControl : public BidirectionalDFTControl {};

    //! \brief helper to perform initial visit of any DFT
    template< typename DiGraphT, typename ControlT, typename UnorderedSet >
    void dftVisit( DiGraphT& g, const typename DiGraphT::VertexT& v, ControlT& c, UnorderedSet& visited )
    {
	visited.insert( v );
	c.visit( g, v, visited );
    }

    //! \brief helper to perform edge exploration
    //! \param e edge to explore
    //! \param v endpoint of e from which to explore
    //! \param u endpoint of e to which to explore
    template< typename DiGraphT, typename ControlT, typename UnorderedSetT >
    void dftExplore( DiGraphT& g, const typename DiGraphT::VertexT& v, const typename DiGraphT::VertexT& u, const typename DiGraphT::EdgeT& e
		     , ControlT& c, UnorderedSetT& visited
		     , const std::function< ControlT&( DiGraphT&, const typename DiGraphT::VertexT&, ControlT&, UnorderedSetT& ) >& recCall )
    {
	c.explore( g, e );
	if( c.isBacktrack( g, u, visited ) || visited.find( u ) != visited.end() )
	    c.backtrack( g, u, visited );
	else
	{
	    recCall( g, u, c, visited );
	    c.returned( g, e );
	}
    }

    /*!
      \brief performs dfs traversal forward
      \param g graph to traverse
      \param v node at which to start traversal
      \param c is notified of events in search and controls termination
      \param visited set of nodes visited
    */
    template< typename DiGraphT, typename ControlT, typename UnorderedSetT >
    ControlT& forwardDFTRecursive( DiGraphT& g, const typename DiGraphT::VertexT& v, ControlT& c, UnorderedSetT& visited )
    {
	std::function< ControlT&( DiGraphT&, const typename DiGraphT::VertexT&, ControlT&, UnorderedSetT& ) >
	    selfCall( forwardDFTRecursive< DiGraphT, ControlT, UnorderedSetT > );

	dftVisit( g, v, c, visited );

	for( auto iouts = outEdges( g, v ); iouts.first != iouts.second && !c.terminate( g, v, visited ); ++iouts.first )
	{
	    const typename DiGraphT::VertexT& t = target( g, *iouts.first );
	    dftExplore( g, v, t, *iouts.first, c, visited, selfCall );
	}

	c.leave( g, v, visited );
	return c;
    }

    /*!
      \brief performs dfs traversal backward
      \param g graph to traverse
      \param v node at which to start traversal
      \param c is notified of events in search and controls termination
      \param visited set of nodes visited
    */
    template< typename DiGraphT, typename ControlT, typename UnorderedSetT >
    ControlT& backwardDFTRecursive( DiGraphT& g, const typename DiGraphT::VertexT& v, ControlT& c, UnorderedSetT& visited )
    {
	std::function< ControlT&( DiGraphT&, const typename DiGraphT::VertexT&, ControlT&, UnorderedSetT& ) >
	    selfCall( backwardDFTRecursive< DiGraphT, ControlT, UnorderedSetT > );
	dftVisit( g, v, c, visited );

	for( auto iIn = inEdges( g, v ); iIn.first != iIn.second && !c.terminate( g, v, visited ); ++iIn.first )
	{
	    const typename DiGraphT::VertexT& s = source( g, *iIn.first );
	    dftExplore( g, v, s, *iIn.first, c, visited, selfCall );
	}

	c.leave( g, v, visited );
	return c;
    }

    template< typename DiGraphT, typename ControlT
	      , typename HashT = std::hash< typename DiGraphT::ValueT >
	      , typename CompareT = std::equal_to< typename DiGraphT::ValueT > >
    ControlT& forwardDFT( DiGraphT& g, const typename DiGraphT::VertexT& v, ControlT& c
			  , const HashT& hsh = HashT()
			  , const CompareT& cmp = CompareT() )
    {
	std::unordered_set< typename DiGraphT::VertexT, HashT, CompareT > visitSet( graph::size( g ), hsh, cmp );

	c.init( g, v );
	forwardDFTRecursive( g, v, c, visitSet );
	c.finish( g, v );
	return c;
    }

    template< typename DiGraphT, typename ControlT
	      , typename HashT = std::hash< typename DiGraphT::ValueT >
	      , typename CompareT = std::equal_to< typename DiGraphT::ValueT > >
    ControlT& backwardDFT( DiGraphT& g, const typename DiGraphT::VertexT& v, ControlT& c
				    , const HashT& hsh = HashT()
				    , const CompareT& cmp = CompareT() )
    {
	std::unordered_set< typename DiGraphT::VertexT, HashT, CompareT > visitSet( graph::size( g ), hsh, cmp );

	c.init( g, v );
	backwardDFTRecursive( g, v, c, visitSet );
	c.finish( g, v );
	return c;
    }

} // namespace

#endif // incl guard
