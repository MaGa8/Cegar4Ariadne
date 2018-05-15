#ifndef DI_GRAPH_INTERFACE_HPP
#define DI_GRAPH_INTERFACE_HPP

#include <utility>
#include <vector>

namespace graph
{

    /*!
      Template for directed graphs on unique values as keys to vertices
    */
    template< typename G >
    struct DiGraphTraits
    {
	// things
	// non const should be convertible to const
	typedef typename G::ValueT ValueT;
	// typedef typename G::ConstValueT ConstValueT
	typedef typename G::VertexT VertexT;
	// typedef typename G::ConstVertex ConstVertex
	typedef typename G::EdgeT EdgeT;
	// typedef typename G::ConstEdgeT ConstEdgeT;
	// iterators: should comply with ForwardIterator
	// should be iterators to const => any change in graph routed through graph class
	typedef typename G::VIterT VIterT;
	// typedef typename G::ConstVIterT ConstVIterT;
	// edge iterators should dereference to EdgeT
	typedef typename G::OutIterT OutIterT;
	// typedef typename G::ConstOutIterT ConstOutIterT;
	typedef typename G::InIterT InIterT;
	// typedef typename G::ConstInIterT ConstInIterT;
	//ranges
	template< typename IterT >
	using RangeT = std::pair< IterT, IterT >;

	typedef RangeT< VIterT > VRangeT;
	// typedef RangeT< ConstVIterT > ConstVRangeT;
	typedef RangeT< OutIterT > OutRangeT;
	// typedef RangeT< ConstOutIterT > ConstOutRangeT;
	typedef RangeT< InIterT > InRangeT;
	// typedef RangeT< ConstInIterT > ConstInRangeT;
    };

    //! \return a reference to the value stored in v
    //! \note values are not meant to be modified as nodes are accessed by their values
    template< typename G >
    const typename DiGraphTraits< G >::ValueT& value( const G& g, const typename DiGraphTraits< G >::ConstVertexT& v );

    template< typename G >
    typename DiGraphTraits< G >::VertexT source( const G& g, typename DiGraphTraits< G >::EdgeT& e );

    // template< typename G >
    // const typename DiGraphTraits< G >::ConstVertexT source( const G& g, const typename DiGraphTraits< G >::ConstEdgeT& e )
    // {
    // 	auto& eMutable = static_cast< typename DiGraphTraits< G >::EdgeT& >( const_cast< typename DiGraphTraits< G >::ConstEdgeT& >( e ) );
    // 	typename DiGraphTraits< G >::VertexT& nonConstSrc = source( const_cast< G& >( g ), eMutable );
    // 	return return nonConstSrc;
    // 	// const_cast< const typename DiGraphTraits< G >::ConstVertexT& >
    // 	//     ( static_cast< typename DiGraphTraits< G >::ConstVertexT& >( nonConstSrc ) );
    // }

    template< typename G >
    typename DiGraphTraits< G >::VertexT target( const G& g, typename DiGraphTraits< G >::EdgeT& e );

    // template< typename G >
    // const typename DiGraphTraits< G >::ConstVertexT target( const G& g, const typename DiGraphTraits< G >::EdgeT& e )
    // {
    // 	auto& eMutable = static_cast< typename DiGraphTraits< G >::EdgeT& >( const_cast< typename DiGraphTraits< G >::ConstEdgeT& >( e ) );
    // 	typename DiGraphTraits< G >::VertexT& nonConstTrg = target( const_cast< G& >( g ), eMutable );
    // 	return nonConstTrg;
    // 	// return const_cast< const typename DiGraphTraits< G >::ConstVertexT& >
    // 	//     ( static_cast< typename DiGraphTraits< G >::ConstVertexT& >( nonConstTrg ) );
    // }

    //! \return the vertices as a vertex range
    template< typename G >
    typename DiGraphTraits< G >::VRangeT vertices( const G& g );

    // //! \return the vertices as a const vertex range
    // template< typename G >
    // typename DiGraphTraits< G >::ConstVRangeT vertices( const G& g )
    // {
    // 	typename DiGraphTraits< G >::VRangeT nonConstR = vertices( const_cast< G& >( g ) );
    // 	return nonConstR;
    // }

    //! \return the outgoing edges of v as a out edge range
    template< typename G >
    typename DiGraphTraits< G >::OutRangeT outEdges( const G& g, const typename DiGraphTraits< G >::VertexT& v );

    // //! \return the outgoing edges of v as a const out edge range
    // template< typename G >
    // typename DiGraphTraits< G >::ConstOutRangeT outEdges( const G& g, const typename DiGraphTraits< G >::VertexT& v )   
    // {
    // 	auto& vMutable = static_cast< typename DiGraphTraits< G >::VertexT& >( const_cast< typename DiGraphTraits< G >::ConstVertexT& >( v ) );
    // 	typename DiGraphTraits< G >::OutRangeT nonConstR = outEdges( const_cast< G& >( g ), vMutable );
    // 	// return const_cast< const typename DiGraphTraits< G >::ConstOutRange& >
    // 	//     ( static_cast< typename DiGraphTraits< G >::ConstOutRangeT& >( nonConstR ) );
    // 	return nonConstR;
    // }

    //! \return the incoming edges of v as a in edge range
    template< typename G >
    typename DiGraphTraits< G >::InRangeT inEdges( const G& g, const typename DiGraphTraits< G >::VertexT& v );

    // //! \return the incoming edges of v as a const in edge range
    // template< typename G >
    // typename DiGraphTraits< G >::ConstInRangeT inEdges( const G& g, const typename DiGraphTraits< G >::VertexT& v )
    // {
    // 	auto& vMutable = static_cast< typename DiGraphTraits< G >::VertexT& >( const_cast< typename DiGraphTraits< G >::ConstVertexT& >( v ) );
    // 	typename DiGraphTraits< G >::InRangeT nonConstR = inEdges( const_cast< G& >( g ), vMutable );
    // 	return nonConstR;
    // 	// return const_cast< const typename DiGraphTraits< G >::ConstInRange& >
    // 	//     ( static_cast< typename DiGraphTraits< G >::ConstInRangeT& >( nonConstR ) );
    // }

    template< typename G >
    typename DiGraphTraits< G >::VIterT findVertex( const G& g, const typename DiGraphTraits< G >::ValueT& val );

    // template< typename G >
    // typename DiGraphTraits< G >::ConstVIterT findVertex( const G& g, const typename DiGraphTraits< G >::ConstValueT& val )
    // {
    // 	typename DiGraphTraits< G >::VIterT nonConstV = findVertex( const_cast< G& >( g ), val );
    // 	return nonConstV;
    // }

    //! \return out iterator to edge from s to t or one past end if none exists
    template< typename G >
    typename DiGraphTraits< G >::OutIterT findEdgeTo( const G& g
						      , const typename DiGraphTraits< G >::VertexT& s
						      , const typename DiGraphTraits< G >::VertexT& t );

    // //! \return const out iterator to edge from s to t or one past end if none exists
    // template< typename G >
    // typename DiGraphTraits< G >::ConstOutIterT findEdgeTo( const G& g
    // 							   , const typename DiGraphTraits< G >::ConstVertexT& s
    // 							   , const typename DiGraphTraits< G >::ConstVertexT& t )
    // {
    // 	auto nonConstS = static_cast< typename DiGraphTraits< G >::VertexT& >( const_cast< typename DiGraphTraits< G >::ConstVertexT& >( s ) );
    // 	auto nonConstT = static_cast< typename DiGraphTraits< G >::VertexT& >( const_cast< typename DiGraphTraits< G >::ConstVertexT& >( t ) );
    // 	typename DiGraphTraits< G >::OutIterT nonConstOut = findEdgeTo( const_cast< G& >( g ), nonConstS, nonConstT );
    // 	return nonConstOut;
    // }

    //! \return in iterator to edge into t from s or one past end if none exists
    template< typename G >
    typename DiGraphTraits< G >::InIterT findEdgeFrom( G& g
						       , const typename DiGraphTraits< G >::VertexT& s
						       , const typename DiGraphTraits< G >::VertexT& t );

    // //! \return const in iterator to edge into t from s or one past end if none exists
    // template< typename G >
    // typename DiGraphTraits< G >::ConstInIterT findEdgeFrom( const G& g
    // 							    , const typename DiGraphTraits< G >::ConstVertexT& s
    // 							    , const typename DiGraphTraits< G >::ConstVertexT& t )
    // {
    // 	auto nonConstS = static_cast< typename DiGraphTraits< G >::VertexT& >( const_cast< typename DiGraphTraits< G >::ConstVertexT& >( s ) );
    // 	auto nonConstT = static_cast< typename DiGraphTraits< G >::VertexT& >( const_cast< typename DiGraphTraits< G >::ConstVertexT& >( t ) );
    // 	typename DiGraphTraits< G >::OutIterT nonConstIn = findEdgeFrom( const_cast< G& >( g ), nonConstS, nonConstT );
    // 	return nonConstIn;
    // }

    //! \brief add new node for value v to g assuming v is unique
    //! \return iterator to vertex added
    template< typename G >
    typename DiGraphTraits< G >::VIterT& addVertex( G& g, const typename DiGraphTraits< G >::ValueT& v );

    //! \brief add edge between s and t or keep existing edge
    //! \return pair of out and in iterators to edges added
    template< typename G >
    std::pair< typename DiGraphTraits< G >::OutEdgeT
	       , typename DiGraphTraits< G >::InEdgeT > addEdge( G& g
								 , const typename DiGraphTraits< G >::VertexT& s
								 , const typename DiGraphTraits< G >::VertexT& t );

    //! \brief remove vertex v and invalidate it after the operation
    //! \return iterator to vertex after removed vertex
    template< typename G >
    typename DiGraphTraits< G >::VIterT removeVertex( G& g, const typename DiGraphTraits< G >::VertexT& v );

    //! \brief remove edge between vertices s and t or do nothing if none exists
    //! \return pair of out and in iterators past edges removed
    template< typename G >
    std::pair< typename DiGraphTraits< G >::OutIterT
	       , typename DiGraphTraits< G >::InIterT > removeEdge( G& g
								    , const typename DiGraphTraits< G >::EdgeT& e );

    //! \brief remove edge between vertices s and t or do nothing if none exists
    //! \return pair of out and in iterators past edges removed
    template< typename G >
    std::pair< typename DiGraphTraits< G >::OutIterT
	       , typename DiGraphTraits< G >::InIterT > removeEdge( G& g
								    , const typename DiGraphTraits< G >::VertexT& s
								    , const typename DiGraphTraits< G >::VertexT& t );
}

#endif
