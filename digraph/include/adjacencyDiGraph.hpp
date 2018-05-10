#ifndef ADJACENCY_DI_GRAPH_HPP
#define ADJACENCY_DI_GRAPH_HPP

#include "diGraphInterface.hpp"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <functional>
#include <cmath>

#include <assert.h>

#define ADJ_GRAPH_TEMPLATE template< typename T,    template< typename K, typename V > class VCT,    template< typename S > class OECT,    template< typename S > class IECT >

#define AGRAPH AdjacencyDiGraph< T, VCT, OECT, IECT >

namespace graph
{

    //! \class insert vector: provide insert and find for OECT, IECT edge containers
    template< typename T >
    class InVec : public std::vector< T >
    {
      public:
	typedef std::vector< T > BaseT;

	InVec() {}

	InVec( typename BaseT::const_iterator istart
	       , typename BaseT::const_iterator iend )
	    : BaseT( istart, iend )
	{}

	typename BaseT::const_iterator insert( const T& val )
	{
	    BaseT::push_back( val );
	    return BaseT::end() - 1;
	}

	typename BaseT::const_iterator find( const T& val ) const
	{
	    return std::find( BaseT::begin(), BaseT::end(), val );
	}
    };

    //! \class vector map: simpler replacement for vertex map
    template< typename K, typename V  >
    class VecMap : public std::vector< std::pair< K, V > >
    {
      public:
	typedef std::vector< std::pair< K, V > > BaseT;
	typedef K key_type;
	typedef V mapped_type;
	typedef const mapped_type& reference;
	typedef const mapped_type* pointer;

	struct ValueIterator : public BaseT::const_iterator
	{
	    ValueIterator() = default;

	    ValueIterator( const typename BaseT::const_iterator& i ) : BaseT::const_iterator( i ) {}

	    ValueIterator( const ValueIterator& ci ) : BaseT::const_iterator( ci ) {}

	    ValueIterator& operator =( const ValueIterator& ci ) {  BaseT::const_iterator::operator =( ci ); return *this; }

	    const mapped_type& operator *() { return (BaseT::const_iterator::operator *()).second; }

	    const mapped_type* operator ->() { return &(BaseT::const_iterator->second); } // hope this is not a bug due to temp address
	};

	VecMap() {}
	
	VecMap( typename BaseT::const_iterator start
		, typename BaseT::const_iterator end )
	    : BaseT( start, end )
	{}

	ValueIterator find( const K& key ) const
	{
	    return ValueIterator( std::find_if( BaseT::begin(), BaseT::end()
						, [&key] (const std::pair< key_type, mapped_type >& kv ) { return kv.first == key; } ) );
	}
	
	ValueIterator insert( const K& key, const V& val )
	{
	    BaseT::push_back( std::make_pair( key, val ) );
	    return ValueIterator( BaseT::end() - 1 );
	}

	ValueIterator erase( const K& key )
	{
	    return ValueIterator( BaseT::erase( this->find( key ) ) );
	}
    };

    /*!
      directed graph implemented by means of adjacency containers for incoming and outgoing edges
      values stored are used as a basis for comparison and therefore are immutable
      \param T type to be stored as value, has to be copy constructible, assignable and comparable w.r.t. equality
      \param VCT container type for storing the vertices; should support interface alike to std::map with methods insert, at, find and whose iterators should dereference to key-value pairs
      \param OECT container template type supporting find, insert and remove to use for storing out-edges and whose contents can be modified using iterators
      \param IECT container template type supporting find, insert and remove to use for storing in-edges and whose contents can be modified using iterators

      \note should be Container::iterator insert( const T& x ) and erase( Container::iterator i )
    */
    template< typename T
	      , template< typename K, typename V > class VCT
	      , template< typename S > class OECT
	      , template< typename S > class IECT >
    class AdjacencyDiGraph
    {
      private:
	struct InternalNode;
      public:
	// class ConstEdge;
	struct Edge;
	// class ConstNode;
	struct Node;
	// template< typename WrapIterT, typename CustomConstT >
	// struct ConstIteratorAdapterTemplate;

	// things
	typedef AdjacencyDiGraph< T, VCT, OECT, IECT > Agraph;
	// typedef ConstT ConstValueT;
	typedef T ValueT;
	// typedef ConstNode ConstVertexT;
	typedef Node VertexT;
	// typedef ConstEdge ConstEdgeT;
	typedef Edge EdgeT;

	// containers
	typedef VCT< T, Node > VertexContainerT;
	typedef OECT< EdgeT > OutEdgeContainerT;
	typedef IECT< EdgeT > InEdgeContainerT;

	// iterators
	// required because VertexT is not just value type of VCT
	typedef typename VertexContainerT::ValueIterator VIterT;
	// typedef typename VertexContainerT::const_iterator ConstVIterT;
	// typedef typename ConstIteratorAdapterTemplate< VertexContainerT::iterator, ConstVertexT > ConstVIterT;
	typedef typename OutEdgeContainerT::const_iterator OutIterT;
	// typedef typename OutEdgeContainerT::const_iterator ConstOutIterT;
	// typedef typename ConstIteratorAdapterTemplate< OutEdgeContainerT::iterator, ConstEdgeT > ConstOutIterT;
	typedef typename InEdgeContainerT::const_iterator InIterT;
	// typedef typename InEdgeContainerT::const_iterator ConstInIterT;
	// typedef typename ConstIteratorAdapterTemplate< InEdgeContainerT::iterator, ConstEdgeT > ConstInIterT;

	// represent node
	// class ConstNode
	// {
	//     friend class AdjacencyDiGraph< T, VCT, OECT, IECT, ConstT >;
	//   public:
	//     ConstNode() = default;

	//     ConstNode( std::shared_ptr< InternalNode > pnode ) : mPtr( pnode ) {}

	//     ConstNode( const ConstNode& orig ) : mPtr( orig.mPtr ) {}

	//     ConstNode& operator =( const ConstNode& orig ) { this->mPtr = orig.mPtr; }

	//     bool operator ==( const ConstNode& other ) { return *mPtr == orig.mPtr; }
	//   private:
	//     std::shared_ptr< InternalNode > mPtr;
	// };

	class Node
	{
	    friend class AdjacencyDiGraph< T, VCT, OECT, IECT >;
	  public:
	    Node() = default;

	    Node( std::shared_ptr< InternalNode > pnode ) : mPtr( pnode ) {}

	    Node( const Node& orig ) : mPtr( orig.mPtr ) {}

	    Node& operator =( const Node& orig ) { this->mPtr = orig.mPtr; return *this; }

	    bool operator ==( const Node& other ) const { return *mPtr == *other.mPtr; }
	  private:
	    std::shared_ptr< InternalNode > mPtr;
	};

	class Edge
	{
	    friend class AdjacencyDiGraph< T, VCT, OECT, IECT >;
	  public:
	    Edge() = default;
	    
	    Edge( const VertexT& src, const VertexT& trg ) : mSource( src ), mTarget( trg ) {}

	    Edge( const Edge& orig ) = default;

	    Edge& operator =( const Edge& orig ) = default;

	    bool operator ==( const Edge& other ) const
	    {
		return this->mSource == other.mSource &&
		    this->mTarget == other.mTarget;
	    }
	  private:
	    VertexT mSource, mTarget;
	};

	// struct Edge : public ConstEdge
	// {
	//     Edge() = default;

	//     Edge( VertexT& src, VertexT& trg ) : ConstEdge( src, trg ) {}

	//     Edge( const Edge& orig ) : ConstEdge( orig ) {}

	//     Edge& operator =( const Edge& orig ) { ConstEdge::operator =( orig ); }
	// };

	//unpacking methods: const stuff is done by casting in interface functions
	const ValueT& value( const VertexT& v ) const
	{
	    return v.mPtr->mValue;
	}
	
	typename DiGraphTraits< Agraph >::VRangeT vertices() const
	{
	    return std::make_pair( mVertices.begin(), mVertices.end() );
	}

	VIterT findVertex( const ValueT& v ) const
	{ 
	    return mVertices.find( v );
	}

	typename DiGraphTraits< Agraph >::OutRangeT outEdges( const VertexT& v ) const
	{
	    return std::make_pair( v.mPtr->mOuts.begin(), v.mPtr->mOuts.end() );
	}

	//! \return edge in src to trg
	OutIterT findEdgeTo( const VertexT& src, const VertexT& trg ) const
	{
	    return src.mPtr->mOuts.find( EdgeT( src, trg ) );
	}

	typename DiGraphTraits< Agraph >::InRangeT inEdges( const VertexT& v ) const
	{
	    return std::make_pair( v.mPtr->mIns.begin(), v.mPtr->mIns.end() );
	}
	
	//! \rturn edge in trg from src
	InIterT findEdgeFrom( const VertexT& src, const VertexT& trg ) const
	{
	    return trg.mPtr->mIns.find( EdgeT( src, trg ) );
	}

	const VertexT& source( const EdgeT& e ) const
	{
	    return e.mSource;
	}

	const VertexT& target( const EdgeT& e ) const
	{
	    return e.mTarget;
	}

	//! \return iterator to vertex added
	VIterT addVertex( const ValueT& v )
	{
	    VIterT ivAdd = mVertices.find( v );
	    if( ivAdd == mVertices.cend() )
		ivAdd = mVertices.insert( v, VertexT( std::shared_ptr< InternalNode >( new InternalNode( v ) ) ) );
	    return ivAdd; 
	}

	//! \return iterator to vertex following v
	VIterT removeVertex( const VertexT& v ) // v may be const, because v still points to the same vertex, only the context (graph) does
	{
	    const ValueT& vval = value( v );
	    // remove all edges leading to v
	    typename DiGraphTraits< AGRAPH >::InRangeT ins = inEdges( v );
	    for( InIterT in = ins.first; in != ins.second; ++in )
	    {
		const VertexT& src = source( *in );
		OutIterT iout = findEdgeTo( src, v );

		if( iout != outEdges( src ).second )
		    src.mPtr->mOuts.erase( iout ); // this should work, because const mPtr derefs to non const reference
	    }

	    // remove all edges going out from v
	    typename DiGraphTraits< AGRAPH >::OutRangeT outs = outEdges( v );
	    for( OutIterT out = outs.first; out != outs.second; ++out )
	    {
		const VertexT& trg = target( *out );
		InIterT iin = findEdgeFrom( v, trg );
		if( iin != inEdges( trg ).second )
		    trg.mPtr->mIns.erase( findEdgeFrom( v, trg ) );
	    }

	    v.mPtr->mIns.clear();
	    v.mPtr->mOuts.clear();

	    return mVertices.erase( vval );
	}

	std::pair< OutIterT, InIterT > addEdge( const VertexT& src, const VertexT& trg )
	{
	    const EdgeT newEdge( src, trg ) ;
	    OutIterT iout = findEdgeTo( src, trg );
	    InIterT iin = findEdgeFrom( src, trg );
	    if( iout == src.mPtr->mOuts.end() )
		iout = src.mPtr->mOuts.insert( newEdge );
	    if( iin == trg.mPtr->mIns.end() )
		iin = trg.mPtr->mIns.insert( newEdge );

	    return std::make_pair( iout, iin );
	}

	std::pair< OutIterT, InIterT > removeEdge( const VertexT& src, const VertexT& trg )
	{
	    OutIterT ieSrc = findEdgeTo( src, trg );
	    InIterT ieTrg = findEdgeFrom( src, trg );
	    
	    assert( (ieSrc == outEdges( src ).second && ieTrg == inEdges( trg ).second ) ||
		    (ieSrc != outEdges( src ).second && ieTrg != inEdges( trg ).second ) );

	    if( ieSrc != outEdges( src ).second )
		ieSrc = src.mPtr->mOuts.erase( ieSrc );
	    if( ieTrg != inEdges( trg ).second )
		ieTrg = trg.mPtr->mIns.erase( ieTrg );
	    return std::make_pair( ieSrc, ieTrg );
	}
	
      private:
	// store node
	struct InternalNode
	{
	    InternalNode() = default;

	    InternalNode( const ValueT& v ) : mValue( v ) {}
	    
	    InternalNode( const InternalNode& v )
		: mIns( v.mIns.begin(), v.mIns.end() )
		, mOuts( v.mOuts.begin(), v.mOuts.end() )
	    {}

	    InternalNode& operator =( const InternalNode& v ) = delete;

	    bool operator ==( const InternalNode& v ) { return this->mValue == v.mValue; }

	    ValueT mValue;
	    InEdgeContainerT mIns;
	    OutEdgeContainerT mOuts;
	};

	/* template for wrappers around Bidirectional iterators WrapIterT
	   dereferencing to custom const types
	 */
	// template< typename WrapIterT, typename DerefConstT >
	// struct ConstIteratorAdapterTemplate
	// {
	//   public:
	//     typedef ConstIteratorAdapterTemplate< WrapIterT, DerefConstT > ThisT;
	//     typedef typename std::iterator_traits< WrapIterT >::difference_type difference_type;
	//     typedef typename std::iterator_traits< WrapIterT >::iterator_category iterator_category;
	//     typedef typename DerefConstT value_type;
	//     typedef const DerefConstT& reference;
	//     typedef const DerefConstT& pointer;

	//     ConstIteratorAdapterTemplate() = default;

	//     ConstIteratorAdapterTemplate( const WrapIterT& iter ) : mIter( iter ) {}

	//     ConstIteratorAdapterTemplate( const ThisT& orig ) : mIter( orig.mIter ) {}

	//     ThisT& operator =( const ThisT& orig ) = default;

	//     ThisT& operator ++() { ++mIter; return *this; }

	//     ThisT& operator --() { --mIter; return *this; }

	//     ThisT& operator +=( const long int& add ) { std::advance( mIter, add ); return *this; }

	//     ThisT& operator -=( const long int& sub ) { std::advance( miter, -sub ); return *this; }

	//     reference operator *() const { return *mWrapped; }

	//     bool operator ==( const ThisT& other ) const { return this->mIter == other.mIter; }

	//     bool operator !=( const ThisT& other ) const { return this->mIter != other.mIter; }

	//     void swap( ThisT& other ) { swap( this->mIter, other.mIter ); }
	    
	//   protected:
	//     WrapIterT mIter;
	// };
	
	VertexContainerT mVertices;
    };

    // how to handle changing values when using e.g. set as underlying container type
    // where keys may not be altered after insertion
    // while for other containers, e.g. vector this is fine
    ADJ_GRAPH_TEMPLATE
    const typename AGRAPH::ValueT& value( const AGRAPH& ag, const typename AGRAPH::VertexT& v ) { return ag.value( v ); }

    ADJ_GRAPH_TEMPLATE
    const typename AGRAPH::VertexT& source( const AGRAPH& ag, const typename AGRAPH::EdgeT& e ) { return ag.source( e ); }

    ADJ_GRAPH_TEMPLATE
    const typename AGRAPH::VertexT& target( const AGRAPH& ag, const typename AGRAPH::EdgeT& e ) { return ag.target( e ); }

    ADJ_GRAPH_TEMPLATE
    typename DiGraphTraits< AGRAPH >::VRangeT vertices( const AGRAPH& ag ) { return ag.vertices(); }

    ADJ_GRAPH_TEMPLATE
    typename DiGraphTraits< AGRAPH >::OutRangeT outEdges( const AGRAPH& ag, const typename AGRAPH::VertexT& v ) { return ag.outEdges( v ); }

    ADJ_GRAPH_TEMPLATE
    typename DiGraphTraits< AGRAPH >::InRangeT inEdges( const AGRAPH& ag, const typename AGRAPH::VertexT& v ) { return ag.inEdges( v ); }

    ADJ_GRAPH_TEMPLATE
    typename AGRAPH::VIterT findVertex( const AGRAPH& ag, const typename AGRAPH::ValueT& val ) { return ag.findVertex( val ); }

    ADJ_GRAPH_TEMPLATE
    typename AGRAPH::OutIterT findEdgeTo( const AGRAPH& ag, const typename AGRAPH::VertexT& src, const typename AGRAPH::VertexT& trg ) { return ag.findEdgeTo( src, trg ); }

    ADJ_GRAPH_TEMPLATE
    typename AGRAPH::InIterT findEdgeFrom( const AGRAPH& ag, const typename AGRAPH::VertexT& src, const typename AGRAPH::VertexT& trg ) { return ag.findEdgeFrom( src, trg ); }

    ADJ_GRAPH_TEMPLATE
    typename AGRAPH::VIterT addVertex( AGRAPH& ag, const typename AGRAPH::ValueT& v ) { return ag.addVertex( v ); }

    ADJ_GRAPH_TEMPLATE
    std::pair< typename AGRAPH::OutIterT, typename AGRAPH::InIterT > addEdge( AGRAPH& ag, const typename AGRAPH::VertexT& src, const typename AGRAPH::VertexT& trg ) { return ag.addEdge( src, trg ); }

    /*! \note invalidates all VertexT, VIterT, VRangeT and constant derivatives */
    ADJ_GRAPH_TEMPLATE
    typename AGRAPH::VIterT removeVertex( AGRAPH& ag, const typename AGRAPH::VertexT& v ) { return ag.removeVertex( v ); }

    /*! \note invalidates all OutEdgeT of src and InEdgeT of trg as well as corresponding iterators and consts */
    ADJ_GRAPH_TEMPLATE
    std::pair< typename AGRAPH::OutIterT, typename AGRAPH::InIterT >  removeEdge( AGRAPH& ag, const typename AGRAPH::EdgeT& e ) { return ag.removeEdge( ag.source( e ), ag.target( e ) ); }

    /*! \note invalidates all OutEdgeT of src and InEdgeT of trg as well as corresponding iterators and consts */
    ADJ_GRAPH_TEMPLATE
    std::pair< typename AGRAPH::OutIterT, typename AGRAPH::InIterT > removeEdge( AGRAPH& ag, const typename AGRAPH::VertexT& src, const typename AGRAPH::VertexT& trg ) { return ag.removeEdge( src, trg ); }

    template< typename T
	      , template< typename K, typename V > typename VCT
	      , template< typename S > typename OECT
	      , template< typename S > typename IECT
	      , typename CharT, typename Traits
	      , typename P>
    std::basic_ostream< CharT, Traits >& print( std::basic_ostream< CharT, Traits >& os, const AGRAPH& ag, const typename AGRAPH::VertexT& v
						, const std::function< P( const typename AGRAPH::ValueT& ) >& converter = [] (const typename AGRAPH::ValueT& v) { return v; } )
    {
	const int segmentSize = 40;
	    
	typename DiGraphTraits< AGRAPH >::InRangeT ins = graph::inEdges( ag, v );
	typename DiGraphTraits< AGRAPH >::OutRangeT outs = graph::outEdges( ag, v );
	uint maxInOut = std::max( std::distance( ins.first, ins.second )
				  , std::distance( outs.first, outs.second ) );
	
	std::stringstream line;
	for( uint cEdge = 0; cEdge < maxInOut || cEdge == 0; ++cEdge )
	{
	    if( ins.first != ins.second )
		line << converter( graph::value( ag, graph::source( ag, *ins.first ) ) ) << " ->";
	    line << std::setw( std::max( 0, segmentSize - static_cast< const int&& >( line.str().size() ) ) );
	    if( cEdge == std::floor( maxInOut / 2.0 ) )
		line << converter( graph::value( ag, v ) );
	    line << std::setw( std::max( 0, 2 * segmentSize - static_cast< const int&& >( line.str().size() ) ) ) << " ";
	    if( outs.first != outs.second )
	    	line << "-> " << converter( graph::value( ag, graph::target( ag, *outs.first ) ) );
	    os << line.str() << std::endl;
	    line.str( std::string() );
	    
	    if( ins.first != ins.second )
		++ins.first;
	    if( outs.first != outs.second )
		++outs.first;
	}
	return os;
    }

    template< typename T
	      , template< typename K, typename V > typename VCT
	      , template< typename S > typename OECT
	      , template< typename S > typename IECT
	      , typename CharT, typename Traits
	      , typename P>
    std::basic_ostream< CharT, Traits >& print( std::basic_ostream< CharT, Traits >& os, const AGRAPH& ag
						, const std::function< P( const typename AGRAPH::ValueT& ) >& converter = [] (const typename AGRAPH::ValueT& v) {return v;} )
    {
	typename DiGraphTraits< AGRAPH >::VRangeT vs = graph::vertices( ag );
	for( ; vs.first != vs.second; ++vs.first )
	{
	    // std::cout << std::distance( vs.first, vs.second ) << "   :";
	    print( os, ag, *vs.first, converter );
	    os << std::endl;
	}
	return os;
    }

    // somehow this function causes a segfault
    template< typename T
	      , template< typename K, typename V > typename VCT
	      , template< typename S > typename OECT
	      , template< typename S > typename IECT
	      , typename CharT, typename Traits >
    std::basic_ostream< CharT, Traits >& operator <<( std::basic_ostream< CharT, Traits >& os, const AGRAPH& ag )
    {
	const uint offset = 4;
    
	typename DiGraphTraits< AGRAPH >::VRangeT vs = graph::vertices( ag );
	for( typename AGRAPH::VIterT iv = vs.first; iv != vs.second; ++iv )
	{
	    os << "vertex " << graph::value( ag, *iv ) << std::endl;
	    typename DiGraphTraits< AGRAPH >::InRangeT ins = graph::inEdges( ag, *iv );
	    typename DiGraphTraits< AGRAPH >::OutRangeT outs = graph::outEdges( ag, *iv );

	    for( typename AGRAPH::InIterT in = ins.first; in != ins.second; ++in )
		os << std::setw( offset ) << "<- " << graph::value( ag, graph::source( ag, *in ) ) << std::endl;
	    for( typename AGRAPH::OutIterT out = outs.first; out != outs.second; ++out )
		os << std::setw( offset ) << "-> " << graph::value( ag, graph::target( ag, *out ) ) << std::endl;
	    os << std::endl;
	}
	return os;
    }
}

#undef AGRAPH
#undef ADJ_GRAPH_TEMPLATE

#endif
