#ifndef LINKED_N_TREE_HPP
#define LINKED_N_TREE_HPP

#include "fixedBranchTreeInterface.hpp"
#include "childrenIteratorInterface.hpp"

#include <array>
#include <memory>
#include <algorithm>
#include <exception>
#include <iterator>
#include <iostream>
#include <functional>

#define TREE_TEMPLATE template< typename T, size_t NC >
#define LTREE LinkedFixedBranchTree< T, NC >

namespace tree
{
    /*!
      \brief Linked n-ary tree
      \param T type to store in nodes
      \param N number of children per interior node
      \param ConstT custom const type, T should be implicitly convertible to ConstT
      \note rename NodeT to node NodeProxy
    */
    TREE_TEMPLATE
    class LinkedFixedBranchTree : public FixedBranchTreeInterface< T, NC >
    {
      public:
	// forward decls
	struct LinkedInternalNode;
	struct Node;
	// struct ConstNode;
	// struct ConstChildrenIterator;
	
	static const size_t N = NC;
	// trait defs
	typedef T ValueT;
	// typedef ConstT ConstValueT;
	typedef Node NodeT;
	// typedef ConstNode ConstNodeT;
	typedef std::array< NodeT, NC >CListT;
	// typedef ConstChildrenIterator ConstCIterT;
	typedef typename CListT::const_iterator CIterT;
    
	//! \return assignable and copyable node for emulating a const reference
	class Node
	{
	    friend class LinkedFixedBranchTree< T, NC >;
	  public:
	    Node() {}

	    Node( const std::shared_ptr< LinkedInternalNode >& nodePtr ) : mPtr( nodePtr ) {}

	    Node( const Node& orig ) = default;

	    Node& operator =( const Node& orig ) = default;

	    bool operator ==( const Node& other ) const
	    {
		return *(this->mPtr) == *(other.mPtr);
	    }
	    
	  private:
	    std::shared_ptr< LinkedInternalNode > mPtr;
	};

	//! \return assignable and copyable node for emulating a reference
	// class Node : public ConstNode
	// {
	//   public:
	//     Node() {}
	    
	//     Node( const std::shared_ptr< LinkedInternalNode >& nodePtr ) : ConstNode( nodePtr ) {}

	//     Node( const Node& orig ) : ConstNode( orig ) {}

	//     Node& operator =( const Node& orig )
	//     {
	// 	ConstNode::operator =( orig );
	// 	return *this;
	//     }
	// };

	// class std::iterator_traits< ConstChildrenIterator >
	// {
	//     typedef typename CListT::const_iterator WrapT;
	//     typedef typename std::iterator_traits< WrapT >::difference_type difference_type;
	//     typedef typename std::iterator_traits< WrapT >::iterator_category iterator_category;
	//     typedef ConstNode value_type;
	//     typedef ConstNode& reference;
	//     typedef ConstNode* pointer;
	// };

	// class ConstChildrenIterator
	// {
	//   public:
	//     typedef typename CListT::const_iterator WrapT;
	//     typedef typename std::iterator_traits< WrapT >::difference_type difference_type;
	//     // for generalization beyond random access
	//     typedef typename std::iterator_traits< WrapT >::iterator_category iterator_category;
	//     typedef ConstNode value_type;
	//     typedef const ConstNode& reference;
	//     typedef const ConstNode* pointer;

	//     ConstChildrenIterator() : mWrapped() {}

	//     ConstChildrenIterator( typename CListT::const_iterator wrap ) : mWrapped( wrap ) {}
	    
	//     ConstChildrenIterator( const ConstChildrenIterator& orig ) = default;

	//     ConstChildrenIterator& operator =( const ConstChildrenIterator& orig ) = default;

	//     ConstChildrenIterator& operator ++() { ++mWrapped; return *this; }

	//     ConstChildrenIterator& operator --() { --mWrapped; return *this; }

	//     ConstChildrenIterator& operator +=( const long int& add ) { mWrapped += add; return *this; }

	//     ConstChildrenIterator& operator -=( const long int& add ) { mWrapped -= add; return *this; }

	//     difference_type operator -( const ConstChildrenIterator& other ) const { return mWrapped - other.mWrapped; }

	//     reference operator *() const { return *mWrapped; }

	//     bool operator ==( const ConstChildrenIterator& other ) const { return this->mWrapped == other.mWrapped; }

	//     bool operator !=( const ConstChildrenIterator& other ) const { return this->mWrapped != other.mWrapped; }
		
	//     void swap( ConstChildrenIterator other ) { swap( this->mWrapped, other.mWrapped ); }
	    
	//   private:
	//     WrapT mWrapped;
	// };

	//! \brief stores local structure of tree
	struct LinkedInternalNode
	{
	    LinkedInternalNode( const NodeT& pp, const ValueT& value ) : mParent( pp ), mValue( value )
	    {}

	    LinkedInternalNode( const LinkedInternalNode& n )
		: mChildren( n.mChildren.begin(), n.mChildren.end() )
		, mParent( n.mParent )
		, mValue( n.mValue ) {}

	    // LinkedInternalNode operator =( const LinkedInternalNode& src )
	    // 	: mChildren( src.mChildren )
	    // 	, mParent( src.mParent )
	    // 	, mValue( src.mValue ) {}

	    bool operator ==( const LinkedInternalNode& other )
	    {
		return this->mValue == other.mValue;
	    }

	    CListT mChildren;
	    NodeT mParent;
	    T mValue;
	};

	//! exception class for linked n-ary tree
	class LinkedFixedBranchTreeException : public std::runtime_error
	{
	  public:
	    LinkedFixedBranchTreeException( std::string what ) : std::runtime_error( what ) {}
	};

	LinkedFixedBranchTree( const T& rootValue )
	    : mRoot( std::shared_ptr< LinkedInternalNode > ( new LinkedInternalNode( NodeT( std::shared_ptr< LinkedInternalNode >() ), rootValue ) ) )
	    , mSize( 1 )
	    , mHeight( 1 ) {}

	~LinkedFixedBranchTree()
	{
	    delChildren( mRoot );
	}

	const NodeT& root() const { return mRoot; }

	/*! \return number of nodes in tree */
	size_t size() const
	{
	    return mSize;
	}

	/*! \return length of longest path - 1 */
	size_t height() const
	{
	    return mHeight;
	}

	// accessors to node
	ValueT& value( const NodeT& n ) { return n.mPtr->mValue; }

	const NodeT& parent( const NodeT& n ) const { return n.mPtr->mParent; }

	std::pair< CIterT, CIterT > children( const NodeT& n ) const
	{
	    if( isLeaf( n ) )
		return std::make_pair( n.mPtr->mChildren.end(), n.mPtr->mChildren.end() );
	    return std::make_pair( n.mPtr->mChildren.begin(), n.mPtr->mChildren.end() );
	}

	bool isRoot( const NodeT& n ) const { return !parent( n ).mPtr; }

	bool isLeaf( const NodeT& n ) const
	{
	    if( !n.mPtr )
		throw std::runtime_error( "bad node" );
	    if( n.mPtr->mChildren[ 0 ].mPtr && std::any_of( n.mPtr->mChildren.begin(), n.mPtr->mChildren.end()
							    , [] (const NodeT& n ) { return !n.mPtr; } ) )
		throw std::logic_error( "first child good, another bad" );
	    
	    return !(n.mPtr->mChildren[ 0 ].mPtr);
	}
	
	void expand( const NodeT& n, const std::array< T, NC >& vals )
	{
	    if( !isLeaf( n ) )
		throw LinkedFixedBranchTreeException( "is interior node, cannot expand" );
	    for( uint i = 0; i < NC; ++i )
		n.mPtr->mChildren[ i ].mPtr.reset( new LinkedInternalNode( n.mPtr, vals[ i ] ) );
	    mSize += NC;
	    mHeight = std::max( mHeight, depth( *this, n ) + 2 ); // depth incremented and +1 to count bottom
	}

	void delChildren( const NodeT& n )
	{
	    size_t subSize = subtreeSize( *this, n );
	    for( NodeT& pn : n.mPtr->mChildren )
	    {
		if( pn.mPtr )
		{
		    pn.mPtr->mParent.mPtr.reset();
		    delChildren( pn );
		    pn.mPtr.reset();
		}
	    }
	    mSize -= subSize - 1;
	    mHeight = subtreeHeight( *this, NodeT( mRoot ) );
	}
      private:
	NodeT mRoot;
	size_t mSize, mHeight;
    };

    TREE_TEMPLATE
    typename LTREE::ValueT& value( LTREE& lt, const typename LTREE::NodeT& n )
    {
	return lt.value( n );
    }

    TREE_TEMPLATE
    typename LTREE::NodeT root( const LTREE& lt )
    {
	return lt.root();
	// return typename LTREE::NodeT( lt.root() );
    }

    TREE_TEMPLATE
    typename LTREE::NodeT parent( const LTREE& lt, const typename LTREE::NodeT & n )
    {
	return lt.parent( n );
    }

    TREE_TEMPLATE
    std::pair< typename LTREE::CIterT
    	       , typename LTREE::CIterT >
    children( const LTREE& lt, const typename LTREE::NodeT& n )
    {
	return lt.children( n );
    }

    //! \note add to interface
    //! replace by visitor code
    //! \return the depth of the node i.e. length of path until root - 1
    TREE_TEMPLATE
    size_t depth( const LTREE& lt, const typename LTREE::NodeT& n )
    {
    	typename LTREE::NodeT x = n;
    	size_t d = 0;
    	while( !isRoot( lt, x ) )
    	{
    	    ++d;
    	    x = parent( lt, x );
    	}
    	return d;
    }

    // replace by visitor code
    //! \note add to interface
    //! \return the height of the subtree in t at n
    TREE_TEMPLATE
    size_t subtreeHeight( const LTREE& lt, const typename LTREE::NodeT& n )
    {
    	size_t h = 1;
	typename FixedBranchTreeTraits< LTREE >::CRangeT crange = children( lt, n );
    	for( ; crange.first != crange.second; ++crange.first )
    	    h = std::max( h, 1 + subtreeHeight( lt, *crange.first ) );
    	return h;
    }

    //! \note add to interface
    //! replace by visitor code
    TREE_TEMPLATE
    size_t subtreeSize( const LTREE& lt, const typename LTREE::NodeT& n )
    {
    	if( isLeaf( lt, n ) )
    	    return 1;

    	size_t size = 1;
    	typename FixedBranchTreeTraits< LTREE >::CRangeT crange = children( lt, n );
    	for( ; crange.first != crange.second; ++crange.first )
    	    size += subtreeSize( lt, *crange.first );
    
    	return size;
    }

    TREE_TEMPLATE
    bool isRoot( const LTREE& lt, const typename LTREE::NodeT& n )
    {
	return lt.isRoot( n );
    }

    TREE_TEMPLATE
    bool isLeaf( const LTREE& lt, const typename LTREE::NodeT& n )
    {
	// every node needs to have exactly n children
	return lt.isLeaf( n );
    }

    TREE_TEMPLATE
    void delChildren( LTREE& lt, const typename LTREE::NodeT& n )
    {
	lt.delChildren( n );
    }

    TREE_TEMPLATE
    void expand( LTREE& t, const typename LTREE::NodeT& n, const std::array< T, NC >& vals = std::array< T, NC >() )
    {
	t.expand( n, vals );
    }

    // template< typename T, size_t NC
    // 	      , typename CharT, typename Traits
    // 	      , typename P>
    // std::basic_ostream< CharT, Traits >& print( st::basic_ostream< CharT, Traits >& os
    // 						, const LTREE& t, const typename LTREE::NodeT& n
    // 						, std::function< P( const typename LTREE::ValueT& ) > converter = [] ( const typename LinkedFixedBranchTree< T, NC, >::ValueT& v ) {return v;} )
    // {
	
    // }
}

#endif
