#ifndef TREE_INTERFACE_HPP
#define TREE_INTERFACE_HPP

#include <array>
#include <iostream>

namespace tree
{
    blablbabla this is a bug ohohohohnonononoo
    /*!
      \brief specifies types to be defined by a tree of fixed number of branches
      \note non-const types should be implicitly convertible to corresponding const types
     */
    template< typename NTreeT >
    struct FixedBranchTreeTraits
    {
	typedef typename NTreeT::ValueT ValueT;
	// typedef typename NTreeT::ConstValueT ConstValueT;
	
	typedef typename NTreeT::NodeT NodeT;
	// typedef typename NTreeT::ConstNodeT ConstNodeT;
	
	typedef typename NTreeT::CListT CListT;
	typedef typename NTreeT::CIterT CIterT;
	// typedef typename NTreeT::ConstCIterT ConstCIterT;

	static const size_t N = NTreeT::N;

	typedef std::pair< CIterT, CIterT > CRangeT;
	// typedef std::pair< ConstCIterT, ConstCIterT > ConstCRangeT;
    };

    /*! 
      \interface for trees whose nodes are either leafs or have exactly N children at any time
      \param T type of value stored at nodes
      \param N number of children
      \param NodeT type of node
      \param ChildIt type of iterator over chilren of a node
      !*/
    template< typename T, size_t N >
    struct FixedBranchTreeInterface
    {
	//! \return number of nodes in tree
	virtual size_t size() const = 0;

	//! \return length of longest path in tree - 1
	virtual size_t height() const = 0;    
    };

    //! \return root node of tree
    template< typename NTreeT >
    typename FixedBranchTreeTraits< NTreeT >::NodeT root( const NTreeT& t );

    // //! \return const root node of tree t
    // template< typename NTreeT >
    // const typename FixedBranchTreeTraits< NTreeT >::ConstNodeT root( const NTreeT& t )
    // {
    // 	// avoid copying step
    // 	typename NTreeT::NodeT nonConst = root( const_cast< NTreeT& >( t ) );
    // 	return static_cast< typename NTreeT::ConstNodeT& >( nonConst );
    // }

    //! \return parent node in tree t of node n
    template< typename NTreeT >
    typename FixedBranchTreeTraits< NTreeT >::NodeT parent( const NTreeT& t, const typename FixedBranchTreeTraits< NTreeT >::NodeT& n );

    //! \return const parent node in tree t of node n
    // cannot be implemented by cast, because casting ConstNodeT to NodeT may not be feasible
    // FIXED! NodeT now derived from ConstNodeT so can be downcasted
    // template< typename NTreeT >
    // typename FixedBranchTreeTraits< NTreeT >::ConstNodeT parent( const NTreeT& t, const typename FixedBranchTreeTraits< NTreeT >::ConstNodeT& n )
    // {
    // 	typename NTreeT::NodeT nonConst = parent( const_cast< NTreeT& >( t )
    // 						  , static_cast< typename NTreeT::NodeT& > ( const_cast< typename NTreeT::ConstNodeT& >( n ) ) );
    // 	return static_cast< typename NTreeT::ConstNodeT& >( nonConst );
    // }

    //! \return value  of node n in tree t
    template< typename NTreeT >
    typename FixedBranchTreeTraits< NTreeT >::ValueT& value( NTreeT& t, const typename FixedBranchTreeTraits< NTreeT >::NodeT& n );

    //! \return value  of node n in tree t
    template< typename NTreeT >
    const typename FixedBranchTreeTraits< NTreeT >::ValueT& value( const NTreeT& t, const typename FixedBranchTreeTraits< NTreeT >::NodeT& n )
    {
	typename NTreeT::ValueT& nonConst = value( const_cast< NTreeT& >( t ), n );
	return nonConst;
    }

    //! \return pair of iterators (begin/end) over children of node n in tree t which should be an empty range if n has no children
    template< typename NTreeT >
    typename FixedBranchTreeTraits< NTreeT >::CRangeT children( const NTreeT& t, const typename FixedBranchTreeTraits< NTreeT >::NodeT& n );

    //! \return pair of iterators to const (begin/end) over children of node n in tree t which should be an empty range if n has no children
    // template< typename NTreeT >
    // typename FixedBranchTreeTraits< NTreeT >::ConstCRangeT children( const NTreeT& t
    // 								     , const typename FixedBranchTreeTraits< NTreeT >::ConstNodeT& n )
    // {
    // 	NTreeT& nonConstTree = const_cast< NTreeT& >( t );
    // 	typename NTreeT::NodeT& nonConstNode = static_cast< typename NTreeT::NodeT& >( const_cast< typename NTreeT::ConstNode& >( n ) );
    // 	typename FixedBranchTreeTraits< NTreeT >::CRangeT nonConst = children( nonConstTree, nonConstNode );
    // 	return static_cast< typename FixedBranchTreeTraits< NTreeT >::ConstCRangeT >( nonConst );
    // }

    //! \return true if n is the root node in t
    template< typename NTreeT >
    bool isRoot( const NTreeT& t, const typename FixedBranchTreeTraits< NTreeT >::NodeT& n );

    //! \return true if n is a leaf in t
    template< typename NTreeT >
    bool isLeaf( const NTreeT& t, const typename FixedBranchTreeTraits< NTreeT >::NodeT& n );

    //! \brief delete subtree rooted at n, not including n
    template< typename NTreeT >
    void delChildren( NTreeT& t, const typename FixedBranchTreeTraits< NTreeT >::NodeT& n );

    //! \brief expand node n and initialize the children to value vals
    //! \note vals should be optional and should use the default constructor if not given
    template< typename NTreeT >
    void expand( NTreeT& t
		 , const typename FixedBranchTreeTraits< NTreeT >::NodeT& n
		 , const typename FixedBranchTreeTraits< NTreeT >::CListT& vals );
}

#endif 
