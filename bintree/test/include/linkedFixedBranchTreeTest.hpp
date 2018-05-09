#ifndef LINKED_FIXED_BRANCH_TREE_TEST_HPP
#define LINKED_FIXED_BRANCH_TREE_TEST_HPP

#include "linkedFixedBranchTree.hpp"
#include "testInterface.hpp"
#include "testGroupInterface.hpp"

#include <iterator>
#include <algorithm>
#include <random>
#include <memory>
#include <utility>

using namespace tree;

class LinkedFixedBranchTreeTest : public ITestGroup
{
  public:
    static const uint mNoChildren = 10;
    typedef LinkedFixedBranchTree< int, 10 > TestTreeT;
    static std::default_random_engine mRandom;
    static std::uniform_int_distribution<> mIntDist, mChildDist;
    static std::uniform_real_distribution<> mProbDist;

    template< typename DistT, typename TreeT >
    static typename TreeT::NodeT randomExpandTree( TreeT& tree, DistT valueDist )
    {
    	std::uniform_int_distribution< int > childrenD( 0, TreeT::N - 1 );

	typename TreeT::NodeT expandNode = root( tree );
	while( !isLeaf( tree, expandNode ) )
	{
	    // D( std::cout << "current node " << value( tree, expandNode ) << ", obtain children" <<  std::endl; );
	    typename FixedBranchTreeTraits< TreeT >::CRangeT crange = children( tree, expandNode );
	    // D( std::cout << "go to child to expand" << std::endl; );
	    std::advance( crange.first, childrenD( mRandom ) );
	    // D( std::cout << "deref expand node" << std::endl; );
	    expandNode = *crange.first;
	}
	expand( tree, expandNode );
	return expandNode;
    }

    // test leaf
    // generate random tree, traverse down from root to some leaf and test
    class LeafTest : public ITest
    {
	STATEFUL_TEST( LeafTest );
      private:
	std::unique_ptr< TestTreeT > mpTree;
    };

    // expand tree and test number of nodes
    class ExpandSizeTest : public ITest
    {
	STATEFUL_TEST( ExpandSizeTest );
      private:
	std::unique_ptr< TestTreeT > mpTree;
	uint mNoExpansions;
    };

    // expand tree and test height
    class ExpandHeightTest : public ITest
    {
	STATEFUL_TEST( ExpandHeightTest );
      private:
	std::unique_ptr< TestTreeT > mpTree;
	uint mNoExpansions;
    };
    
    // expand tree and test that root is root
    class RootTest : public ITest
    {
	STATEFUL_TEST( RootTest );
      private:
	std::unique_ptr< TestTreeT > mpTree;
	typename TestTreeT::NodeT mRoot;
    };

    // expand tree and test that no descendent of root is root
    class NotRootTest : public ITest
    {
	STATEFUL_TEST( NotRootTest );
      private:
	std::unique_ptr< TestTreeT > mpTree;
    };

    // set up tree by some initial expansions; iterate: first expand then delete children; check: ensure node is now leaf
    class DeleteTest : public ITest
    {
	STATEFUL_TEST( DeleteTest );
      private:
	std::unique_ptr< TestTreeT > mpTree;
	std::optional< typename TestTreeT::NodeT > mDeletedAt;
	const double mTraverseThresh = 1.0 / std::log2( mTestSize );
    };

    // test size following random expansion and delete
    class DeleteSizeTest : public ITest
    {
	STATEFUL_TEST( DeleteSizeTest );
      private:
	std::unique_ptr< TestTreeT > mpTree;
	const double mTraverseThresh = 1.0 / std::log2( mTestSize );
    };

    class DeleteHeightTest : public ITest
    {
	STATEFUL_TEST( DeleteHeightTest );
      private:
	std::unique_ptr< TestTreeT > mpTree;
	const double mTraverseThresh = 1.0 / std::log2( mTestSize );
	uint mHeightBefore;
    };

    // test delete height
    // generate tree of b=m with n expansions, then delete root
    
    GROUP_CTOR_DECL( LinkedFixedBranchTreeTest );

    void init();
    
}; // test group

#endif
