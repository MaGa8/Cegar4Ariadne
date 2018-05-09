
#include <cmath.h>
#include <vector>


/*!
  \brief perfect binary tree, sequentially stored
  \param T type to store as nodes
 */
template< class T >
class SeqBinTree
{
  public:
    /*!
      \brief Tree iterator allowing for branching. To be used for accessing and modifying the tree.
    */
    class Iterator
    {
	Iterator( size_t index ) : mIndex( index ) {}

	isRoot() { return mIndex == 0; }

	// perfect bin tree contains one less interior than leaf node
	isLeaf() { return mIndex < pow( 2, depth() - 1 ); }

      private:
	size_t mIndex;
    }; // iterator

    BinTree() {}

    size_t size() const
    {
	return mNodes.size();
    }

    size_t depth() const
    {
	return std::ceil( std::log2( this->size() + 1 ) - 1 );
    }

  private:
    std::vector< T > mNodes;
}; // BinTree
