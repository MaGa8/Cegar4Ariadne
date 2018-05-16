#ifndef GRAPH_VALUE_HPP
#define GRAPH_VALUE_HPP

//! \class abstract base class for value to store in graph
struct IGraphValue
{
    virtual bool isInside() const = 0;
};

//! \class value to store in graph of region inside first initial abstraction
template< typename TreeNodeT >
class InsideGraphValue : public IGraphValue
{
  public:
    InsideGraphValue( TreeNodeT& treeNode )
	: mTreeNode( treeNode )
    {}

    InsideGraphValue( const InsideGraphValue& orig ) = default;

    InsideGraphValue& operator =( const InsideGraphValue& orig ) = default;

    const TreeNodeT& treeNode() const { return mTreeNode; }

    bool isInside() const { return true; }
  private:
    TreeNodeT mTreeNode;
};

//! \class value to store in graph of region outside first initial abstractio
struct OutsideGraphValue : public IGraphValue
{
    bool isInside() const { return false; }
};

#endif
