#ifndef TREE_VALUE_HPP
#define TREE_VALUE_HPP

#include "geometry/box.hpp"
#include "numeric/logical.hpp"

//! \class value to store in refinement tree
//! \param EnclosureT type used to represent abstract states
//! \note needs to be specialized for each enclosure type, as default construction is case dependent
//  made a base class even though it can be instantiated so that methods don't have to be duplicated
template< typename EnclosureT >
class AbstractInteriorTreeValue
{
  public:
    AbstractInteriorTreeValue( const unsigned long& id, const EnclosureT& e, const Ariadne::ValidatedKleenean& safe )
	: mId( id )
	, mEnclosure( e )
	, mSafe( safe )
    {}

    AbstractInteriorTreeValue( const AbstractInteriorTreeValue< EnclosureT >& orig ) = default;

    AbstractInteriorTreeValue& operator =( const AbstractInteriorTreeValue< EnclosureT >& orig ) = default;

    //! \return unique identifier of this value
    const unsigned long& id() const { return mId; }

    //! \return box stored
    const EnclosureT& getEnclosure() const { return mEnclosure; }

    //! \return true if the box stored is completely covered by all constraints, indeterminate if it lies within the initial abstraction but is not covered completely by all constraints and false if it lies outside of the initial abstraction
    Ariadne::ValidatedKleenean isSafe() const { return mSafe; }

    bool operator ==( const AbstractInteriorTreeValue< EnclosureT >& tv ) const
    {
	return this->mId == tv.mId;
    }

    virtual bool isLeaf() const { return false; }
	
  private:
    unsigned long mId;
    EnclosureT mEnclosure;
    Ariadne::ValidatedKleenean mSafe;
};

template< typename EnclosureT >
class InteriorTreeValue : public AbstractInteriorTreeValue< EnclosureT >
{};

//! \class value to store in leaves of refinement tree, does map to graph
template< typename EnclosureT, typename VertexT >
class LeafTreeValue : public InteriorTreeValue< EnclosureT >
{
  public:
    LeafTreeValue() {}

    LeafTreeValue( const unsigned long& id, const EnclosureT& e, const Ariadne::ValidatedKleenean& safe )
	: InteriorTreeValue< EnclosureT >( id, e, safe )
    {}

    LeafTreeValue( const LeafTreeValue< EnclosureT, VertexT >& orig ) = default;

    LeafTreeValue& operator =( const LeafTreeValue< EnclosureT, VertexT >& orig ) = default;

    const VertexT& graphNode() const { return mGraphNode; }
	
    bool isLeaf() const { return true; }

    void setGraphNode( const VertexT& v ) { mGraphNode = v; }
  private:
    VertexT mGraphNode;
};

template<>
class InteriorTreeValue< Ariadne::ExactBoxType > : public AbstractInteriorTreeValue< Ariadne::ExactBoxType >
{
  public:
    InteriorTreeValue()
	: AbstractInteriorTreeValue< Ariadne::ExactBoxType >( -1, Ariadne::ExactBoxType::zero( 0 ), false )
    {}

    InteriorTreeValue( const unsigned long& id, const Ariadne::ExactBoxType& e, const Ariadne::ValidatedKleenean& safe )
	: AbstractInteriorTreeValue( id, e, safe )
    {}

    InteriorTreeValue( const InteriorTreeValue< Ariadne::ExactBoxType >& orig ) = default;

    InteriorTreeValue& operator =( const InteriorTreeValue< Ariadne::ExactBoxType >& orig ) = default;
};

#endif
