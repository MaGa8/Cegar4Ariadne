#ifndef GRAPH_VALUE_HPP
#define GRAPH_VALUE_HPP

#include "geometry/box.hpp"
#include "numeric/logical.hpp"

//! \class abstract base class for value to store in graph
struct IGraphValue
{
    virtual bool isInside() const = 0;
};

//! \class value to store in graph of region inside first initial abstraction
template< typename EnclosureT >
class InsideGraphValue : public IGraphValue
{
  public:
    InsideGraphValue( const unsigned long& id, const EnclosureT& e, const Ariadne::ValidatedKleenean& safe )
	: mId( id )
	, mEnclosure( e )
	, mSafe( safe )
    {}

    InsideGraphValue( const InsideGraphValue& orig ) = default;

    InsideGraphValue& operator =( const InsideGraphValue& orig ) = default;

    bool isInside() const { return true; }

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

  private:
    unsigned long mId;
    EnclosureT mEnclosure;
    Ariadne::ValidatedKleenean mSafe;
};

//! \class value to store in graph of region outside first initial abstractio
struct OutsideGraphValue : public IGraphValue
{
    bool isInside() const { return false; }
};

#endif
