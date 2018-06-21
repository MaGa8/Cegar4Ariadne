#ifndef GRAPH_VALUE_HPP
#define GRAPH_VALUE_HPP

#include "geometry/box.hpp"
#include "numeric/logical.hpp"

//! \class abstract base class for value to store in graph
struct IGraphValue
{
    virtual ~IGraphValue() = default;
    
    virtual bool isInside() const = 0;
};

//! \class value to store in graph of region inside first initial abstraction
template< typename EnclosureT >
class InsideGraphValue : public IGraphValue
{
  public:
    InsideGraphValue()
	: mId( -1 )
	, mEnclosure()
	, mSafe( false )
	, mTransSafe( Ariadne::indeterminate )
    {}
    
    InsideGraphValue( const unsigned long& id, const EnclosureT& e, const Ariadne::ValidatedKleenean& safe )
	: mId( id )
	, mEnclosure( e )
	, mSafe( safe )
	, mTransSafe( Ariadne::indeterminate )
    {}

    InsideGraphValue( const InsideGraphValue& orig ) = default;

    virtual ~InsideGraphValue() = default;

    InsideGraphValue& operator =( const InsideGraphValue& orig ) = default;

    bool isInside() const { return true; }

    //! \return unique identifier of this value
    const unsigned long& id() const { return mId; }

    //! \return box stored
    const EnclosureT& getEnclosure() const { return mEnclosure; }

    //! \return true if the box stored is completely covered by all constraints, indeterminate if it lies within the initial abstraction but is not covered completely by all constraints and false if it lies outside of the initial abstraction
    Ariadne::ValidatedKleenean isSafe() const { return mSafe; }

    //! \return true if no unsafe state can be reached, false if unsafe state can be reached, indeterminate if not set
    Ariadne::ValidatedKleenean isTransSafe() const { return mTransSafe; }

    bool operator ==( const InsideGraphValue< EnclosureT >& tv ) const
    {
	return this->mId == tv.mId;
    }

    //! \brief initializes graph value, required for pooling
    void init( const unsigned long& id, const EnclosureT& e, const Ariadne::ValidatedKleenean& safe )
    {
	this->mId = id;
	this->mEnclosure = e;
	this->mSafe = safe;
	this->mTransSafe = Ariadne::indeterminate;
    }

    void resetTransSafe()
    {
	mTransSafe = Ariadne::indeterminate;
    }

    void transSafe()
    {
	// if( definitely( mTransSafe ) )
	    // throw std::logic_error( "attempt to set node safe, even though node is already set safe" );
	mTransSafe = true;
    }

    void transUnsafe()
    {
	if( definitely( mTransSafe ) )
	    throw std::logic_error( "attempt to set node unsafe, even though node is already set safe" );
	mTransSafe = false;
    }

  private:
    unsigned long mId;
    EnclosureT mEnclosure;
    Ariadne::ValidatedKleenean mSafe;
    Ariadne::ValidatedKleenean mTransSafe;
};

//! \class value to store in graph of region outside first initial abstractio
struct OutsideGraphValue : public IGraphValue
{
    virtual ~OutsideGraphValue() = default;

    bool isInside() const { return false; }
};

template< typename E, typename CharT, typename TraitsT >
std::basic_ostream< CharT, TraitsT >& operator <<( std::basic_ostream< CharT, TraitsT >& os, const InsideGraphValue< E >& val )
{
    os << val.id() << ": ";
    os << val.getEnclosure() << " (";
    os << val.isSafe() << ", ";
    os << val.isTransSafe() << ") ";
    return os;
    // return os << val.id() << ": " << val.getEnclosure() << " (" << val.isSafe() << ", " << val.isTransSafe() << ") ";
}

template< typename CharT, typename TraitsT >
std::basic_ostream< CharT, TraitsT >& operator <<( std::basic_ostream< CharT, TraitsT >& os, const OutsideGraphValue& val )
{
    return os << "[outside] ";
}

#endif
