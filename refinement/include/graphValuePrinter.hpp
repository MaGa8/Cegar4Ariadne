#ifndef GRAPH_VALUE_PRINTER_HPP
#define GRAPH_VALUE_PRINTER_HPP

#include "refinementTree.hpp"

template< typename E >
struct GraphValuePrinter
{

    static std::function< GraphValuePrinter< E >( const typename RefinementTree< E >::MappingT::ValueT& ) >
    makeFun( bool printId, bool printEnc, bool printLoc, bool printTrans )
    {
	return [=] (const typename RefinementTree< E >::MappingT::ValueT& val ) {
		   return GraphValuePrinter( *val, printId, printEnc, printLoc, printTrans ); };
    }
    
    GraphValuePrinter( const IGraphValue& gval, bool printId, bool printEnc, bool printLoc, bool printTrans )
	: mEnc( Ariadne::Vector< typename E::IntervalType >() )
	, mPrintEnc( printEnc ), mPrintLoc( printLoc ), mPrintTrans( printTrans ), mPrintId( printId )
    {
	if( gval.isInside() )
	{
	    const InsideGraphValue< E >& pIn = static_cast< const InsideGraphValue< E > &>( gval );
	    mEnc = pIn.getEnclosure();
	    mLocSafety = pIn.isSafe();
	    mTransSafety = pIn.isTransSafe();
	    mId = pIn.id();
	}
	else
	{
	    mLocSafety = false;
	    mTransSafety = false;
	    mId = 0;
	}
    }
	
    E mEnc;
    Ariadne::ValidatedKleenean mLocSafety, mTransSafety;
    uint mId;
    bool mPrintEnc, mPrintLoc, mPrintTrans, mPrintId;
};

template< typename E, typename CharT, typename TraitsT >
std::basic_ostream< CharT, TraitsT >& operator <<( std::basic_ostream< CharT, TraitsT >& os, const GraphValuePrinter< E >& info )
{
    if( info.mPrintId )
    	os << info.mId << ":  ";
    if( info.mPrintEnc )
    	os << info.mEnc << " ";
    if( info.mPrintLoc )
    	os << "(" << info.mLocSafety << ") ";
    if( info.mPrintId )
    	os << "(" << info.mTransSafety << ")";

    return os;
}

#endif
