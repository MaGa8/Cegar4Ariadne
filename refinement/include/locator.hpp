#ifndef LOCATOR_HPP
#define LOCATOR_HPP

#include "refinementTree.hpp"

#include <vector>
#include <functional>
#include <algorithm>

/*!
  \interface for locator function objects which shall be used to determine the state to refine
  \param EnclosureT type of enclosure to refine
*/
// template< typename E >
// struct ILocator
// {
//     typedef std::vector< std::reference_wrapper< typename RefinementTree< E >::NodeT > > NodeRefVec;

//     template< typename IterT >
//     virtual NodeRefVec operator ()( IterT& i1, const IterT& i2 ) const = 0;
// };

//! \class represents policy of refining all nodes in a counterexample
template< typename E >
struct CompleteCounterexample
{
    typedef std::vector< std::reference_wrapper< const typename RefinementTree< E >::NodeT > > NodeRefVec;

    template< typename IterT >
    NodeRefVec operator ()( IterT i1, const IterT& i2 ) const
    {
	NodeRefVec ns;
	ns.reserve( std::distance( i1, i2 ) );
	std::transform( i1, i2, std::back_inserter( ns )
			, [] (const typename RefinementTree< E >::NodeT& n) {return std::ref( n );} );
	return ns;
    }
};




#endif
