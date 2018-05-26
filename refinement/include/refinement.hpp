#ifndef REFINER_HPP
#define REFINER_HPP

#include "geometry/box.hpp"

/*
  Refinement classes represent a procedure for subdividing enclosures of one type into smaller enclosures of the same type.
  They should provide
  template< typename E > operator ()( const E& enclosure )
 */

//! \class refine along the coordinate of the largest interval
struct LargestSideRefiner
{

    template< typename E >
    std::vector< E > operator ()( const E& b ) const;

    template< typename IntervalT >
    std::vector< Ariadne::Box< IntervalT > > operator ()( const Ariadne::Box< IntervalT >& b ) const
    {
	Ariadne::Pair< Ariadne::Box< IntervalT >, Ariadne::Box< IntervalT > > splits = b.split();
	return std::vector< Ariadne::Box< IntervalT > >( {splits.first, splits.second} );
    }
};

#endif
