#ifndef REFINER_HPP
#define REFINER_HPP

#include "geometry/box.hpp"

//! \interface for refinement strategies subdividing a box into smaller boxes
template< typename IntervalT >
struct IRefinementStrategy
{
    virtual std::vector< Ariadne::Box< IntervalT > > refine( const Ariadne::Box< IntervalT >& b ) const = 0;
};

//! \class refine along the coordinate of the largest interval
template< typename IntervalT >
struct LargestSideRefiner : public IRefinementStrategy< IntervalT >
{
    std::vector< Ariadne::Box< IntervalT > > refine( const Ariadne::Box< IntervalT >& b ) const
    {
	Ariadne::Pair< Ariadne::Box< IntervalT >, Ariadne::Box< IntervalT > > splits = b.split();
	return std::vector< Ariadne::Box< IntervalT > >( {splits.first, splits.second} );
    }
};

#endif
