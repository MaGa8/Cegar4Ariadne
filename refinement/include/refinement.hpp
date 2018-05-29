#ifndef REFINER_HPP
#define REFINER_HPP

#include "geometry/box.hpp"

#include <random>

/*
  Refinement classes represent a procedure for subdividing enclosures of one type into smaller enclosures of the same type.
  They should provide
  template< typename E > operator ()( const E& enclosure )

  return string description
  constexpr string name() const
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

    std::string name() const { return std::string( "largest_side_refiner" ); }
};

class RandomRefiner
{
  public:
    RandomRefiner( const double& minFrac, const double& maxFrac )
	: mFracDist( minFrac, maxFrac )
	, mRandom( std::random_device()() )
    {}

    std::string name() const {return std::string( "random_position_refiner" ); }

    template< typename E >
    std::vector< E > operator ()( const E& e);

    template< typename IntervalT >
    uint widestDimension( const Ariadne::Box< IntervalT >& bx ) const
    {
	uint nLargest = 0;
	typename IntervalT::WidthType wlargest = bx[ nLargest ].width();
	for( uint i = 1; i < bx.dimension(); ++i )
	{
	    typename IntervalT::WidthType wi = bx[ i ].width();
	    if( definitely( wlargest < wi ) )
	    {
		nLargest = i;
		wlargest = wi;
	    }
	}
	return nLargest;
    }

    std::vector< Ariadne::ExactBoxType > operator ()( const Ariadne::ExactBoxType& b )
    {
	uint widest = widestDimension( b );
	Ariadne::Box< Ariadne::ExactIntervalType > lower = b, upper = b;
	Ariadne::Value< Ariadne::FloatDP > splitPoint = cast_exact( b[ widest ].lower() + mFracDist( mRandom ) * b[ widest ].width() );
	lower[ widest ] = Ariadne::ExactIntervalType( b[ widest ].lower(), splitPoint );
	upper[ widest ] = Ariadne::ExactIntervalType( splitPoint, b[ widest ].upper() );
	return {lower, upper};
    }

  private:
    std::uniform_real_distribution<> mFracDist;
    std::default_random_engine mRandom;
};

#endif
