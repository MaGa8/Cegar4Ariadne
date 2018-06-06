#ifndef TERMINATION_HPP
#define TERMINATION_HPP

#include <chrono>

/*! 
  \brief termination should support
  template< typename Rtree > void start( const Rtree& )
  template< typename Rtree > bool operator ()( const Rtree& )
  to indicate whether cegar should terminate
*/

//! \brief terminate after a given amount of time has elapsed
//! \param D duration type
template< typename D, typename C = std::chrono::high_resolution_clock >
class LimitedTime
{
  public:
    LimitedTime( const D& given )
	: mGiven( given )
	, mStart()
    {}

    template< typename Rtree >
    void start( const Rtree& rtree )
    {
	mStart = C::now();
    }

    template< typename Rtree >
    bool operator ()( const Rtree& rtree ) const
    {
	return (C::now() - mStart ) > mGiven;
    }
    
  private:
    const D mGiven;
    typename C::time_point mStart;
};

class LimitedIterations
{
  public:
    LimitedIterations( const uint& iterations )
	: mIterations( iterations )
	, mCount( 0 )
    {}

    template< typename Rtree >
    void start( const Rtree& rtree )
    {
	mCount = 0;
    }

    template< typename Rtree >
    bool operator ()( const Rtree& rtree )
    {
	++mCount;
	return mCount > mIterations;
    }

  private:
    const uint mIterations;
    uint mCount;
};


#endif
