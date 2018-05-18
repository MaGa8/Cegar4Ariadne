#include "testRunnerInterface.hpp"

bool OnlyOnceRunner::run( ITest* pTest ) const
{
    pTest->init();
    pTest->iterate();
    return pTest->check();
}

bool StatelessRunner::run( ITest* pTest ) const
{
    for( uint cRep = 0; cRep < pTest->mRepetitions; ++cRep )
    {
	pTest->iterate();
	if( !pTest->check() )
	    return false;
    }
    return true;
}

bool ContinuousRunner::run( ITest* pTest ) const
{
    for( uint cRep = 0; cRep < pTest->mRepetitions; ++cRep )
    {
	pTest->init();
	for( uint cSize = 0; cSize < pTest->mTestSize; ++cSize )
	    pTest->iterate();
	if( !pTest->check() )
	    return false;
    }
    return true;
}

std::default_random_engine ContinuousRandomRunner::mRandom = std::default_random_engine( std::random_device()() );

bool ContinuousRandomRunner::run( ITest* pTest ) const
{
    std::uniform_int_distribution<> sizeDist( 0, pTest->mTestSize );
    for( uint cRep = 0; cRep < pTest->mRepetitions; ++cRep )
    {
	uint testSize = sizeDist( mRandom );
	pTest->init();
	for( uint cSize = 0; cSize < testSize; ++cSize )
	    pTest->iterate();
	if( !pTest->check() )
	    return false;
    }
    return true;
}

bool InterleaveRunner::run( ITest* pTest ) const
{
    for( uint cRep = 0; cRep < pTest->mRepetitions; ++cRep )
    {
	pTest->init();
	for( uint cSize = 0; cSize < pTest->mTestSize; ++cSize )
	{
	    pTest->iterate();
	    if( !pTest->check() )
		return false;
	}
    }
    return true;
}

std::default_random_engine InterleaveRandomRunner::mRandom = std::default_random_engine( std::random_device()() );

bool InterleaveRandomRunner::run( ITest* pTest ) const
{
    std::uniform_int_distribution<> sizeDist( 0, pTest->mTestSize );
    for( uint cRep = 0; cRep < pTest->mRepetitions; ++cRep )
    {
	uint testSize = sizeDist( mRandom );
	pTest->init();
	for( uint cSize = 0; cSize < testSize; ++cSize )
	{
	    pTest->iterate();
	    if( !pTest->check() )
		return false;
	}
    }
    return true;
}
