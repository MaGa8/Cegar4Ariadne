#include "testGroupInterface.hpp"

#include <iomanip>
#include <chrono>

ITestGroup::ITestGroup( std::string name, uint size, uint repetitions, uint level, std::ostream& out )
    : ITest( name, size, repetitions )
    , mLevel( level )
    , mOut( out )
    , mResult( true )
{}

void ITestGroup::iterate()
{
    for( TestPolicyPair& tp : mTests )
    {
	std::chrono::time_point< std::chrono::high_resolution_clock > tstart = std::chrono::high_resolution_clock::now();
	uint noSpacesInsert = RESULT_PRINT_POS - mLevel*LEVEL_INDENT - tp.first->mName.size();
	bool result = tp.second->run( tp.first.get() );
	mResult = mResult && result;
	std::chrono::milliseconds testDuration = std::chrono::duration_cast< std::chrono::milliseconds >( std::chrono::high_resolution_clock::now() - tstart );
	mOut << std::setw( mLevel*LEVEL_INDENT )
	     << tp.first->mName << " :"
	     << std::setw( noSpacesInsert )
	     << (result ? "success" : "failure")
	     << std::setw( 4 + LEVEL_INDENT ) << testDuration.count() << " ms"
	     << std::endl;
    }
}

bool ITestGroup::check() const
{
    return mResult;
}
	
void ITestGroup::addTest( ITest* pTest, std::shared_ptr< ITestRunner > pRunner )
{
    mTests.push_back( std::make_pair( std::unique_ptr< ITest >( pTest ), pRunner ) );
}
