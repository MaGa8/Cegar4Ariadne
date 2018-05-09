#include "allTests.hpp"

GROUP_CTOR( AllTests, "all" )

void AllTests::init()
{
    std::shared_ptr< OnlyOnceRunner > poor( new OnlyOnceRunner() );
    addTest( new LinkedFixedBranchTreeTest( mTestSize, mRepetitions ), poor );
    addTest( new AdjacencyDiGraphTest( mTestSize, mRepetitions ), poor );
}
