#include "refinementTreeTest.hpp"

int main()
{
    std::shared_ptr< OnlyOnceRunner > oor( new OnlyOnceRunner() );
    oor->run( new RefinementTreeTest( 100, 1000 ) );
    return 0;
}
