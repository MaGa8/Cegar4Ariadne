#include "refinementTreeTest.hpp"

int main()
{
    std::shared_ptr< OnlyOnceRunner > oor( new OnlyOnceRunner() );
    oor->run( new RefinementTreeTest( 50, 100 ) );
    return 0;
}
