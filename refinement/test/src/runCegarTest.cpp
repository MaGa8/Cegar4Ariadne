#include "cegarTest.hpp"

int main()
{
    std::shared_ptr< OnlyOnceRunner > oor( new OnlyOnceRunner() );
    oor->run( new CegarTest( 100, 1000 ) );
    return 0;
}
