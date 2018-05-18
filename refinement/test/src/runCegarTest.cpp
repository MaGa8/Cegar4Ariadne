#include "cegarTest.hpp"

int main()
{
    std::shared_ptr< OnlyOnceRunner > oor( new OnlyOnceRunner() );
    oor->run( new CegarTest( 50, 100 ) );
    return 0;
}
