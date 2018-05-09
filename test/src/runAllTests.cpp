#include "allTests.hpp"

int main()
{
    OnlyOnceRunner poor;
    AllTests allTests( 25, 1000 );
    poor.run( &allTests );
}
