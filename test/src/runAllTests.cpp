#include "allTests.hpp"

int main()
{
    OnlyOnceRunner poor;
    AllTests allTests( 100, 1000 );
    poor.run( &allTests );
}
