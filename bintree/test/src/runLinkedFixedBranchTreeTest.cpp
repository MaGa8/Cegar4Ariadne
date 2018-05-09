#include "linkedFixedBranchTreeTest.hpp"

int main()
{
    OnlyOnceRunner oor;
    LinkedFixedBranchTreeTest tst( 100, 1000 );
    oor.run( &tst );
    return 0;
}
