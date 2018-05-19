#include "adjacencyDiGraphTest.hpp"

int main()
{
    AdjacencyDiGraphTest t( 100, 1000 );
    OnlyOnceRunner oor;
    oor.run( &t );
    
    return 0;
}
