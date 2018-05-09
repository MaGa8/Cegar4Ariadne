#include "adjacencyDiGraphTest.hpp"

int main()
{
    AdjacencyDiGraphTest t( 50, 1000 );
    OnlyOnceRunner oor;
    oor.run( &t );
    
    return 0;
}
