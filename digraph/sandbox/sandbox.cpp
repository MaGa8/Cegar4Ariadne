#include "adjacencyDiGraph.hpp"

#include <iostream>

using namespace graph;

typedef AdjacencyDiGraph< int, VecMap, InVec, InVec > G;

int main()
{
    G g;
    const uint val1 = 0, val2 = 1;
    addVertex( g, val1 );
    addVertex( g, val2 );

    typename G::VIterT iv = vertices( g ).first;
    std::cout << iv->first << std::endl;
    
    iv = findVertex( g, val2 );
    std::cout << iv->first << std::endl;

    return 0;
}
