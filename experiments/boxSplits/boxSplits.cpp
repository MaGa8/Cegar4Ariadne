#include "geometry/box.hpp"
#include "geometry/box.decl.hpp"

#include <iostream>

int main()
{
    Ariadne::ExactBoxType b( { {0,10}, {0,8} } );
    Ariadne::Pair< Ariadne::ExactBoxType, Ariadne::ExactBoxType > splits = b.split();
    std::cout << "split " << b << " along widest coordinate is " << splits.first << "  and  " << splits.second << std::endl;

    return 0;
}
