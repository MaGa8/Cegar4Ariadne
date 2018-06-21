#include "geometry/geometry.hpp"
#include "geometry/function_set.hpp"
#include "geometry/box.hpp"

// compiled by
// g++ -I ../../ariadne/ -I ../../ariadne/source/ -L ../../ariadne/build/ -l ariadne-kernel -l ariadne -w -o inside inside.cpp

int main()
{
    Ariadne::BoundedConstraintSet *bs1, *bs2;

    inside( *bs1, *bs2, 0.00001 );
    return 0;
}
