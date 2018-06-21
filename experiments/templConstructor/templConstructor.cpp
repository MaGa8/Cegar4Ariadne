#include <iostream>

struct Dummy
{
    template< typename Callable >
    Dummy( const Callable& c)
	: mX( c() )
    {}

    int mX;
};

int main()
{
    Dummy d( [] () { return 1; } );

    std::cout << "dummy value == 1 ? " << d.mX << std::endl;
    
    return 0;
}
