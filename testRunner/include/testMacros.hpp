#ifndef TEST_MACROS_HPP
#define TEST_MACROS_HPP

#define D(x) do{ if( DEBUG ) x } while( false )

#define Dout(x) D( std::cout << x << std::endl; )

#endif
