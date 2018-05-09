#include "testInterface.hpp"

ITest::ITest( std::string name, uint size, uint repetitions )
    : mName( name )
    , mTestSize( size )
    , mRepetitions( repetitions )
{}

// do nothing
void ITest::init() {}
