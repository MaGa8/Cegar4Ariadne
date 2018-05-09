#include <sys/ioctl.h> // linux only!
#include <unistd.h>
#include <iomanip>
#include <algorithm>
#include "testRunner.hpp"

ITestRunner::ITestRunner( std::string name ) : mName( name ) {}

ITestGroup::ITestGroup( std::string name, uint level, std::ostream& out) : ITestRunner( name ), mOut( out ), mLevel( level ) {}

uint ITestGroup::getLevel() const
{
  return mLevel;
}

bool ITestGroup::run()
{
  struct winsize termSize;
  ioctl( STDOUT_FILENO, TIOCGWINSZ, &termSize );
  uint termcols = termSize.ws_col, shift = -40;
  bool success = true;

  for(std:: unique_ptr< ITestRunner >& pTest : mQueue )
    {
      bool successTest = pTest->run();
      success = success && successTest;

      uint midFill =  termcols + shift - pTest->mName.size();
      mOut << std::setw( getLevel() * LEVEL_INDENT ) << " "
	   << pTest->mName << " : ";
      mOut << std::setw( midFill ) << " "
	   << ( successTest ? " success" : " failure")
	   << std::endl;
    }
  return success;
}

void ITestGroup::addTest( ITestRunner* pTest )
{
  mQueue.push_back( std::unique_ptr< ITestRunner >( pTest ) );
}
