#ifndef TEST_RUNNER_H
#define TEST_RUNNER_H

#include <iostream>
#include <iomanip>
#include <queue>
#include <memory>

#define TEST_STRUCT(name) struct name : public ITestRunner { name(); bool run(); };

#define D(x) do{ if( DEBUG ) x } while( false )

/**
   \interface encapsulating a test that can be run
*/
class ITestRunner
{
  public:
    ITestRunner( std::string name );

    /**
       \brief runs the test
    */
    virtual bool run() = 0;

    /** name of the test */
    const std::string mName;

  private:
    ITestRunner( const ITestRunner& copy ) = delete;
};

/**
   \class base class for a group of tests
*/
class ITestGroup : public ITestRunner
{
  public:
    const uint LEVEL_INDENT = 4;
    
    ITestGroup( std::string name, uint level = 0, std::ostream& out = std::cout );

    uint getLevel() const;

    /**
       \brief runs the tests stored in order, printing the status, and clears the test queue afterwards
       \return true if all tests succeeded, false if at least one failed
       \note to be invoked by derived classes after they populated the queue. Populating should happen anew for each run because the queue is cleared
    */
    virtual bool run();

    /**
       \brief adds test to the queue
       \param test pointer to test to perform. Will take ownership of test
    */
    void addTest( ITestRunner* pTest );

  protected:
    std::ostream& mOut;
    
  private:
    ITestGroup( const ITestGroup& copy ) = delete;

    ITestGroup operator =( const ITestGroup& copy ) = delete;

    std::vector< std::unique_ptr< ITestRunner > > mQueue;
    uint mLevel;
};

  
#endif
