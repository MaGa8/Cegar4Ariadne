#ifndef TEST_GROUP_INTERFACE
#define TEST_GROUP_INTERFACE

#include "testInterface.hpp"
#include "testRunnerInterface.hpp"

#include <map>
#include <memory>
#include <iostream>

#define GROUP_CTOR_DECL(name) name( uint size, uint repetitions = 1, uint level = 0, std::ostream& out = std::cout )

#define GROUP_CTOR(name,description) name::name( uint size, uint repetitions, uint level, std::ostream& out ) : ITestGroup( description, size, repetitions, level, out ) {}

/**
   \class base class for a group of tests
*/
class ITestGroup : public ITest
{
  public:
    // column at which to print result of tests
    const uint RESULT_PRINT_POS = 100;
    const uint LEVEL_INDENT = 4;
    const uint mLevel;
    
    /*!
      \param name string to print to identify group
      \param size complexity of tests to pass down to children
      \param repetitions number of times to repeat all tests
      \param level level in hierarchy (for indentation)
      \param out stream to write output to
    */
    ITestGroup( std::string name, uint size, uint repetitions = 1, uint level = 0, std::ostream& out = std::cout );

    //! \brief set up tests
    virtual void init() = 0;

    //! \brief execute all tests once according to policy, prints and stores result
    void iterate();

    //! \return outcome of tests
    bool check() const;

    /**
       \brief adds test to the queue
       \param test pointer to test to perform
       \param 
    */
    void addTest( ITest* pTest, std::shared_ptr< ITestRunner > pRunner );
    
  protected:
    std::ostream& mOut;
    
  private:
    typedef std::pair< std::unique_ptr< ITest >
		       , std::shared_ptr< ITestRunner > > TestPolicyPair;
    
    ITestGroup( const ITestGroup& copy ) = delete;

    std::vector< TestPolicyPair > mTests;
    bool mResult;
};

#endif
