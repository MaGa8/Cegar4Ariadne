#ifndef TEST_INTERFACE_HPP
#define TEST_INTERFACE_HPP

#include <string>

#define STATELESS_METHODS 

#define STATELESS_TEST(name) public: name( uint size, uint repetitions ); void iterate(); bool check() const

#define STATEFUL_TEST(name) public: name( uint size, uint repetitions ); void init(); void iterate(); bool check() const

#define TEST_CTOR(name,description) name::name( uint size, uint repetitions ) : ITest(description, size, repetitions) {}

/**
   \interface encapsulating a test that can be run
*/
class ITest
{
  public:
    /*!
       \param name string to print to identify the test
       \param reps number of times to repeat the test
       \param testSize complexity of the test to be carried out
    */
    ITest( std::string name, uint size, uint repetitions );

    //! \brief restores initial state which should be such that check returns true
    virtual void init();

    //! \brief given an arbitrary present state, execute an iteration of the test
    virtual void iterate() = 0;

    //! \brief checks the state for validity
    virtual bool check() const = 0;

    /** name of the test */
    const std::string mName;
    const uint mTestSize, mRepetitions;

  private:
    ITest( const ITest& copy ) = delete;
};

#endif
