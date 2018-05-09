#ifndef TEST_RUNNER_INTERFACE_HPP
#define TEST_RUNNER_INTERFACE_HPP

#include "testInterface.hpp"

#include <memory>
#include <random>

/*!
   \interface for policy of running a test i.e. specifying the order and frequency in which the methods of a test will be called
*/
struct ITestRunner
{
    /**
       \brief runs the test without taking ownership
    */
    virtual bool run( ITest* pTest ) const = 0;
};

//! \class initializes and runs the test only once, useful for e.g. test groups
struct OnlyOnceRunner : public ITestRunner
{
    // static const std::shared_ptr< const OnlyOnceRunner > mpSelf = make_shared( new OnlyOnceRunner() );
    
    virtual bool run( ITest* pTest ) const;
};

//! \class carries out required number of repetitions; sequence per repetition: initialize then iterate according to test size and finally check once
struct ContinuousRunner : public ITestRunner
{
    // static const std::shared_ptr< const ContinuousRunner > mpSelf = make_shared( new ContinuousRunner() );

    //! \return return result of final check
    virtual bool run( ITest* pTest ) const;
};

//! \class carries out required number of repetitions; sequence per repetition: initialize then iterate any number of times between zero or test size and finally check once
struct ContinuousRandomRunner : public ITestRunner
{
    // static const std::shared_ptr< const ContinuousRandomRunner > mpSelf = make_shared( new ContinuousRandomRunner() );

    static std::default_random_engine mRandom;

    //! \return return result of final check
    virtual bool run( ITest* pTest ) const;
};

//! \class carries out required number of repetitions; sequence per repetition: initialize once, then alternatingly call iterate and check as often as test size
struct InterleaveRunner : public ITestRunner
{
    // static const std::shared_ptr< const InterleaveRunner > mpSelf = make_shared( new InterleaveRunner );

    //! \return true if check always held, false otherwise
    virtual bool run( ITest* pTest ) const;
};

//! \class carries out required number of repetitions; sequence per repetition: initialize once, then alternatingly call iterate and check any number of times between zero and test size
struct InterleaveRandomRunner : public ITestRunner
{
    // static const std::shared_ptr< const InterleaveRandomRunner > mpSelf = make_shared( new InterleaveRandomRunner );

    static std::default_random_engine mRandom;

    //! \return true if check always held, false otherwise
    virtual bool run( ITest* pTest ) const;
};

#endif
