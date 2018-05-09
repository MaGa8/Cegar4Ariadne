#ifndef ALL_TESTS_HPP
#define ALL_TESTS_HPP

#include "linkedFixedBranchTreeTest.hpp"
#include "adjacencyDiGraphTest.hpp"

struct AllTests : public ITestGroup
{
    GROUP_CTOR_DECL( AllTests );
    
    void init();
};






#endif







