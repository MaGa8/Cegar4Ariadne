
#include "linkedFixedBranchTree.hpp"

int main()
{
    typedef LinkedFixedBranchTree< int, 2 > Bintree;
    Bintree binTree( 0 );

    typename FixedBranchTreeTraits< Bintree >::NodeT rt = root( binTree );
    std::cout << "value " << value( binTree, rt ) << std::endl;
    expand( binTree, rt );

    // auto crange = children( binTree, rt );

    std::pair< FixedBranchTreeTraits< Bintree >::CIterT, FixedBranchTreeTraits< Bintree >::CIterT > crange = children( binTree, rt );

    for( auto ci = crange.first; ci < crange.second; ++ci )
    {
	FixedBranchTreeTraits< Bintree >::NodeT n = *ci;
    	std::cout << value( binTree, n ) << std::endl;
    }

    uint cnt = 0;
    for( FixedBranchTreeTraits< Bintree >::CIterT& ci = crange.first; ci < crange.second; ++ci )
    {
	FixedBranchTreeTraits< Bintree >::NodeT n = *ci;
	int& v = value( binTree, n );
	v = ++cnt;
    }

    // make iterators copyable!
    crange = children( binTree, rt );
    
    for( FixedBranchTreeTraits< Bintree >::CIterT& ci = crange.first; ci < crange.second; ++ci )
    	std::cout << value( binTree, *ci ) << std::endl;

    crange = children( binTree, rt );
    
    FixedBranchTreeTraits< Bintree >::NodeT n = *(children( binTree, rt ).first );
    FixedBranchTreeTraits< Bintree >::NodeT npt = parent( binTree, n );
    std::cout << "parent " << value( binTree, npt ) << std::endl;

    return 0;
}
