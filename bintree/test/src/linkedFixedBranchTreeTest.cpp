#include "linkedFixedBranchTreeTest.hpp"
#include "testRunnerInterface.hpp"
#include "testMacros.hpp"

#ifndef DEBUG
#define DEBUG false
#endif

using namespace tree;

std::default_random_engine LinkedFixedBranchTreeTest::mRandom = std::default_random_engine( std::random_device()() );
std::uniform_int_distribution<> LinkedFixedBranchTreeTest::mIntDist = std::uniform_int_distribution<>();
std::uniform_int_distribution<> LinkedFixedBranchTreeTest::mChildDist = std::uniform_int_distribution<>( 0, LinkedFixedBranchTreeTest::mNoChildren - 1 );
std::uniform_real_distribution<> LinkedFixedBranchTreeTest::mProbDist = std::uniform_real_distribution<>( 0, 1 );

LinkedFixedBranchTreeTest::TEST_CTOR( LeafTest, "test isLeaf" )

void LinkedFixedBranchTreeTest::LeafTest::init()
{
    mpTree.reset( new TestTreeT( mIntDist( mRandom ) ) );
}

void LinkedFixedBranchTreeTest::LeafTest::iterate()
{
    randomExpandTree( *mpTree, mIntDist );
    D( std::cout << "expanded: size = " << mpTree->size() << ", height = " << mpTree->height() << std::endl; );
}

// look for some leaf and test that it does not have any children
bool LinkedFixedBranchTreeTest::LeafTest::check() const
{
    TestTreeT::NodeT n = root( *mpTree );
    while( !isLeaf( *mpTree, n ) )
    {
    	auto crange = children( *mpTree, n );
    	std::advance( crange.first, mChildDist( mRandom ) );
    	n = *crange.first;
    }
    typename FixedBranchTreeTraits< TestTreeT >::CRangeT leafRange = children( *mpTree, n );
    if( leafRange.first != leafRange.second )
    {
    	D( std::cout << "erroneous children: " << std::distance( leafRange.first, leafRange.second ) << std::endl; );
    	return false;
    }
    return true;
}

LinkedFixedBranchTreeTest::TEST_CTOR( ExpandSizeTest, "size of tree" )

void LinkedFixedBranchTreeTest::ExpandSizeTest::init()
{
    mNoExpansions = 0;
    mpTree.reset( new TestTreeT( mIntDist( mRandom ) ) );
}

void LinkedFixedBranchTreeTest::ExpandSizeTest::iterate()
{
    randomExpandTree( *mpTree, mIntDist );
    ++mNoExpansions;
}
    
bool LinkedFixedBranchTreeTest::ExpandSizeTest::check() const
{
    if( mpTree->size() != mNoChildren * mNoExpansions + 1 )
	return false;
    return true;
}

LinkedFixedBranchTreeTest::TEST_CTOR( ExpandHeightTest, "height of tree" )

void LinkedFixedBranchTreeTest::ExpandHeightTest::init()
{
    mNoExpansions = 0;
    mpTree.reset( new TestTreeT( mIntDist( mRandom ) ) );
}

void LinkedFixedBranchTreeTest::ExpandHeightTest::iterate()
{
    randomExpandTree( *mpTree, mIntDist );
    ++mNoExpansions;
}

bool LinkedFixedBranchTreeTest::ExpandHeightTest::check() const
{
    /* min: perfect n-ary tree => size( t ) = b^height( t ) + 1
       max: number of expansions */
    if( mpTree->height() < std::log2( mpTree->size() - 1 ) / std::log2( mNoChildren )
	&& mpTree->height() > mNoExpansions )
	return false;

    return true;
}

LinkedFixedBranchTreeTest::TEST_CTOR( RootTest,  "test isRoot on deep tree" )

void LinkedFixedBranchTreeTest::RootTest::init()
{
    mpTree.reset( new TestTreeT( mIntDist( mRandom ) ) );
    mRoot = root( *mpTree );
}

void LinkedFixedBranchTreeTest::RootTest::iterate()
{
    randomExpandTree( *mpTree, mIntDist );
}

bool LinkedFixedBranchTreeTest::RootTest::check() const
{
    if( !isRoot( *mpTree, mRoot ) )
    {
	return false;
    }
    
    return true;
}

LinkedFixedBranchTreeTest::TEST_CTOR( NotRootTest, "test complement of isRoot on deep tree" )

void LinkedFixedBranchTreeTest::NotRootTest::init()
{
    mpTree.reset( new TestTreeT( mIntDist( mRandom ) ) );
}

void LinkedFixedBranchTreeTest::NotRootTest::iterate()
{
    randomExpandTree( *mpTree, mIntDist );
}

bool LinkedFixedBranchTreeTest::NotRootTest::check() const
{
    typename TestTreeT::NodeT n = root( *mpTree );
    while( !isLeaf( *mpTree, n ) )
    {
	auto crange = children( *mpTree, n );
	std::advance( crange.first, mChildDist( mRandom ) );
	n = *crange.first;
	if( isRoot( *mpTree, n ) )
	{
	    D( std::cout << "non-root node " << value( *mpTree, n ) << " identified as root " << std::endl; );
	    return false;
	}
    }
    return true;
}

LinkedFixedBranchTreeTest::TEST_CTOR( DeleteTest, "deletion of children" )

void LinkedFixedBranchTreeTest::DeleteTest::init()
{
    mpTree.reset( new TestTreeT( mIntDist( mRandom ) ) );
    mDeletedAt.reset();
    for( uint cExpand = 0; cExpand < mTestSize; ++cExpand )
	randomExpandTree( *mpTree, mIntDist );
}

void LinkedFixedBranchTreeTest::DeleteTest::iterate()
{
    randomExpandTree( *mpTree, mIntDist );
    TestTreeT::NodeT n = root( *mpTree );
    while( !isLeaf( *mpTree, n ) && mProbDist( mRandom ) > mTraverseThresh )
    {
	auto crange = children( *mpTree, n );
	std::advance( crange.first, mChildDist( mRandom ) );
	n = *crange.first;
    }
    mDeletedAt = std::make_optional( n );
    delChildren( *mpTree, n );
}

bool LinkedFixedBranchTreeTest::DeleteTest::check() const
{
    if( !isLeaf( *mpTree, mDeletedAt.value_or( root( *mpTree ) ) ) )
    {
	std::cout << value( *mpTree, mDeletedAt.value() ) << " is not a child even though deletion took place here" << std::endl;
	return false;
    }
    return true;
}

LinkedFixedBranchTreeTest::TEST_CTOR( DeleteSizeTest, "effect of deletion on size" )

void LinkedFixedBranchTreeTest::DeleteSizeTest::init()
{
    mpTree.reset( new TestTreeT( mIntDist( mRandom ) ) );
    for( uint cExpand = 0; cExpand < mTestSize; ++cExpand )
	randomExpandTree( *mpTree, mIntDist );
}

void LinkedFixedBranchTreeTest::DeleteSizeTest::iterate()
{
    randomExpandTree( *mpTree, mIntDist );
    TestTreeT::NodeT n = root( *mpTree );
    while( !isLeaf( *mpTree, n ) && mProbDist( mRandom ) > mTraverseThresh )
    {
	auto crange = children( *mpTree, n );
	std::advance( crange.first, mChildDist( mRandom ) );
	n = *crange.first;
    }
    delChildren( *mpTree, n );
}

bool LinkedFixedBranchTreeTest::DeleteSizeTest::check() const
{
    size_t subSize = subtreeSize( *mpTree, root( *mpTree ) );
    if( subSize != mpTree->size() )
    {
	D( std::cout << "incorrect size after deletion " << mpTree->size() << " conflicts with size of complete tree " << subSize << std::endl; );
	return false;
    }
    return true;
}

LinkedFixedBranchTreeTest::TEST_CTOR( DeleteHeightTest, "effect of deletion on height" )

void LinkedFixedBranchTreeTest::DeleteHeightTest::init()
{
    mpTree.reset( new TestTreeT( mIntDist( mRandom ) ) );
    for( uint cExpand = 0; cExpand < mTestSize; ++cExpand )
	randomExpandTree( *mpTree, mIntDist );
    mHeightBefore = mpTree->height();
    D( std::cout << "init delete height test: tree height " << mHeightBefore << std::endl; );
}

void LinkedFixedBranchTreeTest::DeleteHeightTest::iterate()
{
    mHeightBefore = mpTree->height();
    TestTreeT::NodeT n = root( *mpTree );
    while( !isLeaf( *mpTree, n ) && mProbDist( mRandom ) > mTraverseThresh )
    {
	auto crange = children( *mpTree, n );
	std::advance( crange.first, mChildDist( mRandom ) );
	n = *crange.first;
    }
    delChildren( *mpTree, n );
    D( std::cout << "iterate delete height test: tree height was " << mHeightBefore << std::endl; );
}

bool LinkedFixedBranchTreeTest::DeleteHeightTest::check() const
{
    size_t heightNow = mpTree->height();
    if( heightNow > mHeightBefore )
    {
	D( std::cout << "incorrect height after deletion " << mpTree->height() << " whereas before this was " << mHeightBefore << std::endl; );
	return false;
    }
    return true;
}

LinkedFixedBranchTreeTest::MemoryFreed::Dummy::Dummy()
    : mDefaultConstructed( true ) {}

LinkedFixedBranchTreeTest::MemoryFreed::Dummy::Dummy( uint *counterRef, const uint& id )
    : mCounterRef( counterRef ), mId( id ), mDefaultConstructed( false ) { ++(*mCounterRef); }

LinkedFixedBranchTreeTest::MemoryFreed::Dummy::Dummy( const Dummy& orig )
    : mCounterRef( orig.mCounterRef ), mId( orig.mId ), mDefaultConstructed( orig.mDefaultConstructed )
{
    if( !mDefaultConstructed )
	++(*mCounterRef);
}

LinkedFixedBranchTreeTest::MemoryFreed::Dummy& LinkedFixedBranchTreeTest::MemoryFreed::Dummy::operator =( const Dummy& orig )
{
    mId = orig.mId; // no reference rebinding, so leave counter as is
    return *this;
}

LinkedFixedBranchTreeTest::MemoryFreed::Dummy::~Dummy()
{
    if( !mDefaultConstructed )
	--(*mCounterRef);
}

bool LinkedFixedBranchTreeTest::MemoryFreed::Dummy::operator ==( const Dummy& other ) { return this->mId == other.mId; }

LinkedFixedBranchTreeTest::MemoryFreed::MemoryFreed( const uint& testSize, const uint& repetitions )
    : ITest( "test tree nodes are freed on destruction of tree", testSize, repetitions )
    , mSizeDist( 0, testSize )
{}

void LinkedFixedBranchTreeTest::MemoryFreed::MemoryFreed::iterate()
{
    mObjectCounter = 0;
    mCreationCounter = 0;
    mpTree.reset( new LinkedFixedBranchTree< Dummy, mNoChildren >( Dummy( &mObjectCounter, ++mCreationCounter ) ) );

    for( uint i = 0; i < mSizeDist( mRandom ); ++i )
    {
	auto node = tree::root( *mpTree );
	while( !tree::isLeaf( *mpTree, node ) )
	{
	    auto crange = tree::children( *mpTree, node );
	    std::advance( crange.first, mIntDist( mRandom ) % mNoChildren );
	    node = *crange.first;
	}
	std::array< Dummy, mNoChildren > exvals = {Dummy( &mObjectCounter, ++mCreationCounter )
						   , Dummy( &mObjectCounter, ++mCreationCounter )
						   , Dummy( &mObjectCounter, ++mCreationCounter )
						   , Dummy( &mObjectCounter, ++mCreationCounter ) };
	tree::expand( *mpTree, node, exvals );
    }
    mpTree.reset();
}
    
bool LinkedFixedBranchTreeTest::MemoryFreed::MemoryFreed::check() const
{
    if( mObjectCounter != 0 )
    {
	std::cout << "out of " << mCreationCounter << " dummy values created "
		  << mObjectCounter << " are still not deallocated after deletion of tree" << std::endl;
	return false;
    }
    return true;
}

GROUP_CTOR( LinkedFixedBranchTreeTest, "linked fixed branch tree" )

void LinkedFixedBranchTreeTest::init()
{
    std::shared_ptr< InterleaveRandomRunner > printerleave( new InterleaveRandomRunner() );
    std::shared_ptr< StatelessRunner > pStateless( new StatelessRunner() );

    addTest( new LeafTest( mTestSize, mRepetitions ), printerleave );
    addTest( new ExpandSizeTest( mTestSize, mRepetitions ), printerleave );
    addTest( new ExpandHeightTest( mTestSize, mRepetitions ), printerleave );
    addTest( new RootTest( mTestSize, mRepetitions ), printerleave );
    addTest( new NotRootTest( mTestSize, mRepetitions ), printerleave );
    addTest( new DeleteTest( 0.25 * mTestSize, mRepetitions ), printerleave );
    addTest( new DeleteSizeTest( 0.5 * mTestSize, mRepetitions ), printerleave );
    addTest( new DeleteHeightTest( 0.25 * mTestSize, mRepetitions ), printerleave );
    addTest( new MemoryFreed( 8 * mTestSize, mRepetitions ), pStateless );
}
