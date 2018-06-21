#ifndef OBJECT_POOL_HPP
#define OBJECT_POOL_HPP

#include <stack>
#include <vector>

template< typename T > class ObjectPool;

/*!
  \class bulk allocates objects and hands out individual objects as shared pointers
  \param T type of object to manage. Requires T to be default constructible
 */
template< typename T >
class ObjectPoolRaw
{
  public:

    //!\param chunkSize number of objects to allocate at once
    ObjectPoolRaw( size_t chunkSize )
	: mSize( chunkSize )
	, mNewSize( chunkSize )
	, mPosUnused( 0 )
	, mChunks()
	, mUsedObjects()
    {
	allocateNewChunk();
    }

    // \note pool requires have longer lifetime than pointers
    ~ObjectPoolRaw()
    {
	while( !mUsedObjects.empty() )
	    mUsedObjects.pop();
	
	while( !mChunks.empty() )
	{
	    delete [] mChunks.top();
	    mChunks.pop();
	}

	// only deallocate unmanaged, new memory manually
	// for( ; mPosUnused < mSize; ++mPosUnused )
	    // delete (mNewObjects + mPosUnused);
    }

    //! \return pool pointer owning object of type T initialized to undef value
    T* handOut()
    {
	if( mUsedObjects.empty() )
	{
	    if( mPosUnused >= mSize )
		allocateNewChunk();

	    ++mPosUnused;
	    return (mChunks.top() + mPosUnused - 1 );
	}
	
	T* ret = mUsedObjects.top();
	mUsedObjects.pop();
	return ret;
    }

    //! \brief hands pptr back to pool allowing it to be handed out again
    template< typename TT >
    void handBack( TT* pptr )
    {
	mUsedObjects.push( static_cast< T* >( pptr ) );
    }
    
  private:
    void allocateNewChunk()
    {
	mChunks.push( new T[ mSize ]() );
	mSize = mNewSize;
	mPosUnused = 0;
    }

    uint mSize, mNewSize, mPosUnused;
    std::stack< T* > mChunks, mUsedObjects;
};

/*!
  \brief smart pointer with shared semantics akin to std::shared_ptr< T >.  Instead of deallocating memory when the last reference goes out of scope, the memory is returned to the home pool.  The use count is continuously adapted on construction, copy construction, move construction, copy assignment, move assignment and deletion.
 */
template< typename T >
class PoolPointer
{
  public:
    /*!
      \param resource pointer to object to take ownership of
      \param pRefCount owns
      \param home pool to return to if reference count reaches zero
    */
    PoolPointer( T* resource, uint* pRefCount, ObjectPool< T >& home )
	: mRes( resource )
	, mUseCount( pRefCount )
	, mHome( home )
    {}

    PoolPointer( const PoolPointer< T >& orig )
	: mRes( orig.mRes )
	, mUseCount( orig.mUseCount )
	, mHome( orig.mHome )
    {
	incrementCount();
    }

    PoolPointer( const PoolPointer< T >&& orig )
	: mRes( orig.mRes )
	, mUseCount( orig.mUseCount )
	, mHome( orig.mHome )
    {
	orig.mRes = nullptr;
	mUseCount = nullptr;
    }

    ~PoolPointer()
    {
	decrementCount();
    }

    PoolPointer< T >& operator =( const PoolPointer< T >& orig )
    {
	decrementCount();
	
	this->mRes = orig.mRes;
	mUseCount = orig.mUseCount;
	mHome( orig.mHome );

	incrementCount();
	return *this;
    }

    PoolPointer< T >& operator =( const PoolPointer< T >&& orig )
    {
	decrementCount();
	
	this->mRes = orig.mRes;
	mUseCount = orig.mUseCount;
	mHome( orig.mHome );

	orig.mRes = nullptr;
	orig.mUseCount = nullptr;
	
	return *this;
    }

    T& operator *() const
    {
	return *mRes;
    }

    T* operator ->() const
    {
	return mRes;
    }

    //! \return pointer stored
    T* get() const
    {
	return mRes;
    }

    void swap( PoolPointer& other )
    {
	if( this->mRes != other.mRes ) // otherwise: they're the same and swapping makes no difference
	{
	    T* tmpRes = this->mRes;
	    uint* tmpUseCount = this->mUseCount;
	    std::reference_wrapper< ObjectPool< T > > tmpHome( this->mHome );

	    this->mRes = other.mRes;
	    this->mUseCount = other.mUseCount;
	    this->mHome = other.mHome; //ref rebinding possible with reference_wrapper

	    other.mRes = tmpRes;
	    other.mUseCount = tmpUseCount;
	    other.mHome = tmpHome;
	}
    }

    //! \brief releases resources
    void free()
    {
	delete mUseCount;
	delete mRes;
	mUseCount = nullptr;
	mRes = nullptr;
    }
    
  private:

    void decrementCount()
    {
	if( mUseCount != nullptr )
	{
	    --mUseCount;
	    if( mUseCount == 0 )
		mHome.hand_back( mRes );
	}
    }

    void incrementCount()
    {
	if( mUseCount != nullptr )
	    ++(*mUseCount);
    }
    
    T* mRes;
    uint* mUseCount;
    std::reference_wrapper< ObjectPool< T > > mHome;
};

template< typename T >
void swap( PoolPointer< T >& ptr1, PoolPointer< T >& ptr2 )
{
    ptr1.swap( ptr2 );
}

/*!
  \class bulk allocates objects and hands out individual objects as shared pointers
  \param T type of object to manage. Requires T to be default constructible
 */
template< typename T >
class ObjectPool
{
  public:

    //!\param chunkSize number of objects to allocate at once
    ObjectPool( size_t chunkSize )
	: mSize( chunkSize )
	, mNewSize( chunkSize )
	, mPosUnused( 0 )
	, mNewObjects()
	, mUsedObjects()
    {
	allocateNewChunk();
    }

    // \note pool requires have longer lifetime than pointers
    ~ObjectPool()
    {
	// deallocate recycled pointers
	while( !mUsedObjects.empty() )
	{
	    mUsedObjects.top().free();
	    mUsedObjects.pop();
	}
	    
	// only deallocate unmanaged, new memory manually
	for( ; mPosUnused < mSize; ++mPosUnused )
	{
	    delete (mNewObjects + mPosUnused);
	    delete (mNewRefCounts + mPosUnused);
	}
    }

    template< typename U >
    ObjectPool& cast()
    {
	return *this;
    }

    //! \return pool pointer owning object of type T initialized to undef value
    PoolPointer< T > handOut()
    {
	if( mUsedObjects.empty() )
	{
	    if( mPosUnused >= mSize )
		allocateNewChunk();

	    ++mPosUnused;
	    return PoolPointer< T >( mNewObjects + mPosUnused - 1, mNewRefCounts + mPosUnused - 1, std::ref( *this ) );
	}
	
	PoolPointer< T > ret = mUsedObjects.top();
	mUsedObjects.pop();
	return ret;
    }

    //! \brief hands pptr back to pool allowing it to be handed out again
    void handBack( const PoolPointer< T >& pptr )
    {
	mUsedObjects.push( pptr );
    }
    
  private:
    void allocateNewChunk()
    {
	mNewObjects = new T[ mSize ]();
	mNewRefCounts = new uint[ mSize ]();
	mSize = mNewSize;
	mPosUnused = 0;
    }

    uint mSize, mNewSize, mPosUnused;
    T* mNewObjects;
    uint* mNewRefCounts;
    std::stack< PoolPointer< T > > mUsedObjects;
};

#endif

    
