#ifndef SIBLING_ITERATOR
#define SIBLING_ITERATOR

template< typename NodeT >
class ChildrenIteratorInterface
{
  public:
    virtual NodeT operator * () = 0;

    // virtual const NodeT operator * () const = 0;

    virtual void operator ++ () = 0;

    virtual void operator -- () = 0;
};

#endif
