
#include <map>


template< typename T >
class MapNTree
{
  public:
    typename long KeyT;





  private:
    std::map< KeyT, T > mNodeMap;
};
