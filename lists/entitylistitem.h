#ifndef __ENTITYLISTITEM
#define __ENTITYLISTITEM

#include "../entities/entity.h" 

template<class T> class EntityListItem{ 
  public:
    T * item;
    int count;
    EntityListItem(){
      count = 0;
    }
};

#endif
