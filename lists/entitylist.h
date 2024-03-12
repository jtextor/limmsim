#ifndef __ENTITYLIST

#define __ENTITYLIST

#include "entitylistitem.h"
#include "../settings/settings.h"
#include "../threadlocals/threadlocals.h"

template<class T> class EntityList
{
  public: 
    EntityList(){
		items = new EntityListItem<T> * [Settings::initial_listsize];
		for( int i = 0 ; i < Settings::initial_listsize ; i ++ ) items[i]=0;
		max_size = Settings::initial_listsize;
		current_size = 0;
	};
    ~EntityList(){
		delete[] items;
	};

	// find uses biological identity
    T * find( T * e ){
		for( int i = 0 ; i < current_size ; i ++ ){
			if( *(items[i]->item) == *e )
			return items[i]->item;
		}
		return 0;
	};

	// all methods below use pointer identity. 
    void add( T * e , int _count = 1 ){
		for( int i = 0 ; i < current_size ; i++ )
			if( items[i]->item == e )
			{ items[i]->count += _count; return; }

		if( current_size == max_size ) grow();
		if( items[current_size] == 0 ) items[current_size] = new EntityListItem<T>();

		items[current_size]->item = e;
		items[current_size]->count = _count;

		current_size ++;
		return;
	};

    int count_elements(){
		int ret = 0;
		for( int i = 0 ; i < current_size ; i ++ )
			ret += items[i]->count;
		return ret;
	};

    inline int length(){ return current_size; }

    inline EntityListItem<T> * item( int i ){ return items[i]; }

    // shuffle list
    void shuffle( int thread_nr ){
		EntityListItem<T> * tmp; int dest;
		for( int i = 0 ; i < current_size - 1 ; i ++ ){
			dest = ThreadLocals::rng[thread_nr].randInt( current_size - 1 );
			tmp = items[dest];
			items[dest] = items[i];
			items[i] = tmp;
		}
	};

    // if the removed entity was the last one of its kind in the list,
    // a pointer to this entity is returned. the caller can then decide if 
    // the entity should be deleted (normally only when called from the global
    // hashtables) 
    T * remove( T * e , int _count = 1 ){
		for( int i = 0 ; i < current_size ; i++ )
			if( items[i]->item == e )
			{ 
				items[i]->count -=_count;
				if( items[i]->count > 0 )
				return 0;
				items[i]->item = items[--current_size]->item;
				items[i]->count = items[current_size]->count;
				return e;
			}
		return 0;
	};

  private:
    EntityListItem<T> **items;
    int max_size;
    int current_size;

    // reallocate array if it is full
    void grow(){
		int i;
		max_size += Settings::initial_listsize;
		EntityListItem<T> ** _newitems = new EntityListItem<T> * [max_size];
		for( i = 0 ; i < current_size ; i ++ )
			_newitems[i] = items[i];
		for( ; i < max_size ; i ++ )
			_newitems[i] = 0;
		delete [] items;
		items = _newitems;
	};
};

#endif
