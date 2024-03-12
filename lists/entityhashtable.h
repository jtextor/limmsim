#ifndef __ENTITY_HASHTABLE

#define __ENTITY_HASHTABLE

#include "entitylist.h"
#include "../entities/entity.h"
#include <pthread.h>

template<class T> class EntityHashtable
{
  public: 
    EntityHashtable(){
		for( int i = 0 ; i < size ; i ++ ){
			pthread_mutex_init (&list_mutex[i], 0 );
		}
	};

    T * find( T * e ){
		return lists[e->hashValue()%size].find( e );
	}

    // Methods below are not threadsafe. 
	void add( T * e, int _count = 1 ){
		lists[e->hashValue()%size].add( e, _count );
	};
    void remove( T * e , int _count = 1 ){
		int h = e->hashValue()%size;
		
		// ** begin critical region
		pthread_mutex_lock(&list_mutex[h]);
		
		T * e2 = lists[e->hashValue()%size].remove( e, _count );
		if( e2 != 0 ) delete e2;
		
		pthread_mutex_unlock(&list_mutex[h]);
		// ** end critical region
	};

    // store is a threadsafe add. If the entity is not yet in the hastable,
    // a new entity with equal attributes is created and stored. A pointer to the new
    // entity is then returned.
    // store uses a mutex on the internal list corresponding to the Entity's hashValue.
    T * store( T * e, int _count = 1 ){
		int h = e->hashValue()%size;
		
		/** begin critical region */
		pthread_mutex_lock(&list_mutex[h]);
		
		T * e2 = lists[h].find( e );
		if( e2 == 0 )
			e2 = (T * ) e -> clone();
		lists[h].add( e2, _count );
		
		pthread_mutex_unlock(&list_mutex[h]);
		/** end critical region */
		
		return e2;
	};

    // count ALL elements in this hashtable.
	int count_elements(){
		int c = 0;
		for( int i = 0; i < size ; i ++ ) 
			c += lists[i].count_elements();
		return c;
	};

    // iterating through the hashtable is not threadsafe.
    void init_iterator(){
		_current_hashvalue = 0;
		while( _current_hashvalue < size && lists[_current_hashvalue].length() == 0 ) _current_hashvalue++;
		if( _current_hashvalue < size ){
			_current_entity = lists[_current_hashvalue].item(0)->item;
			_current_count = lists[_current_hashvalue].item(_current_position)->count;
		}
		else _current_entity = 0;
		_current_position = 0;
	};

    T * next(){
		T * e = _current_entity;
		if( e != 0 ){
			if( _current_position < lists[_current_hashvalue].length() - 1 ){
				_current_position ++;
				_current_entity = lists[_current_hashvalue].item(_current_position)->item;
				_current_count = lists[_current_hashvalue].item(_current_position)->count;
			}
			else{
				_current_hashvalue++;
				_current_position = 0;
				while( _current_hashvalue < size && lists[_current_hashvalue].length() == 0 ) _current_hashvalue++;	
				if( _current_hashvalue < size ){ 
					_current_entity =  lists[_current_hashvalue].item(0)->item;
					_current_count = lists[_current_hashvalue].item(_current_position)->count;
				}
				else _current_entity = 0;
			}
		}
		return e;
	};

    inline T * current_entity(){ return _current_entity; } 
    inline int current_count(){ return _current_count; }

    // the hash table should not habe too many cells
	// since a pthread_mutex has to be created for each one. 
    // on linux, libpthread crashes for size = 1 << 12 
    const static int size = 1 << 8;
    pthread_mutex_t list_mutex[size];
    EntityList<T> lists[size];

  private: 
    int _current_hashvalue;
    int _current_position;
    T * _current_entity; 
    int _current_count;
};

#endif
