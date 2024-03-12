#include "tcell.h"

TCell::TCell( int _cd4 ){
  state = NAIVE;
  cd4 = _cd4;
  age = 0;
  dupcycles = 0;
}

TCell::TCell(){
  state = NAIVE;
  cd4 = 0;
  age = 0;
  dupcycles = 0;
}

bool TCell::operator==(const Entity &other) {
  if( cd4 != ((TCell *) &other) -> cd4 ) return false;
  if( state != ((TCell *) &other) -> state ) return false;
  if( dupcycles != ((TCell *) &other) -> dupcycles ) return false;
  if( age != ((TCell *) &other) -> age ) return false;
 
  return true;
}

Entity * TCell::clone(){
  TCell * ret = new TCell();
  ret -> state = state;
  ret -> cd4 = cd4;
  ret -> age = age;
  ret -> dupcycles = dupcycles;

  return (Entity *) ret;
}
