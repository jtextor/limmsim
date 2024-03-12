#include "antibody.h"

Antibody::Antibody(){
  paratope = 0;
}

Antibody::Antibody( int _paratope ){
  paratope = _paratope;
}

bool Antibody::operator==(const Entity &other) {
  return paratope == ((Antibody *) &other) -> paratope;
}

inline Entity * Antibody::clone(){
  return (Entity *) new Antibody( paratope );
}
