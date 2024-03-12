#include "immunecomplex.h"

inline bool ImmuneComplex::operator==(const Entity &other) {
  return true;
}

Entity * ImmuneComplex::clone(){
  return (Entity * ) new ImmuneComplex();
}
