#include "bcell.h"

#include <iostream>

BCell::BCell(){
  state = NAIVE;
  receptor = 0;
  isMem = false;
  dupcycles = 0;
  MHCIIpep = 0;
}

BCell::BCell(int _receptor){
  state = NAIVE;
  receptor = _receptor;
  isMem = false;
  dupcycles = 0;
  MHCIIpep = 0;
}

bool BCell::operator==(const Entity &other) {

  if( state != ((BCell *) &other) -> state ) return false;
  if( receptor != ((BCell *) &other) -> receptor ) return false;
  if( isMem != ((BCell *) &other) -> isMem ) return false;
  if( dupcycles != ((BCell *) &other) -> dupcycles ) return false;
  if( MHCIIpep != ((BCell *) &other) -> MHCIIpep ) return false;

  return true;
}

Entity * BCell::clone(){
  BCell * ret = new BCell();
  ret -> receptor = receptor;
  ret -> state = state;
  ret -> MHCIIpep = MHCIIpep;
  ret -> isMem = isMem;
  ret -> ag = ag;
  ret -> dupcycles = dupcycles;

  return (Entity *) ret;
}

