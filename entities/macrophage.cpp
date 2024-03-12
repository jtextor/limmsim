#include "macrophage.h"

#include <iostream>

bool Macrophage::operator==(const Entity &other) {
  if( state != ((Macrophage *) &other) -> state ) return false;
  if( state == PRES_II && (MHCIIpep != ((Macrophage *) &other) -> MHCIIpep) ) return false;
  if( state == INTERNALIZED && !(ag == (((Macrophage *) &other)->ag)) ) return false;

  return true;
}

Macrophage::Macrophage(){
  state = ACTIVE;
  ag = 0;
  MHCIIpep = 0;
}

Entity * Macrophage::clone() {
  Macrophage * ret = new Macrophage();
  ret -> state = state;
  ret -> ag = ag;
  ret -> MHCIIpep = MHCIIpep;

  return (Entity *) ret;
}
