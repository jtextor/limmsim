#include "gridpoint.h"
#include "iostream"

GridPoint::GridPoint(){
  neighbours = new GridPoint * [7];
}

GridPoint::~GridPoint(){
  delete [] neighbours;
}
