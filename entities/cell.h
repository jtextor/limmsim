#ifndef __CELL_H

#define __CELL_H

#include "entity.h"

class Cell : public Entity {
  public:
    // every cell has an internal state.
    int state; 
};

#endif
