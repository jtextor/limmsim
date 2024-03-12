#include "../lists/entitylist.h"

#ifndef __GRIDPOINT_H

#define __GRIDPOINT_H

class GridPoint
{
  public:
    GridPoint ** neighbours;

    // entities at this grid point
    // all pointers in these lists are _shared_ with global lists ! 
    EntityList<Antigen> ag;
    EntityList<Antibody> ab;
    EntityList<ImmuneComplex> ic;
    EntityList<TCell> t;
    EntityList<BCell> b;
    EntityList<Macrophage> ma;

    GridPoint();
    ~GridPoint();
};

typedef EntityList<Entity> GridPoint::* listPtr;


#endif
