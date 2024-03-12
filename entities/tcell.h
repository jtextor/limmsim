#ifndef __TCELL_H

#define __TCELL_H

#include "cell.h"

class TCell : public Cell {
  public:
    // bitstring of t-cell mhcII-ag-binding (cd4) receptor
    int cd4;
    // number of times this b-cell has divided
    // relevant iff state == DUPLICATING
    int dupcycles; 
    // age of this cell
    // relevant for memory cells only, non-memory cells die uniformly
    int age;
    // thus, a cell is a memory cell iff its age is > 0
    inline bool isMemory(){ return age > 0; }
    
    virtual bool operator==(const Entity &);
    virtual int hashValue(){ return cd4; }

    TCell();
    TCell( int );

    virtual Entity * clone();

    // STATES 
    const static short NAIVE = 0;
    const static short DUPLICATING = 1;
    const static short EFFECTOR = 2;
    const static short MEMORY = 3;
};

#endif
