#ifndef __BCELL_H

#define __BCELL_H

#include "cell.h"
#include "antigen.h"

class BCell : public Cell {
  public:
    // bitstring of ag-binding receptor
    int receptor;
    // bitstring of peptide bound to mhc II - molecule 
    // relevant iff state == PRES_II
    int MHCIIpep;
    // number of times this b-cell has divided
    // relevant iff state == DUPLICATING
    int dupcycles; 
    // age of this cell
    // relevant for memory cells only, non-memory cells die uniformly
    // *** change: don't keep track of age. increasing half life 
    // already has the desired effect. *** 
    //int age;
    bool isMem;
    // thus, a cell is a memory cell iff its age is > 0
    inline bool isMemory(){ return isMem; }
    
    // Antigen bound to this cell
    // relavant iff state == INTERNALIZED
    Antigen * ag;

    virtual bool operator==(const Entity &);
    virtual int hashValue(){ return receptor; }

    BCell();
    BCell(int);

    virtual Entity * clone();

    // STATES 
    const static short NAIVE = 0;
    const static short INTERNALIZED = 1;
    const static short PRES_II = 2;
    const static short DUPLICATING = 3;
    const static short PLASMA = 4;
};

#endif 
