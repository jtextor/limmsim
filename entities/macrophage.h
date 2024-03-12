#ifndef __MACROPHAGE

#define __MACROPHAGE

#include "cell.h"
#include "antigen.h"

class Macrophage : public Cell {
  public:
    // bitstring of peptide bound to mhc II - molecule 
    // relevant iff state == PRES_II
    int MHCIIpep;
    // Antigen bound to this cell
    // relavant iff state == INTERNALIZED
    Antigen * ag;

    virtual bool operator==(const Entity &);
    virtual int hashValue(){ return MHCIIpep; }

    Macrophage();

    virtual Entity * clone();

    // STATES 
    const static short ACTIVE = 0;
    const static short INTERNALIZED = 1;
    const static short PRES_II = 2;
    const static short MEMORY = 3;
};

#endif
