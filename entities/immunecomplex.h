#include "molecule.h"

class ImmuneComplex : public Molecule 
{
  public:
    // all immunecomplexes are equal :) their sole destiny is to get eaten by a macrophage 
    virtual bool operator==(const Entity &other) ;
    virtual int hashValue(){ return 0; }

    virtual Entity * clone();
};

