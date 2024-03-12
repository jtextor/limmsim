#include "molecule.h"

class Antibody : public Molecule 
{
  public:
    int paratope;

    virtual bool operator==(const Entity &other) ;
    virtual int hashValue(){ return paratope; }

    Antibody();
    Antibody(int);

    virtual Entity * clone();

};

