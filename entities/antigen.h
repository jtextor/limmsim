#ifndef __ANTIGEN_H

#define __ANTIGEN_H

#include "molecule.h"

class Antigen : public Molecule 
{
  public:
    int * epitopes;
    int * peptides;

    virtual bool operator==(const Entity &other);
    virtual int hashValue(){ return epitopes[0]; }

    // default constructor will create ag with existing, yet undefined epitopes and peptides
    Antigen();
    // this constructor will deep-copy the values of the given arrays
    Antigen(int *, int *);
    ~Antigen();

    virtual Entity * clone();
};

#endif
