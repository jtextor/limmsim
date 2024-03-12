#ifndef __THREADLOCALS_H

#define __THREADLOCALS_H

#include "../mersennetwister/mersennetwister.h" 
#include "../settings/settings.h" 
#include "../entities/tcell.h" 
#include "../entities/bcell.h" 
#include "../entities/antibody.h" 
#include "../entities/macrophage.h" 
#include "../entities/immunecomplex.h" 

class ThreadLocals{
  public: 
    static pthread_t * threads;
    static int * message;

    static MTRand * rng;
    static TCell * tcell;
    static Macrophage * ma;
    static BCell * bcell;
    static Antibody * ab;
    static ImmuneComplex * ic;

    static int **peptideperm;
    static int **reactionperm;

    static void initialize();
};

#endif
