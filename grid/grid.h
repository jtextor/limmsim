#ifndef __GRID_H

#define __GRID_H

#include "../settings/settings.h"
#include "../lists/entityhashtable.h"
#include "../threadlocals/threadlocals.h"
#include "gridpoint.h"
#include <pthread.h>
#include <string>

class Grid
{
  public:
    GridPoint * gridpoints1;
    GridPoint * gridpoints2;
    GridPoint * gridpoints;

    Grid();
    ~Grid();

    // global lists. hold global object counts and the objects themselves to save memory.
    EntityHashtable<Antigen> ag;
    EntityHashtable<Antibody> ab;
    EntityHashtable<ImmuneComplex> ic;
    EntityHashtable<TCell> t;
    EntityHashtable<BCell> b;
    EntityHashtable<Macrophage> ma;

    // IMPORTANT: all static method expect globar var "grid" to exist and point to an 
    // initialized grid instance. I did this because one cannot multithread member methods,
    // only static methods. 

    // generate cells of immune system 
    static void generate_is( );

    // inject antigen to grid 
    static void inject_ag( Antigen * );

    // dump ag concentration as bmp 
    static void dump_bmps( char * );

    // global statistics per cell state
    static std::string ma_stats();
    static std::string b_stats();
    static std::string t_stats();

    // detailed statistics (complete dump of global lists)
    static std::string ab_details();

    int size; 

  private:
    static void * generate_is_mod( void * );
    static void generate_is_line( int, int m );

    // polymorph bmp dumping 
    static void dump_bmp( listPtr, char *  );


    // color arrays for bmp output
    unsigned char * red;
    unsigned char * green;
    unsigned char * blue;

};

typedef EntityHashtable<Entity> Grid::* hashPtr;

#endif
