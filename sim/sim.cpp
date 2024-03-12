#include "sim.h"

extern Grid * grid;

int mod, blockofthree; 

pthread_t * threads;
int * message;

void housekeeping(){
  // nothing to do yet
}

void sim(){

  threads = new pthread_t[Settings::threads];
  message = new int[Settings::threads];

  blockofthree = Settings::gridsize - Settings::gridsize % 3;
  for( mod = 0 ; mod < 3 ; mod ++ ){
    for( int b = 0 ; b < Settings::threads ; b ++ ){ 
        message[b] = b;
        pthread_create( &threads[b], 0, sim_block, (void *) &message[b]);
    }
    for( int b = 0 ; b < Settings::threads ; b ++ ){ 
        pthread_join( threads[b], 0 );
    }
  }

  for( int l = blockofthree ; l < Settings::gridsize ; l ++ ){
      // perform simulation on line l 
      for( int i = Settings::gridsize * l ; i < Settings::gridsize * (l+1) ; i ++ ){
          sim_point( i, 0 ); 
      }
  }
  
  /** perform tasks not directly related to this simulation step */
  housekeeping();

  /** simstep done. switch grids **/ 
  grid->gridpoints = grid->gridpoints == grid->gridpoints1 ? grid->gridpoints2 : grid->gridpoints1;

  delete[] threads;
  delete[] message;
}

void *sim_block( void *n ){
  int b = *((int*)n);

  for( int l = mod + 3*b ; l < blockofthree ; l += 3*Settings::threads ){
      // perform simulation on line l 
      for( int i = Settings::gridsize * l ; i < Settings::gridsize * (l+1) ; i ++ ){
          sim_point( i, b ); 
      }
  }

  return 0;
}

void sim_point( int i, int b ){
  /** LIST SHUFFLIN' */
  grid->gridpoints[i].t.shuffle( b );
  grid->gridpoints[i].b.shuffle( b );
  grid->gridpoints[i].ag.shuffle( b );
  grid->gridpoints[i].ab.shuffle( b );

  /** REACTION */
  react_point( &grid->gridpoints[i], b );

  /** DIFFUSION */ 
  // diffuse macrophages
  diffuse_point(&grid->gridpoints[i], (EntityList<Entity> GridPoint::*) &GridPoint::ma, Settings::vdiff_ma, b );
  // diffuse t-cells
  diffuse_point(&grid->gridpoints[i], (EntityList<Entity> GridPoint::*) &GridPoint::t, Settings::vdiff_t, b );
  // diffuse b-cells
  diffuse_point(&grid->gridpoints[i], (EntityList<Entity> GridPoint::*) &GridPoint::b, Settings::vdiff_b, b );
  // diffuse antigen
  diffuse_point(&grid->gridpoints[i], (EntityList<Entity> GridPoint::*) &GridPoint::ag, Settings::vdiff_ag, b );
  // diffuse antibody
  diffuse_point(&grid->gridpoints[i], (EntityList<Entity> GridPoint::*) &GridPoint::ab, Settings::vdiff_ab, b );
  // diffuse immune complexes
  diffuse_point(&grid->gridpoints[i], (EntityList<Entity> GridPoint::*) &GridPoint::ic, Settings::vdiff_ic, b );
}
