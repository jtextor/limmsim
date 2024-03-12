#include "sim.h"

extern Grid * grid;

void diffuse_point( GridPoint * gp, listPtr loff, double diff_speed, int thread_nr ){
  EntityList<Entity> * l; int move;
  l = &(gp->*loff);
  // approximate diffusion for large entity numbers 
  while( l->length() > 0 ){
    if( l->item(0)->count >= Settings::th_diffusion ){
      move = (int) ((double)l->item(0)->count*(1-diff_speed));
      if( move > 0 )
        ((gp->neighbours[0])->*loff).add( l->item(0)->item, move );
      move = l->item(0)->count-move;
      if( move / 6 > 0 ){
        for( int j = 1 ; j < 7 ; j ++ ){
          ((gp->neighbours[j])->*loff).add( l->item(0)->item, move / 6 );
        }
        move = move % 6;
      }
      if( move > 0 ){
        for( int i = 0 ; i < move ; i ++ )
          ((gp->neighbours[1+RANDINT(6)])->*loff).add( l->item(0)->item );
      }
    }
    // perform exact diffusion for smaller entity numbers 
    else { 
      for( int i = 0 ; i < l->item(0)->count ; i ++ )
        if( COIN(diff_speed) ){
          ((gp->neighbours[1+RANDINT(6)])->*loff).add( l->item(0)->item );
        }
        else{
          ((gp->neighbours[0])->*loff).add( l->item(0)->item );
        }
    }
    l->remove(l->item(0)->item,l->item(0)->count);
  }

}
