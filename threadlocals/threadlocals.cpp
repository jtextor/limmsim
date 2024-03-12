#include "threadlocals.h"

pthread_t * ThreadLocals::threads = 0;
int * ThreadLocals::message = 0;
MTRand * ThreadLocals::rng = 0;
BCell * ThreadLocals::bcell = 0;
TCell * ThreadLocals::tcell = 0;
Macrophage * ThreadLocals::ma = 0;
Antibody * ThreadLocals::ab = 0;
ImmuneComplex * ThreadLocals::ic = 0;
int ** ThreadLocals::peptideperm = 0;
int ** ThreadLocals::reactionperm = 0;

void ThreadLocals::initialize(){
  threads = new pthread_t[Settings::threads];
  message = new int[Settings::threads];
  rng = new MTRand[Settings::threads];
  bcell = new BCell[Settings::threads];
  tcell = new TCell[Settings::threads];
  ma = new Macrophage[Settings::threads];
  ab = new Antibody[Settings::threads];
  ic = new ImmuneComplex[Settings::threads];

  peptideperm = new int*[Settings::threads];
  for( int i = 0 ; i < Settings::threads ; i ++ )
    peptideperm[i] = new int[Settings::ag_peptides];
  for( int i = 0 ; i < Settings::threads ; i ++ )
    for( int j = 0 ; j < Settings::ag_peptides; j ++ ) 
      peptideperm[i][j] = j;

  reactionperm = new int*[Settings::threads];
  for( int i = 0 ; i < Settings::threads ; i ++ )
    reactionperm[i] = new int[Settings::nr_reactions];
  for( int i = 0 ; i < Settings::threads ; i ++ )
    for( int j = 0 ; j < Settings::nr_reactions; j ++ ) 
      reactionperm[i][j] = j;
}
