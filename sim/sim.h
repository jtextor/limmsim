#ifndef __SIM_H

#define __SIM_H

#include "../settings/settings.h"
#include "../grid/grid.h"
#include "../grid/gridpoint.h"
#include "../mersennetwister/mersennetwister.h"
#include "../entities/bcell.h"
#include "../entities/tcell.h"
#include "../entities/macrophage.h"
#include "../threadlocals/threadlocals.h"
#include "../affinity/affinity.h"
#include <iostream>
#include <pthread.h>

#define MIN(n,m) (n<m?n:m)
#define MAX(n,m) (n<m?m:n)

#define COIN(p) (ThreadLocals::rng[thread_nr].randExc() < (p))
#define RANDINT(n) (ThreadLocals::rng[thread_nr].randInt( n - 1 ))

/** Helper class for common simulation tasks */ 
class Sim{

public:
	/** Shuffle the integer list * ar */ 
	static inline void permutate( int * ar, int n, int thread_nr ){
		int tmp,dst;
		for( int i = 0 ; i < n-1 ; i ++ ){
			dst = RANDINT(n-i);
			tmp = ar[dst];
			ar[dst]=ar[i];
			ar[i] = tmp;
		}
	};

	/** Flip n coins with probability p, return the number of successful flips.
		 If n >= th, do not acually flip the coins, but return the expected value
		 instead. If that is not an integer, flip one coin with "residual"
		 success probability. */ 
	static inline int monte_carlo( int n, double p, int th, int thread_nr ){
		if( n >= th ){
			double d = p*(double)n;
			if( d >= n ) return n;
			// since the expected outcome is not always integer, perform
			// additional coin throw to achieve outcome on average
			return (floor(d)+COIN(d-(double)floor(d))?1:0);
		} else {
			int r = 0;
			for( int i = 0; i < n ; i ++ ) if( COIN(p) ) r++;
			return r;
		}
	};
};

// helper function for the common task of replacing an entity on the grid
// with another one 

// no classes here because of openmp.

// master sim entry function.
void sim();

// general functions to perform between simulation steps (e.g. aging of cells)
void housekeeping();

// sim function for all lines mod m. executed in thread
// from here on, the last param is always the nr of the current thread
// (used to choose the right rng) 
void *sim_block( void * );

// sim function for one gridpoint.
void sim_point( int , int  );

// diffusion for one particle type at one gridpoint
void diffuse_point( GridPoint * , listPtr, double, int  );

// reaction for one particle type at one gridpoint 
void react_point( GridPoint *, int  );

// birth and death of entities
void entity_death( GridPoint *, int );
// insertion of new t and b cells
void lymphocyte_insertion( GridPoint *, int );

// unspecific phagocytosis of ag by macrophages
void ma_ag_phagocytosis( GridPoint *, int );
void ma_ag_processing( GridPoint *, int );

// specific phagocytosis of ag by bcells
void b_ag_phagocytosis( GridPoint *, int );
void b_ag_processing( GridPoint *, int );

// removal of presented ag
void apc_remove_ag( GridPoint *, int );

// t cell help (costimulation)
void ma_t_stimulate( GridPoint *, int );
void b_t_stimulate( GridPoint *, int );

// proliferation of b and t cells 
void b_t_proliferate( GridPoint *, int );
// helper function for somatic hypermutation
int mutate_bitwise( int bitstring, int thread_nr );

// ab secretion by plasma b cells
void secrete_ab( GridPoint *, int );

// ag-ab binding
void ag_ab_binding( GridPoint *, int );

// ic removal
void ma_ic_phagocytosis( GridPoint *, int );

// state transition of effector t cells to memory or back to naive t cells
void t_eff_become_mem( GridPoint *, int );
void t_eff_become_naive( GridPoint *, int );

// state transition of macrophages into memory cells which carry persisting antigen
void ma_become_mem( GridPoint *, int );

#endif
