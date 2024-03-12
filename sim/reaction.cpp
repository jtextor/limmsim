#include "sim.h"

extern Grid * grid;

#include "phagocytosis.cpp"
#include "proliferation.cpp"
#include "costimulation.cpp"
#include "response.cpp"

void react_point( GridPoint * gp, int thread_nr ){

  Sim::permutate( ThreadLocals::reactionperm[thread_nr], Settings::nr_reactions, thread_nr );

  for( int i = 0 ; i < Settings::nr_reactions ; i ++ ){
    switch( ThreadLocals::reactionperm[thread_nr][i] ){
      case 0: entity_death( gp, thread_nr ); break;
      case 1: lymphocyte_insertion( gp, thread_nr ); break;
      case 2: ag_ab_binding( gp, thread_nr ); break;
      case 3: b_t_proliferate( gp, thread_nr ); break;
      case 4: ag_proliferate( gp, thread_nr ); break;
      case 5: b_t_stimulate( gp, thread_nr ); break;
      case 6: ma_t_stimulate( gp, thread_nr ); break;
      case 7: apc_remove_ag( gp, thread_nr ); break;
      case 8: b_ag_processing( gp, thread_nr ); break;
      case 9: ma_ag_processing( gp, thread_nr ); break;
      case 10: b_ag_phagocytosis( gp, thread_nr ); break;
      case 11: ma_ic_phagocytosis( gp, thread_nr ); break;
      case 12: ma_ag_phagocytosis( gp, thread_nr ); break;
      case 13: secrete_ab( gp, thread_nr ); break; 
      case 14: t_eff_become_mem( gp, thread_nr ); break; 
      case 15: t_eff_become_naive( gp, thread_nr ); break;
      case 16: ma_become_mem( gp, thread_nr ); break;
    }
  }
}

/** death of cells and disintegration of molecules */ 
void entity_death_list( EntityList<Entity> * l, hashPtr hoff, double tau, int thread_nr ){
  int tote;
  Entity * e;
  for( int i = l->length()-1 ; i >=0 ; i -- ){
    e = l->item(i)->item; 
    tote = Sim::monte_carlo( l->item(i)->count, 1.0/tau, Settings::th_reaction, thread_nr );

    if( tote > 0 ){
      l->remove( e , tote );
      (grid->*hoff).remove( e , tote );
    }

  }
}

/** death of b-cells (it is a little more involved) */ 
void bcell_death( GridPoint * gp, int thread_nr ){
  int tote;
  BCell * bc;
  for( int i = gp->b.length()-1 ; i >=0 ; i -- ){
    bc = gp->b.item(i)->item; 
    // Do not kill cells which are currently duplicating.
    // Also do not kill bcells that have just internalized ag, since it
    // would result in the antigen not being freed.
    if( bc->state == BCell::INTERNALIZED || bc->state == BCell::DUPLICATING ) continue;

		tote = 0;
		double tau = 1.0/Settings::tau_b; // non-memory, non-plasma cells
		if( bc->state == BCell::PLASMA ) tau = 1.0/Settings::tau_plasma; // plasma cells 
		else if( bc->isMemory() ) tau = 1.0/Settings::tau_memory; // memory cells

	 tote = Sim::monte_carlo( gp->b.item(i)->count, tau, Settings::th_reaction, thread_nr );
    if( tote > 0 ){
      gp->b.remove( bc, tote );
      grid->b.remove( bc, tote );
    }
  }
}

/** death of t-cells (it is a little more involved) */ 
void tcell_death( GridPoint * gp, int thread_nr ){
  int tote;
  TCell * tc;
  for( int i = gp->t.length()-1 ; i >=0 ; i -- ){
	tc = gp->t.item(i)->item; 
    // Do not kill cells which are currently duplicating.
    if( tc->state == TCell::DUPLICATING ) continue;

	 double tau = 1.0/Settings::tau_t; // non-effector or memory T cells
	 if( tc->state == TCell::EFFECTOR ) tau = 1.0/Settings::tau_t_eff;
	 else if( tc->state == TCell::MEMORY ) tau = 1.0/Settings::tau_memory;
    tote = Sim::monte_carlo( gp->t.item(i)->count, tau, Settings::th_reaction, thread_nr );

    if( tote > 0 ){
		gp->t.remove( tc, tote );
		grid->t.remove( tc, tote );
    }
  }
}

void entity_death( GridPoint * gp, int thread_nr ){
  tcell_death( gp, thread_nr );
  bcell_death( gp, thread_nr );
  entity_death_list( (EntityList<Entity> *) &gp->ab, (EntityHashtable<Entity> Grid::*) &Grid::ab, Settings::tau_ab, thread_nr );
  //antigen currently does not die, it is only killed by the IS.
  //entity_death_list( &gp->ag, &Grid::ag, Settings::tau_ag, thread_nr );
  entity_death_list( (EntityList<Entity> *) &gp->ic, (EntityHashtable<Entity> Grid::*) &Grid::ic, Settings::tau_ic, thread_nr );
  // macrophages lifecycle is not modelled
}

/** appearance of b-cells (from bone marrow) and t-cells (from thymus) */ 
void lymphocyte_insertion( GridPoint * gp, int thread_nr ){
  BCell * bc2; TCell * tc2;
  ThreadLocals::bcell[thread_nr].state = BCell::NAIVE;
  ThreadLocals::bcell[thread_nr].ag = 0;
  ThreadLocals::bcell[thread_nr].MHCIIpep = 0;
  ThreadLocals::bcell[thread_nr].dupcycles = 0;
  ThreadLocals::bcell[thread_nr].isMem = false;
  ThreadLocals::tcell[thread_nr].state = TCell::NAIVE;
  ThreadLocals::tcell[thread_nr].dupcycles = 0;
  ThreadLocals::tcell[thread_nr].age = 0;
  // create new b cells.
  for( int i = 0 ; i < floor( Settings::new_bcells ) ; i ++ ){
    ThreadLocals::bcell[thread_nr].receptor = RANDINT( 1<<Settings::nbitstr );
    bc2 = grid->b.store( &ThreadLocals::bcell[thread_nr] );
    gp->b.add( bc2 );
  }
  if( COIN( Settings::new_bcells - floor( Settings::new_bcells ) ) ){
    ThreadLocals::bcell[thread_nr].receptor = RANDINT( 1<<Settings::nbitstr );
    bc2 = grid->b.store( &ThreadLocals::bcell[thread_nr] );
    gp->b.add( bc2 ) ;
  }
  // create new t cells.
  for( int i = 0 ; i < floor( Settings::new_tcells ) ; i ++ ){
    ThreadLocals::tcell[thread_nr].cd4 = RANDINT( 1<<(Settings::nbitstr>>1) );
    tc2 = grid->t.store( &ThreadLocals::tcell[thread_nr] );
    gp->t.add( tc2 );
  }
  if( COIN( Settings::new_tcells - floor( Settings::new_tcells ) ) ){
    ThreadLocals::tcell[thread_nr].cd4 = RANDINT( 1<<(Settings::nbitstr>>1) );
    tc2 = grid->t.store( &ThreadLocals::tcell[thread_nr] );
    gp->t.add( tc2 );
  }
}

/** State switch of effector T cells to memory cells **/
void t_eff_become_mem( GridPoint * gp, int thread_nr ){
  ThreadLocals::tcell[thread_nr].state = TCell::MEMORY;
  ThreadLocals::tcell[thread_nr].dupcycles = 0;
  ThreadLocals::tcell[thread_nr].age = 0;

  TCell * tc;
  int nr_new_memory_cells;
  for( int i = gp->t.length()-1 ; i >=0 ; i -- ){
		tc = gp->t.item(i)->item; 
		if( tc->state == TCell::EFFECTOR ){
			nr_new_memory_cells = Sim::monte_carlo( gp->t.item(i)->count, 
				Settings::p_t_eff_become_mem, Settings::th_reaction, thread_nr );
			if( nr_new_memory_cells > 0 ){
				ThreadLocals::tcell[thread_nr].cd4 = tc->cd4;
				gp->t.remove( tc , nr_new_memory_cells );
				grid->t.remove( tc , nr_new_memory_cells );
				gp->t.add( grid->t.store( &ThreadLocals::tcell[thread_nr], nr_new_memory_cells ) );
			}
		}
	}
}

/** State switch of effector T cells to naive cells **/
void t_eff_become_naive( GridPoint * gp, int thread_nr ){
  ThreadLocals::tcell[thread_nr].state = TCell::NAIVE;
  ThreadLocals::tcell[thread_nr].dupcycles = 0;
  ThreadLocals::tcell[thread_nr].age = 0;

  TCell * tc;
  int nr_new_naive_cells;
	for( int i = gp->t.length()-1 ; i >=0 ; i -- ){
		tc = gp->t.item(i)->item; 
		if( tc->state == TCell::EFFECTOR ){
			nr_new_naive_cells = Sim::monte_carlo( gp->t.item(i)->count, 
				Settings::p_t_eff_become_naive, Settings::th_reaction, thread_nr );
			if( nr_new_naive_cells > 0 ){
				ThreadLocals::tcell[thread_nr].cd4 = tc->cd4;
				gp->t.remove( tc , nr_new_naive_cells );
				grid->t.remove( tc , nr_new_naive_cells );
				gp->t.add( grid->t.store( &ThreadLocals::tcell[thread_nr], nr_new_naive_cells ) );
			}
		}
	}
}

/** State switch of macrophages to memory cells **/
void ma_become_mem( GridPoint * gp, int thread_nr ){
  ThreadLocals::ma[thread_nr].state = Macrophage::MEMORY;
  ThreadLocals::ma[thread_nr].ag = 0;

  Macrophage * ma;
  int nr_new_mem_cells;
  for( int i = gp->ma.length()-1 ; i >=0 ; i -- ){
    	ma = gp->ma.item(i)->item; 
		if( ma->state == Macrophage::PRES_II ){
			nr_new_mem_cells = Sim::monte_carlo( gp->ma.item(i)->count, 
				Settings::p_ma_become_mem, Settings::th_reaction, thread_nr );
			if( nr_new_mem_cells > 0 ){
				ThreadLocals::ma[thread_nr].MHCIIpep = ma->MHCIIpep;
				gp->ma.remove( ma, nr_new_mem_cells );
				grid->ma.remove( ma, nr_new_mem_cells );
				gp->ma.add( grid->ma.store( &ThreadLocals::ma[thread_nr], nr_new_mem_cells ) );
			}
		}
 	}
}

