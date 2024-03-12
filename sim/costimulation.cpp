
// antigen presentation to t cells by macrophages
void ma_t_stimulate( GridPoint * gp, int thread_nr ){
  TCell * tc; Macrophage * ma; bool done; int interactions, max_interactions; double inter_p, saturation;

  ThreadLocals::tcell[thread_nr].state = TCell::DUPLICATING;
  ThreadLocals::tcell[thread_nr].dupcycles = 0;
  ThreadLocals::tcell[thread_nr].age = 0;

  // iterate over pairs of macrophages (antigen presenting short-lived or memory)
  // and naive T cells
  for( int i = gp->ma.length()-1 ; i >= 0 ; i -- ){
    ma = gp->ma.item(i)->item; 
    if( ma->state == Macrophage::PRES_II || 
		ma->state == Macrophage::MEMORY ){
      done = false;
      for( int j = gp->t.length()-1 ; !done && j >= 0 ; j -- ){
        tc = gp->t.item(j)->item;
        if( tc->state == TCell::NAIVE ){
			ThreadLocals::tcell[thread_nr].cd4 = tc->cd4;
			max_interactions = MIN(gp->ma.item(i)->count,gp->t.item(j)->count);
			saturation = ((double)MAX(gp->ma.item(i)->count,gp->t.item(j)->count)) / (double) Settings::capacity;
			inter_p = Settings::p_ma_t * saturation * Affinity::cd4_affinity( tc->cd4, ma->MHCIIpep );
			interactions = Sim::monte_carlo( max_interactions, inter_p, Settings::th_reaction, thread_nr );
			if( interactions == max_interactions ) done = true;

			if( interactions > 0 ){
				gp->t.remove( tc, interactions );
				grid->t.remove( tc, interactions );
				gp->t.add( grid->t.store( &ThreadLocals::tcell[thread_nr], interactions ), interactions );
			}
        }
      }
    }
  }
}

// antigen presentation to t by b 
void b_t_stimulate( GridPoint * gp, int thread_nr ){
  bool done; BCell * bc; TCell * tc; int interactions, max_interactions; double inter_p, saturation;
  ThreadLocals::tcell[thread_nr].state = TCell::DUPLICATING;
  ThreadLocals::tcell[thread_nr].dupcycles = 0;
  ThreadLocals::tcell[thread_nr].age = 0;

  ThreadLocals::bcell[thread_nr].MHCIIpep = 0;
  ThreadLocals::bcell[thread_nr].state = BCell::DUPLICATING;
  ThreadLocals::bcell[thread_nr].ag = 0;
  ThreadLocals::bcell[thread_nr].dupcycles = 0;

  // Iterate over pairs of B and T cells at current grid point
  for( int i = gp->b.length()-1 ; i >= 0 ; i -- ){
    bc = gp->b.item(i)->item; 
	 // Only consider antigen presenting B cells ... 
    if( bc->state == BCell::PRES_II ){
      ThreadLocals::bcell[thread_nr].receptor = bc->receptor;
      ThreadLocals::bcell[thread_nr].isMem = bc->isMem;
      done = false;
      for( int j = gp->t.length()-1 ; !done && j >= 0 ; j -- ){
        tc = gp->t.item(j)->item;
        // ... versus effector T cells
        if( tc->state == TCell::EFFECTOR ){
			ThreadLocals::tcell[thread_nr].cd4 = tc->cd4;
			max_interactions = MIN(gp->b.item(i)->count,gp->t.item(j)->count);
			saturation = ((double)MAX(gp->b.item(i)->count,gp->t.item(j)->count)) / (double) Settings::capacity;
			inter_p = Settings::p_b_t * saturation * Affinity::cd4_affinity( tc->cd4, bc->MHCIIpep );
			interactions = Sim::monte_carlo( max_interactions, inter_p, Settings::th_reaction, thread_nr );
			if( interactions == max_interactions ) done = true;
			
			if( interactions > 0 ){
			gp->t.remove( tc, interactions );
			grid->t.remove( tc, interactions );
			gp->t.add( grid->t.store( &ThreadLocals::tcell[thread_nr], interactions ), interactions );
			
			gp->b.remove( bc, interactions );
			grid->b.remove( bc, interactions );
			gp->b.add( grid->b.store( &ThreadLocals::bcell[thread_nr], interactions ), interactions );
          }
        }
      }
    }
  }
}

