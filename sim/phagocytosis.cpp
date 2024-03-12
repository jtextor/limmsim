/** Phagocytosis by Macrophages */ 

// phagocytosis 
void ma_ag_phagocytosis( GridPoint * gp, int thread_nr ){
  Antigen * ag; Macrophage * ma; int interactions, max_interactions; bool done; double saturation, inter_p;
  ThreadLocals::ma[thread_nr].MHCIIpep = 0;
  ThreadLocals::ma[thread_nr].state = Macrophage::INTERNALIZED;

  for( int i = gp->ma.length()-1 ; i >= 0 ; i -- ){
    ma = gp->ma.item(i)->item; 
    if( ma->state == Macrophage::ACTIVE ){
      done = false;
      for( int j = gp->ag.length()-1 ; !done && j >= 0 ; j -- ){
        ag = gp->ag.item(j)->item;
        ThreadLocals::ma[thread_nr].ag = ag;

        max_interactions = MIN(gp->ag.item(j)->count,gp->ma.item(i)->count);
        saturation = (double)MAX(gp->ag.item(j)->count,gp->ma.item(i)->count)/(double)Settings::capacity;

        inter_p = Settings::p_ma_ag * saturation;

		interactions = Sim::monte_carlo( max_interactions, inter_p, Settings::th_reaction, thread_nr );
        if( interactions == max_interactions ) done = true;

        if( interactions > 0 ){
          gp->ag.remove( ag, interactions );
          // antigen is not removed from grid now, but when the macrophage has finished processing it. 
          gp->ma.remove( ma, interactions );
          grid->ma.remove( ma, interactions );
          gp->ma.add( grid->ma.store( &ThreadLocals::ma[thread_nr], interactions ), interactions );
        }
      }
    }
  }
}

// ag processing
void ma_ag_processing( GridPoint * gp, int thread_nr ){
  int aktpep,aktpep_half,boundpep; Macrophage * ma;
  ThreadLocals::ma[thread_nr].MHCIIpep = 0;
  ThreadLocals::ma[thread_nr].state = Macrophage::PRES_II;
  ThreadLocals::ma[thread_nr].ag = 0;

  for( int i = gp->ma.length()-1 ; i >= 0 ; i -- ){
    ma = gp->ma.item(i)->item; 
    if( ma->state == Macrophage::INTERNALIZED ){

      for( int j = gp->ma.item(i)->count-1 ; j >= 0  ; j -- ){
        // consider ag peptides in random order
        Sim::permutate( ThreadLocals::peptideperm[thread_nr], Settings::ag_peptides, thread_nr );

        boundpep = -1;

        for( int p = 0 ; p < Settings::ag_peptides ; p ++ ){
          aktpep = ma->ag->peptides[ThreadLocals::peptideperm[thread_nr][p]];

          if( RANDINT(2) ){ // left half, then right half;
            aktpep_half = aktpep >> (Settings::nbitstr >> 1);
            if( COIN( Affinity::mhc_affinity( aktpep_half, Settings::mhc ) ) ){
              boundpep = aktpep_half; break;
            }
            aktpep_half = aktpep & ( ( 1 << ( Settings::nbitstr >>1 ) ) - 1 );
            if( COIN( Affinity::mhc_affinity( aktpep_half, Settings::mhc ) ) ){
              boundpep = aktpep_half; break;
            }
          }
          else{ // right half, then left half;
            aktpep_half = aktpep & ( ( 1 << ( Settings::nbitstr >>1 ) ) - 1 );
            if( COIN( Affinity::mhc_affinity( aktpep_half, Settings::mhc ) ) ){
              boundpep = aktpep_half; break;
            }
            aktpep_half = aktpep >> (Settings::nbitstr >> 1);
            if( COIN( Affinity::mhc_affinity( aktpep_half, Settings::mhc ) ) ){
              boundpep = aktpep_half; break;
            }
          }
        }

        // antigen is dead now. 
        grid->ag.remove( ma->ag );

        gp->ma.remove( ma );
        grid->ma.remove( ma );

        // ag binding successful.
        if( boundpep > -1 ){
          ThreadLocals::ma[thread_nr].MHCIIpep = boundpep;
          ThreadLocals::ma[thread_nr].state = Macrophage::PRES_II;
        }
        else{
          ThreadLocals::ma[thread_nr].MHCIIpep = 0;
          ThreadLocals::ma[thread_nr].state = Macrophage::ACTIVE;
        }
        gp->ma.add( grid->ma.store( &ThreadLocals::ma[thread_nr] ) );
      }
    }
  }
}

/* removal of presented ag from b cells and ma */
void apc_remove_ag( GridPoint * gp, int thread_nr ){
  Macrophage * ma;
  ThreadLocals::ma[thread_nr].MHCIIpep = 0;
  ThreadLocals::ma[thread_nr].state = Macrophage::ACTIVE;
  ThreadLocals::ma[thread_nr].ag = 0;
  BCell * bc;
  ThreadLocals::bcell[thread_nr].MHCIIpep = 0;
  ThreadLocals::bcell[thread_nr].state = BCell::NAIVE;
  ThreadLocals::bcell[thread_nr].ag = 0;
  ThreadLocals::bcell[thread_nr].dupcycles = 0;
  for( int i = gp->ma.length()-1 ; i >= 0 ; i -- ){
    ma = gp->ma.item(i)->item; 
    if( ma->state == Macrophage::PRES_II ){
	  int nr_remove = Sim::monte_carlo( gp->ma.item(i)->count, Settings::p_remove_ag, Settings::th_reaction, thread_nr );
      if( nr_remove > 0 ){
        gp->ma.remove( ma, nr_remove );
        grid->ma.remove( ma, nr_remove );
        gp->ma.add( grid->ma.store( &ThreadLocals::ma[thread_nr] ), nr_remove );
      }
    }
  }
  for( int i = gp->b.length()-1 ; i >= 0 ; i -- ){
    bc = gp->b.item(i)->item; 
    ThreadLocals::bcell[thread_nr].receptor = bc->receptor;
    ThreadLocals::bcell[thread_nr].isMem = bc->isMem;
    if( bc->state == BCell::PRES_II ){
	  int nr_remove = Sim::monte_carlo( gp->b.item(i)->count, Settings::p_remove_ag, Settings::th_reaction, thread_nr );
      if( nr_remove > 0 ){
        gp->b.remove( bc, nr_remove );
        grid->b.remove( bc, nr_remove );
        gp->b.add( grid->b.store( &ThreadLocals::bcell[thread_nr], nr_remove ) );
      }
    }
  }
}

/** Phagocytosis by B-Cells */ 

// phagocytosis 
void b_ag_phagocytosis( GridPoint * gp, int thread_nr ){
  Antigen * ag; BCell * bc; int max_interactions, interactions; bool done; double saturation, inter_p;

  ThreadLocals::bcell[thread_nr].MHCIIpep = 0;
  ThreadLocals::bcell[thread_nr].state = BCell::INTERNALIZED;
  ThreadLocals::bcell[thread_nr].dupcycles = 0;

  for( int i = gp->b.length()-1 ; i >= 0 ; i -- ){
    bc = gp->b.item(i)->item; 
    ThreadLocals::bcell[thread_nr].receptor = bc->receptor;
    ThreadLocals::bcell[thread_nr].isMem = bc->isMem;

    if( bc->state == BCell::NAIVE ){
      done = false;
      for( int j = gp->ag.length()-1 ; !done && j >= 0 ; j -- ){
        ag = gp->ag.item(j)->item; 
        ThreadLocals::bcell[thread_nr].ag = ag;

        interactions = 0; 

		max_interactions = MIN(gp->b.item(i)->count,gp->ag.item(j)->count);
		saturation = ((double)MAX(gp->b.item(i)->count,gp->ag.item(j)->count)) / (double) Settings::capacity;
		inter_p = Settings::p_b_ag * saturation * Affinity::affinity( bc->receptor, ag->epitopes[0] );
		interactions = Sim::monte_carlo( max_interactions, inter_p, Settings::th_reaction, thread_nr );
        if( interactions == max_interactions ) done = true;

        if( interactions > 0 ){
          gp->ag.remove( ag, interactions );
          // antigen is not removed from grid now, but when the b-cell has finished processing it. 
          gp->b.remove( bc, interactions );
          grid->b.remove( bc, interactions );
          gp->b.add( grid->b.store( &ThreadLocals::bcell[thread_nr], interactions ), interactions );
        }
      }
    }
  }
}

// ag processing
void b_ag_processing( GridPoint * gp, int thread_nr ){
  int aktpep,aktpep_half,boundpep; BCell * bc;
  ThreadLocals::bcell[thread_nr].ag = 0;
  ThreadLocals::bcell[thread_nr].dupcycles = 0;

  for( int i = gp->b.length()-1 ; i >= 0 ; i -- ){
    bc = gp->b.item(i)->item; 
    ThreadLocals::bcell[thread_nr].receptor = bc->receptor;
    ThreadLocals::bcell[thread_nr].isMem = bc->isMem;
    if( bc->state == BCell::INTERNALIZED ){

      for( int j = gp->b.item(i)->count-1 ; j >= 0  ; j -- ){
        // consider ag peptides in random order
        Sim::permutate( ThreadLocals::peptideperm[thread_nr], Settings::ag_peptides, thread_nr );

        boundpep = -1;

        for( int p = 0 ; p < Settings::ag_peptides ; p ++ ){
          aktpep = ((BCell *)bc)->ag->peptides[ThreadLocals::peptideperm[thread_nr][p]];

          if( RANDINT(2) ){ // left half, then right half;
            aktpep_half = aktpep >> (Settings::nbitstr >> 1);
            if( COIN( Affinity::mhc_affinity( aktpep_half, Settings::mhc ) ) ){
              boundpep = aktpep_half; break;
            }
            aktpep_half = aktpep & ( ( 1 << ( Settings::nbitstr >>1 ) ) - 1 );
            if( COIN( Affinity::mhc_affinity( aktpep_half, Settings::mhc ) ) ){
              boundpep = aktpep_half; break;
            }
          }
          else{ // right half, then left half;
            aktpep_half = aktpep & ( ( 1 << ( Settings::nbitstr >>1 ) ) - 1 );
            if( COIN( Affinity::mhc_affinity( aktpep_half, Settings::mhc ) ) ){
              boundpep = aktpep_half; break;
            }
            aktpep_half = aktpep >> (Settings::nbitstr >> 1);
            if( COIN( Affinity::mhc_affinity( aktpep_half, Settings::mhc ) ) ){
              boundpep = aktpep_half; break;
            }
          }
        }

        // antigen is dead now. 
        grid->ag.remove( bc->ag );

        gp->b.remove( bc );
        grid->b.remove( bc );

        // switch B cell to antigen presenting state if the binding was
        // successful.
        if( boundpep > -1 ){
          ThreadLocals::bcell[thread_nr].MHCIIpep = boundpep;
          ThreadLocals::bcell[thread_nr].state = BCell::PRES_II;
        }
        else{
          ThreadLocals::bcell[thread_nr].MHCIIpep = 0;
          ThreadLocals::bcell[thread_nr].state = BCell::NAIVE;
        }
        gp->b.add( grid->b.store( &ThreadLocals::bcell[thread_nr] ) );
      }
    }
  }
}

