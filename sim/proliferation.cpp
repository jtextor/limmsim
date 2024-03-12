// helper function for somatic hypermutation
int mutate_bitwise( int bitstring, int thread_nr ){
  int ret = bitstring; 
  for( int i = 0 ; i < Settings::nbitstr ; i ++ ){
  	if( COIN( Settings::p_b_rec_mutate ) ){
		ret = ret ^ ( 1 << i );
	}
  }
  return ret;
}

/** proliferation of b and t cells */ 

// contrary to the CS model, proliferation is based on the logistic equation:
// dN/dt = v_prolif N ( 1 - N / capacity ) 

void b_t_proliferate( GridPoint * gp, int thread_nr ){
  ThreadLocals::bcell[thread_nr].MHCIIpep = 0;
  ThreadLocals::bcell[thread_nr].ag = 0;
  ThreadLocals::bcell[thread_nr].dupcycles = 0;

  ThreadLocals::tcell[thread_nr].state = TCell::DUPLICATING;
  ThreadLocals::tcell[thread_nr].dupcycles = 0;
  ThreadLocals::tcell[thread_nr].age = 0;

  double p_prolif; bool birth;
  int neugeborene, memcells;
  BCell * bc; TCell * tc;

  p_prolif = Settings::v_prolif_b * 
    (1.0-(double) gp->b.count_elements() / (double) Settings::capacity);
  birth = (p_prolif >= 0); 
  if( p_prolif < 0 ) p_prolif = -p_prolif;

  for( int i = gp->b.length()-1 ; i >=0 ; i -- ){
    bc = gp->b.item(i)->item; 
    if( bc->state == BCell::DUPLICATING ){
      neugeborene = Sim::monte_carlo( gp->b.item(i)->count, p_prolif, Settings::th_reaction, thread_nr );
      if( neugeborene > 0 ){  
        if( birth ){
			ThreadLocals::bcell[thread_nr].isMem = bc-> isMem;
			if( bc-> dupcycles < Settings::prolif_cycles_b - 1 ){
				// keep track of # of duplications this clone has gone through. (unclear how this would be done in vivo?)
				ThreadLocals::bcell[thread_nr].dupcycles = bc-> dupcycles + 1;
				ThreadLocals::bcell[thread_nr].state = BCell::DUPLICATING;
				// During mitosis, B Cell receptors may mutate. (somatic hypermutation) 
				for( int k = 0 ; k < 2*neugeborene ; k ++ ){
					ThreadLocals::bcell[thread_nr].receptor = mutate_bitwise(bc -> receptor, thread_nr);
					gp->b.add( grid->b.store( &ThreadLocals::bcell[thread_nr], 1 ) , 1 );
	    		}
          	} else {
				// last duplication cycle. B Cells now mature into Plasma and Memory Cells. 
				ThreadLocals::bcell[thread_nr].receptor = bc -> receptor;
				ThreadLocals::bcell[thread_nr].dupcycles = 0;
				memcells = 0;
				for( int j = 0 ; j < neugeborene ; j ++ ){
					if( COIN( Settings::p_b_become_mem ) )
					memcells ++; 
				}
				ThreadLocals::bcell[thread_nr].state = BCell::NAIVE;
				if( neugeborene - memcells > 0 ){
					ThreadLocals::bcell[thread_nr].isMem = false;
					gp->b.add( grid->b.store( &ThreadLocals::bcell[thread_nr], neugeborene - memcells ) , neugeborene - memcells );
				}
				if( memcells > 0 ){
					ThreadLocals::bcell[thread_nr].isMem = true;
					gp->b.add( grid->b.store( &ThreadLocals::bcell[thread_nr], memcells ), memcells );
				}
				ThreadLocals::bcell[thread_nr].isMem = false;
				ThreadLocals::bcell[thread_nr].state = BCell::PLASMA;
				gp->b.add( grid->b.store( &ThreadLocals::bcell[thread_nr], neugeborene ) , neugeborene );
				}
				gp->b.remove( bc , neugeborene );
				grid->b.remove( bc , neugeborene );
			} else {
				gp->b.remove( bc , neugeborene );
				grid->b.remove( bc , neugeborene );
			}
		}
		}
	}

  p_prolif = Settings::v_prolif_t *
    (1.0-(double) gp->t.count_elements() / (double) Settings::capacity);
  birth = (p_prolif >= 0); 
  if( p_prolif < 0 ) p_prolif = -p_prolif;

  for( int i = gp->t.length()-1 ; i >=0 ; i -- ){
    tc = gp->t.item(i)->item; 
    if( tc->state == TCell::DUPLICATING ){
      neugeborene = Sim::monte_carlo( gp->t.item(i)->count, p_prolif, Settings::th_reaction, thread_nr );
      if( neugeborene > 0 ){
        if( birth ){
          ThreadLocals::tcell[thread_nr].cd4 = tc-> cd4;
          if( tc-> dupcycles < Settings::prolif_cycles_t - 1 ){
            ThreadLocals::tcell[thread_nr].dupcycles = tc-> dupcycles + 1;
            ThreadLocals::tcell[thread_nr].state = TCell::DUPLICATING;
          } else  {
            ThreadLocals::tcell[thread_nr].dupcycles = 0;
            ThreadLocals::tcell[thread_nr].state = TCell::EFFECTOR;
          }
          gp->t.remove( tc , neugeborene );
          grid->t.remove( tc , neugeborene );
          gp->t.add( grid->t.store( &ThreadLocals::tcell[thread_nr], 2*neugeborene ) , 2*neugeborene );
        } else {
          gp->t.remove( tc , neugeborene );
          grid->t.remove( tc , neugeborene );
        }
      }
    }
  }
}

void ag_proliferate( GridPoint * gp, int thread_nr ){
  Antigen * ag; int neugeborene; double p_prolif; bool birth;

  p_prolif = Settings::v_prolif_ag *
    (1.0-(double) gp->ag.count_elements() / (double) Settings::capacity);
  birth = (p_prolif >= 0); 
  if( p_prolif < 0 ) p_prolif = -p_prolif;

  for( int i = gp->ag.length()-1 ; i >=0 ; i -- ){ 
    neugeborene = Sim::monte_carlo( gp->ag.item(i)->count, p_prolif, Settings::th_reaction, thread_nr );
    if( neugeborene > 0 ){
	  ag = gp->ag.item(i)->item;
      if( birth ){
        gp->ag.add( grid->ag.store( ag, neugeborene ) , neugeborene );
      }
      else{
        gp->ag.remove( ag, neugeborene );
        grid->ag.remove( ag, neugeborene );
      }
    }
  }
}
