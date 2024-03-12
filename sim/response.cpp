/** humoral response */ 

void secrete_ab( GridPoint * gp, int thread_nr ){
  ThreadLocals::ab[thread_nr].paratope = 0;

  BCell * bc;

  for( int i = gp->b.length()-1 ; i >=0 ; i -- ){
    bc = gp->b.item(i)->item; 
    if( bc->state == BCell::PLASMA ){
      ThreadLocals::ab[thread_nr].paratope = bc->receptor;
      gp->ab.add( grid->ab.store( &ThreadLocals::ab[thread_nr], 
          Settings::ab_secrete*gp->b.item(i)->count ), 
          Settings::ab_secrete*gp->b.item(i)->count );
    }
  }
}

void ag_ab_binding( GridPoint * gp, int thread_nr ){
  Antigen * ag; Antibody * ab; int interactions, max_interactions; bool done; double saturation, inter_p;
  for( int i = gp->ab.length()-1 ; i >=0 ; i -- ){
    ab = gp->ab.item(i)->item; 
    done = false;
    for( int j = gp->ag.length()-1 ; !done && j >=0 ; j -- ){
		ag = gp->ag.item(j)->item; 
		max_interactions =  MIN(gp->ab.item(i)->count,gp->ag.item(j)->count);
		saturation = (double)MAX(gp->ab.item(i)->count,gp->ag.item(j)->count) / (double) Settings::capacity;
		inter_p = Settings::p_ag_ab * saturation * Affinity::affinity( ag->epitopes[0], ab->paratope );
		interactions = Sim::monte_carlo( max_interactions, inter_p, Settings::th_reaction, thread_nr );
		if( interactions == max_interactions ) done = true;
		if( interactions > 0 ){
			gp->ag.remove( ag, interactions );
			grid->ag.remove( ag, interactions );
			gp->ab.remove( ab, interactions );
			grid->ab.remove( ab, interactions );
			gp->ic.add( grid->ic.store( &ThreadLocals::ic[thread_nr], interactions ), interactions );
		}
	}
  }
}

void ma_ic_phagocytosis(  GridPoint * gp, int thread_nr ){
  Macrophage * ma; ImmuneComplex * ic; int interactions, max_interactions; bool done; double saturation, inter_p;

  for( int i = gp->ma.length()-1 ; i >=0 ; i -- ){
    ma = gp->ma.item(i)->item; 
    done = false;
    for( int j = gp->ic.length()-1 ; !done && j >=0 ; j -- ){
      ic = gp->ic.item(j)->item; 

      interactions = 0; 

      max_interactions = MIN(gp->ic.item(j)->count,gp->ma.item(i)->count);
      saturation = (double)MAX(gp->ic.item(j)->count,gp->ma.item(i)->count)/(double)Settings::capacity;

      inter_p = Settings::p_ma_ic * saturation;

      if( max_interactions < Settings::th_reaction ){
        for( int k = 0 ; k < max_interactions ; k ++ )
          if( COIN( inter_p ) ) interactions ++; 
      } else {
        interactions = (int) floor( ((double)max_interactions)*inter_p );
        if( COIN( ((double)max_interactions)*inter_p - floor(((double)max_interactions)*inter_p) ) ) interactions ++; 
      }
      if( interactions == max_interactions ) done = true;

      if( interactions > 0 ){
        gp->ic.remove( ic, interactions );
        grid->ic.remove( ic, interactions );
      }
    }
  }
}

