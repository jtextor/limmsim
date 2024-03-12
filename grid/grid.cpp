#include "grid.h"

#include "../entities/macrophage.h"
#include "../bmpio/bmp_io.h"
#include "../entities/tcell.h"
#include "../entities/bcell.h"

#include <iostream>
#include <string>
#include <sstream>

// there is just one grid, so make it global
Grid * grid;

void * Grid::generate_is_mod( void * msg ){
  int m = *((unsigned int *)msg);
  for( int l = m ; l < Settings::gridsize ; l += Settings::threads )
    grid->generate_is_line( l, m );
  return 0;
}

// generate IS on one line (parallel in thread) 
void Grid::generate_is_line( int l, int m ){
  int rec;
  Macrophage mac; TCell tc; BCell bc;

  // 1. generate Macrophages 
  // distribute macrophages uniformly on grid
  Macrophage * mac2 = grid->ma.store( &mac, Settings::init_macrophages*Settings::gridsize );
  for( int i = l*Settings::gridsize ; i < (l+1)*Settings::gridsize ; i ++ ){
    grid->gridpoints[i].ma.add( mac2, Settings::init_macrophages );
  }


  // 2. generate Th-Cells 
  TCell * tc2;
  for( int i = l*Settings::gridsize ; i < (l+1)*Settings::gridsize ; i ++ ){
    for( int j = 0 ; j < Settings::init_tcells ; j ++ ){
      rec = ThreadLocals::rng[m].randInt( ( 1 << (Settings::nbitstr>>1) )-1 );
      tc.cd4 = rec;
      tc2 = grid->t.store( &tc );
      grid->gridpoints[i].t.add( tc2 );
    }
  }

  // 3. generate B-Cells 
  BCell * bc2;
  for( int i = l*Settings::gridsize ; i < (l+1)*Settings::gridsize ; i ++ ){
    for( int j = 0 ; j < Settings::init_bcells ; j ++ ){
      rec = ThreadLocals::rng[m].randInt( ( 1 << Settings::nbitstr )-1 );
      bc.receptor = rec;
      bc2 = grid->b.store( &bc );
      grid->gridpoints[i].b.add( bc2 );
    }
  }
}


// generate entire IS 
void Grid::generate_is(){
  Antigen ag1;

  for( int m = 0 ; m < Settings::threads ; m++ ){ 
    ThreadLocals::message[m] = m;
    pthread_create( &ThreadLocals::threads[m], 0, &(generate_is_mod), (void *) &ThreadLocals::message[m]);
  }
  for( int m = 0 ; m < Settings::threads ; m ++ ){ 
      pthread_join( ThreadLocals::threads[m], 0 );
  }
}

// inject antigen 
void Grid::inject_ag( Antigen * ag1 ){
  Antigen * a2 = grid->ag.store( ag1, Settings::capacity );
  grid->gridpoints[grid->size / 2+Settings::gridsize/2].ag.add( a2, Settings::capacity );
}

Grid::Grid(){
  size = Settings::gridsize * Settings::gridsize;

  gridpoints1 = new GridPoint[size];
  gridpoints2 = new GridPoint[size];

  red = new unsigned char[size];
  green = new unsigned char[size];
  blue = new unsigned char[size];
  for( int i = 0 ; i < size ; i ++ ){
     green[i] = 0;
     blue[i] = 0;
  }

  gridpoints = gridpoints1;

  int x,y,new_x,new_y;
  // precalculate 6-neighbourhood for each grid point 
  for( int i=0; i<size; i++) {
      x = (int)(i / Settings::gridsize);
      y = i % Settings::gridsize;
      gridpoints1[i].neighbours[0] = &gridpoints2[i];
      gridpoints2[i].neighbours[0] = &gridpoints1[i];
      if( i%2 == 1) { 
        /* odd */
        /* 0 */   
        new_x= (x-1+Settings::gridsize)%Settings::gridsize; new_y=y%Settings::gridsize;
        gridpoints1[i].neighbours[1]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[1]=&gridpoints1[new_x * Settings::gridsize + new_y];
        /* 1 */   new_x= (x-1+Settings::gridsize)%Settings::gridsize; new_y=(y+1)%Settings::gridsize;  
        gridpoints1[i].neighbours[2]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[2]=&gridpoints1[new_x * Settings::gridsize + new_y];
        /* 2 */   new_x= x; new_y=(y+1)%Settings::gridsize; 
        gridpoints1[i].neighbours[3]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[3]=&gridpoints1[new_x * Settings::gridsize + new_y];
        /* 3 */   new_x= (x+1)%Settings::gridsize; new_y=(y+1)%Settings::gridsize;
        gridpoints1[i].neighbours[4]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[4]=&gridpoints1[new_x * Settings::gridsize + new_y];
        /* 4 */   new_x= (x+1)%Settings::gridsize; new_y=y%Settings::gridsize;
        gridpoints1[i].neighbours[5]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[5]=&gridpoints1[new_x * Settings::gridsize + new_y];
        /* 5 */   new_x= x; new_y=(y-1+Settings::gridsize)%Settings::gridsize;
        gridpoints1[i].neighbours[6]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[6]=&gridpoints1[new_x * Settings::gridsize + new_y];
      } else { 
        // even 
        // 0   
        new_x= (x-1+Settings::gridsize)%Settings::gridsize; 
        new_y=(y-1+Settings::gridsize)%Settings::gridsize;  
        gridpoints1[i].neighbours[1]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[1]=&gridpoints1[new_x * Settings::gridsize + new_y];
        // 1 
        new_x= (x-1+Settings::gridsize)%Settings::gridsize; new_y=y%Settings::gridsize;
        gridpoints1[i].neighbours[2]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[2]=&gridpoints1[new_x * Settings::gridsize + new_y];
        // 2 
        new_x= x; new_y=(y+1)%Settings::gridsize;
        gridpoints1[i].neighbours[3]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[3]=&gridpoints1[new_x * Settings::gridsize + new_y];
        // 3 
        new_x= (x+1)%Settings::gridsize; new_y=y%Settings::gridsize;
        gridpoints1[i].neighbours[4]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[4]=&gridpoints1[new_x * Settings::gridsize + new_y];
        // 4
        new_x= (x+1)%Settings::gridsize; new_y=(y-1+Settings::gridsize)%Settings::gridsize;  
        gridpoints1[i].neighbours[5]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[5]=&gridpoints1[new_x * Settings::gridsize + new_y];
        // 5 
        new_x= x; new_y=(y-1+Settings::gridsize)%Settings::gridsize;
        gridpoints1[i].neighbours[6]=&gridpoints2[new_x * Settings::gridsize + new_y];
        gridpoints2[i].neighbours[6]=&gridpoints1[new_x * Settings::gridsize + new_y];
      }
  }
}

Grid::~Grid(){
  delete [] gridpoints1;
  delete [] gridpoints2;
  delete [] red;
  delete [] green;
  delete [] blue;
}

void Grid::dump_bmp( listPtr ptr, char * name ){

  EntityList<Entity> * l; int ac;
  for( int i = 0 ; i < grid->size ; i ++ ){
    grid->red[i] = 0;
    l = &(grid->gridpoints[i].*ptr);
    for( int j = 0 ; j < l->length(); j ++ ){
      ac = l->item(j)->count;
      if(grid->red[i] + ac >= 255){
        grid->red[i] = 255;
        break;
      }
      else{ 
        grid->red[i] += ac;
      }
    }
  }

  bmp_24_write ( name, Settings::gridsize, Settings::gridsize, grid->red, grid->green, grid->blue );
}

void Grid::dump_bmps( char * num ){
  char buf[255];

  strncpy( buf, "output/dumps/ag/", 20 ); 
  strncat( buf, num, 10  );
  strncat( buf, ".bmp" , 4  );

  EntityList<Entity> * l; int ac;
  for( int i = 0 ; i < grid->size ; i ++ ){
    grid->red[i] = 0;
    grid->green[i] = 0;
    grid->blue[i] = 0;
    l = (EntityList<Entity> *)&(grid->gridpoints[i].ag);
    for( int j = 0 ; j < l->length(); j ++ ){
      ac = l->item(j)->count;
      if(grid->red[i] + ac >= 255){
        grid->red[i] = 255;
        break;
      }
      else{ 
        grid->red[i] += ac;
      }
    }
    l = (EntityList<Entity> *)&(grid->gridpoints[i].ab);
    for( int j = 0 ; j < l->length(); j ++ ){
      ac = l->item(j)->count;
      if(grid->green[i] + ac >= 255){
        grid->green[i] = 255;
        break;
      }
      else{ 
        grid->green[i] += ac;
      }
    }
    l = (EntityList<Entity> *)&(grid->gridpoints[i].ic);
    for( int j = 0 ; j < l->length(); j ++ ){
      ac = l->item(j)->count;
      if(grid->blue[i] + ac >= 255){
        grid->blue[i] = 255;
        break;
      }
      else{ 
        grid->blue[i] += ac;
      }
    }
  }

  bmp_24_write ( buf, Settings::gridsize, Settings::gridsize, grid->red, grid->green, grid->blue );

  strncpy( buf, "output/dumps/b/", 20 ); 
  strncat( buf, num, 10  );
  strncat( buf, ".bmp" , 4  );

  for( int i = 0 ; i < grid->size ; i ++ ){
    grid->red[i] = 0;
    grid->green[i] = 0;
    grid->blue[i] = 0;
    l = (EntityList<Entity> *)&(grid->gridpoints[i].b);
    for( int j = 0 ; j < l->length(); j ++ ){
      ac = l->item(j)->count;
      if(grid->red[i] + ac >= 255){
        grid->red[i] = 255;
        break;
      }
      else{ 
        grid->red[i] += ac;
      }
    }
    l = (EntityList<Entity> *)&(grid->gridpoints[i].t);
    for( int j = 0 ; j < l->length(); j ++ ){
      ac = l->item(j)->count;
      if(grid->green[i] + ac >= 255){
        grid->green[i] = 255;
        break;
      }
      else{ 
        grid->green[i] += ac;
      }
    }
    l = (EntityList<Entity> *)&(grid->gridpoints[i].ma);
    for( int j = 0 ; j < l->length(); j ++ ){
      ac = l->item(j)->count;
      if(grid->blue[i] + ac >= 255){
        grid->blue[i] = 255;
        break;
      }
      else{ 
        grid->blue[i] += ac;
      }
    }
  }

  bmp_24_write ( buf, Settings::gridsize, Settings::gridsize, grid->red, grid->green, grid->blue );
}

std::string Grid::ma_stats( ){
  std::stringstream rstream;
  int active_count = 0, internalized_count = 0, pres_count = 0, mem_count = 0;
  // active ma
  for( int l = 0 ; l < grid->ma.size ; l ++ ){
    for( int i = 0 ; i < grid->ma.lists[l].length() ; i ++ ){
		switch( ((Macrophage *)grid->ma.lists[l].item(i)->item)->state ){
			case Macrophage::ACTIVE: active_count += grid->ma.lists[l].item(i)->count; break;
			case Macrophage::INTERNALIZED: internalized_count += grid->ma.lists[l].item(i)->count; break;
			case Macrophage::PRES_II: pres_count += grid->ma.lists[l].item(i)->count; break;
			case Macrophage::MEMORY: mem_count += grid->ma.lists[l].item(i)->count; break;
		}
	}
  }
  rstream << active_count << '\t' << internalized_count << '\t' << pres_count << '\t' << mem_count;
  return rstream.str();
}

std::string Grid::b_stats( ){
  std::stringstream rstream;
  int naive_count = 0, internalized_count = 0,
		pres_ii_count = 0, duplicating_count = 0, plasma_count = 0, mem_count =0;
  // active bcells
  for( int l = 0 ; l < grid->b.size ; l ++ ){
    for( int i = 0 ; i < grid->b.lists[l].length() ; i ++ ){
		switch( ((BCell *)grid->b.lists[l].item(i)->item)->state ){
			case BCell::NAIVE:
				naive_count += grid->b.lists[l].item(i)->count; break;
			case BCell::INTERNALIZED:
				internalized_count += grid->b.lists[l].item(i)->count; break;
			case BCell::PRES_II:
				pres_ii_count += grid->b.lists[l].item(i)->count; break;
			case BCell::DUPLICATING:
				duplicating_count += grid->b.lists[l].item(i)->count; break;
			case BCell::PLASMA:
				plasma_count += grid->b.lists[l].item(i)->count; break;
		}
		if( ((BCell *)grid->b.lists[l].item(i)->item)->isMemory() ){
			mem_count += grid->b.lists[l].item(i)->count;
		}
    }
  }
  rstream << naive_count << '\t' << internalized_count << '\t'
	<< pres_ii_count << '\t' << duplicating_count << '\t'
	<< plasma_count << '\t' << mem_count;
  return rstream.str();
}

std::string Grid::ab_details( ){
  std::stringstream rstream;
  grid->ab.init_iterator();
  while( grid->ab.current_entity() != 0 ){
  	rstream << ((Antibody *) grid->ab.current_entity())->paratope
		<< ": "
		<< grid->ab.current_count()
		<< "\t";
	grid->ab.next();
  }
  return rstream.str();
}

std::string Grid::t_stats( ){
  std::stringstream rstream;
  int naive_count = 0, duplicating_count = 0,
		effector_count = 0, memory_count =0;
  // native tcells
  for( int l = 0 ; l < grid->t.size ; l ++ ){
    for( int i = 0 ; i < grid->t.lists[l].length() ; i ++ ){
		switch( ((TCell *)grid->t.lists[l].item(i)->item)->state ){
			case TCell::NAIVE:
				naive_count += grid->t.lists[l].item(i)->count; break;
			case TCell::DUPLICATING:
				duplicating_count += grid->t.lists[l].item(i)->count; break;
			case TCell::EFFECTOR:
				effector_count += grid->t.lists[l].item(i)->count; break;
			case TCell::MEMORY:
				memory_count += grid->t.lists[l].item(i)->count; break;
		}
    }
  }
  rstream << naive_count << '\t' << duplicating_count << '\t' 
			<< effector_count << '\t' << memory_count;
  return rstream.str();
}
