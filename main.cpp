#include "entities/tcell.h" 
#include "entities/bcell.h" 
#include "lists/entitylist.h"
#include "grid/grid.h"
#include "sim/sim.h"
#include "settings/settings.h"
#include "threadlocals/threadlocals.h"
#include "affinity/affinity.h"
#include <iostream>
#include <fstream>

using namespace std;

extern Grid * grid; 

void affinity_test(){
  for( int i = 0 ; i <= Settings::nbitstr ; i ++ )
    cout << Affinity::affinity( 0, (1 << i) - 1 ) << endl; 
  exit( 0 );
}

int main( int argc, char ** argv ){
  char frame[255]; int i; bool is_interactive_mode = 0; int times_ag_injected = 0;

  if( argc == 2 && strcmp(argv[1],"-h") == 0 ){
    cout << "LImmSim v. 0.2" << endl;
    cout << "Written by Johannes Textor, University of Luebeck, Germany" << endl;
    cout << "Contact: textor@informatik.uni-luebeck.de\n" << endl;

    cout << "This program is free software; you can redistribute it and/or" << endl;
    cout << "modify it under the terms of the GNU General Public License" << endl;
    cout << "as published by the Free Software Foundation; either version 2" << endl;
    cout << "of the License, or (at your option) any later version.\n" << endl;


    cout << "Usage:\n\nlimmsim\nlimmsim -i" << endl;
    cout << "For config options syntax, see settings/default-settings.txt." << endl;
    cout << "Using the \"-i\" switch, limmsim will read options from" << endl;
    cout << "standard input until EOF or a line with a single dot ('.')" << endl;
    cout << "is encountered. When no argument is given, the builtin" << endl;
    cout << "default parameters will be used." << endl;
    exit( 0 );
  } else if( argc == 2 && strcmp(argv[1],"-i") == 0 ){
    Settings::configure();
    is_interactive_mode = 1;
  }

  cout << "Initializing Grid ... ";
  grid = new Grid();
  ThreadLocals::initialize();
  Affinity::initialize();
  cout << "done." << endl;

  cout << "Creating Immune System ... ";
  grid->generate_is();
  cout << "done." << endl;

  cout << "Starting simulation." << endl;
  fflush(stdout);

  // randomly generate infecting ag
  Antigen disease;
  for( int i = 0 ; i < Settings::ag_epitopes ; i ++ )
    disease.epitopes[i] = ThreadLocals::rng[0].randInt()%(1<<Settings::nbitstr);

  for( int i = 0 ; i < Settings::ag_peptides ; i ++ )
    disease.peptides[i] = ThreadLocals::rng[0].randInt()%(1<<Settings::nbitstr);
  
  // open output files
  ofstream log,tdat,bdat,agdat,abdat,madat,icdat,abdetails;
  log.open("output/data/log.dat");
  tdat.open("output/data/t.dat");
  bdat.open("output/data/b.dat");
  agdat.open("output/data/ag.dat");
  abdat.open("output/data/ab.dat");
  madat.open("output/data/ma.dat");
  icdat.open("output/data/ic.dat");
  abdetails.open("output/data/abdetails.dat");

  for( i = 0 ; i <= Settings::timesteps ; i ++ ){
    cout << "step " << i << endl; 

    if( Settings::dump_bmps ){
      snprintf( frame, 255, "%06d", i );
      grid->dump_bmps( frame ); 
    }

    tdat << (i) << "\t" << grid->t_stats() << endl;
    bdat << (i) << "\t" << grid->b_stats() << endl;
    agdat << (i) << "\t" << grid->ag.count_elements() << endl;
    abdat << (i) << "\t" << grid->ab.count_elements() << endl;
    icdat << (i) << "\t" << grid->ic.count_elements() << endl;
    madat << (i) << "\t" << grid->ma_stats() << endl;
    abdetails << (i) << "\t" << grid->ab_details() << endl;

    if( is_interactive_mode ){
      cout << grid->t_stats() << "\n" << grid->b_stats() << "\n" 
        << grid->ag.count_elements() << "\n" << grid->ab.count_elements() << "\n" 
        << grid->ic.count_elements() << "\n" << grid->ma_stats() << "\n";
    }

    if( times_ag_injected < Settings::ag_nr_injections && i == Settings::ag_injection_schedule[times_ag_injected] ){
    	log << "Injecting Ag " << disease.epitopes[0] << " at timestep " << i << endl;
    	grid->inject_ag( &disease );
    	times_ag_injected ++;
    }
    if( i < Settings::timesteps )
      sim();
  }

  log.close(); tdat.close();bdat.close();agdat.close();madat.close();icdat.close();
  abdetails.close();
}
