#include "settings.h"

int Settings::dump_bmps = 0;

int Settings::threads = 2;
int Settings::nbitstr = 16;

int Settings::gridsize = 32;
int Settings::timesteps = 5000;

double Settings::vdiff_ag = 0.1;
double Settings::vdiff_ic = 0.1;
double Settings::vdiff_ab = 0.5;
double Settings::vdiff_t = 0.02;
double Settings::vdiff_b = 0.02;
double Settings::vdiff_ma = 0.01;

double Settings::tau_ag = 10000.0;
double Settings::tau_ab = 690.0;
double Settings::tau_ic = 1000.0;
double Settings::tau_t = 200.0;
double Settings::tau_t_eff = 50.0;
double Settings::tau_b = 200.0;
//double Settings::tau_ma = 50.0;
double Settings::tau_plasma = 100.0;
double Settings::tau_memory = 10000.0;

int Settings::th_diffusion = 32;
int Settings::th_reaction = 32;

// in cimmsim alles 10x mehr 
int Settings::capacity = 64;
int Settings::init_macrophages = 12;
int Settings::init_tcells = 4;
int Settings::init_bcells = 4;

int Settings::ag_epitopes = 1;
int Settings::ag_peptides = 3;

double Settings::new_bcells = Settings::init_bcells / Settings::tau_b;
double Settings::new_tcells = Settings::init_tcells / Settings::tau_t;

int Settings::initial_listsize = 25;
int Settings::nr_reactions = 17;

double Settings::p_ma_ag = 0.0005;
double Settings::p_ma_ic = 1;
double Settings::p_b_ag = 1;
double Settings::p_b_t = 3;
double Settings::p_ma_t = 0.8;
double Settings::p_remove_ag = 0.01;
double Settings::p_ag_ab = 0.5;
double Settings::p_b_become_mem = 1;
double Settings::p_b_rec_mutate = 0.0;
double Settings::p_t_eff_become_mem = 0.0;
double Settings::p_t_eff_become_naive = 0.0;
double Settings::p_ma_become_mem = 0.0;

double Settings::afflevel = 0.05;
int Settings::minmatch = 14;
double Settings::affenhance = 0.77;
int Settings::mhc = 123;

double Settings::v_prolif_b = 0.05;
double Settings::v_prolif_t = 0.025;
double Settings::v_prolif_ag = 0.025;
int Settings::prolif_cycles_t = 6;
int Settings::prolif_cycles_b = 6;
int Settings::ab_secrete = 16;

int Settings::ag_nr_injections = 0;
int * Settings::ag_injection_schedule = 0;

void Settings::configure(){
  char buf[1000];
  char buf2[1000];
  char * tokbuf;
  int i,ipar;
  double dpar; 

  while( fgets( buf, 1000, stdin ) != 0 && strcmp( buf, ".\n" ) != 0 ){
    sscanf( buf, "%s", buf2 );
    if( strcmp( buf2, "dump_bmps" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("dump_bmps: %d\n", ipar);
      dump_bmps = ipar != 0 ? 1 : 0;
    } else if( strcmp( buf2, "p_ma_ag" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_ma_ag: %lf\n", dpar);
      p_ma_ag = dpar;
    } else if( strcmp( buf2, "p_ma_ic" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_ma_ic: %lf\n", dpar);
      p_ma_ic = dpar;
    } else if( strcmp( buf2, "p_b_ag" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_b_ag: %lf\n", dpar);
      p_b_ag = dpar;
    } else if( strcmp( buf2, "p_b_t" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_b_t: %lf\n", dpar);
      p_b_t = dpar;
    } else if( strcmp( buf2, "p_ma_t" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_ma_t: %lf\n", dpar);
      p_ma_t = dpar;
    } else if( strcmp( buf2, "p_remove_ag" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_remove_ag: %lf\n", dpar);
      p_remove_ag = dpar;
    } else if( strcmp( buf2, "p_ag_ab" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_ag_ab: %lf\n", dpar);
      p_ag_ab = dpar;
    } else if( strcmp( buf2, "p_b_become_mem" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_b_become_mem: %lf\n", dpar);
      p_b_become_mem = dpar;
    } else if( strcmp( buf2, "p_ma_become_mem" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_ma_become_mem: %lf\n", dpar);
      p_ma_become_mem = dpar;
    } else if( strcmp( buf2, "p_t_eff_become_mem" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_t_eff_become_mem: %lf\n", dpar);
      p_t_eff_become_mem = dpar;
    } else if( strcmp( buf2, "p_t_eff_become_naive" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_t_eff_become_naive: %lf\n", dpar);
      p_t_eff_become_naive = dpar;
    } else if( strcmp( buf2, "p_b_rec_mutate" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("p_b_rec_mutate: %lf\n", dpar);
      p_b_rec_mutate = dpar;
    } else if( strcmp( buf2, "afflevel" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("afflevel: %lf\n", dpar);
      afflevel = dpar;
    } else if( strcmp( buf2, "affenhance" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("affenhance: %lf\n", dpar);
      affenhance = dpar;
    } else if( strcmp( buf2, "v_prolif_b" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("v_prolif_b: %lf\n", dpar);
      v_prolif_b = dpar;
    } else if( strcmp( buf2, "v_prolif_t" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("v_prolif_t: %lf\n", dpar);
      v_prolif_t = dpar;
    } else if( strcmp( buf2, "v_prolif_ag" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("v_prolif_ag: %lf\n", dpar);
      v_prolif_ag = dpar;
    } else if( strcmp( buf2, "vdiff_ag" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("vdiff_ag: %lf\n", dpar);
      vdiff_ag = dpar;
    } else if( strcmp( buf2, "vdiff_ic" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("vdiff_ic: %lf\n", dpar);
      vdiff_ic = dpar;
    } else if( strcmp( buf2, "vdiff_ab" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("vdiff_ab: %lf\n", dpar);
      vdiff_ab = dpar;
    } else if( strcmp( buf2, "vdiff_t" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("vdiff_t: %lf\n", dpar);
      vdiff_t = dpar;
    } else if( strcmp( buf2, "vdiff_b" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("vdiff_b: %lf\n", dpar);
      vdiff_b = dpar;
    } else if( strcmp( buf2, "vdiff_ma" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("vdiff_ma: %lf\n", dpar);
      vdiff_ma = dpar;
    } else if( strcmp( buf2, "tau_ag" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("tau_ag: %lf\n", dpar);
      tau_ag = dpar;
    } else if( strcmp( buf2, "tau_ab" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("tau_ab: %lf\n", dpar);
      tau_ab = dpar;
    } else if( strcmp( buf2, "tau_ic" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("tau_ic: %lf\n", dpar);
      tau_ic = dpar;
    } else if( strcmp( buf2, "tau_t" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("tau_t: %lf\n", dpar);
      tau_t = dpar;
      new_tcells = init_tcells / tau_t;
    } else if( strcmp( buf2, "tau_t_eff" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("tau_t_armed: %lf\n", dpar);
      tau_t_eff = dpar;
    } else if( strcmp( buf2, "tau_b" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("tau_b: %lf\n", dpar);
      tau_b = dpar;
      new_bcells = init_bcells / tau_b;
    } else if( strcmp( buf2, "tau_memory" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("tau_memory: %lf\n", dpar);
      tau_memory = dpar;
    } else if( strcmp( buf2, "tau_plasma" ) == 0 ){
      sscanf( buf, "%s %lf", buf2, &dpar );
      printf("tau_plasma: %lf\n", dpar);
      tau_plasma = dpar;
    } else if( strcmp( buf2, "threads" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("threads: %d\n", ipar);
      threads = ipar;
    } else if( strcmp( buf2, "nbitstr" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("nbitstr: %d\n", ipar);
      nbitstr = ipar;
    } else if( strcmp( buf2, "gridsize" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("gridsize: %d\n", ipar);
      gridsize = ipar;
    } else if( strcmp( buf2, "timesteps" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("timesteps: %d\n", ipar);
      timesteps = ipar;
    } else if( strcmp( buf2, "th_diffusion" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("th_diffusion: %d\n", ipar);
      th_diffusion = ipar;
    } else if( strcmp( buf2, "th_reaction" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("th_reaction: %d\n", ipar);
      th_reaction = ipar;
    } else if( strcmp( buf2, "capacity" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("capacity: %d\n", ipar);
      capacity = ipar;
    } else if( strcmp( buf2, "init_macrophages" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("init_macrophages: %d\n", ipar);
      init_macrophages = ipar;
    } else if( strcmp( buf2, "init_tcells" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("init_tcells: %d\n", ipar);
      init_tcells = ipar;
      new_tcells = init_tcells / tau_t;
    } else if( strcmp( buf2, "init_bcells" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("init_bcells: %d\n", ipar);
      init_bcells = ipar;
      new_bcells = init_bcells / tau_b;
    } else if( strcmp( buf2, "ag_peptides" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("ag_peptides: %d\n", ipar);
      ag_peptides = ipar;
    } else if( strcmp( buf2, "minmatch" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("minmatch: %d\n", ipar);
      minmatch = ipar;
    } else if( strcmp( buf2, "mhc" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("mhc: %d\n", ipar);
      mhc = ipar;
    } else if( strcmp( buf2, "prolif_cycles_t" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("prolif_cycles_t: %d\n", ipar);
      prolif_cycles_t = ipar;
    } else if( strcmp( buf2, "prolif_cycles_b" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("prolif_cycles_b: %d\n", ipar);
      prolif_cycles_b = ipar;
    } else if( strcmp( buf2, "ab_secrete" ) == 0 ){
      sscanf( buf, "%s %d", buf2, &ipar );
      printf("ab_secrete: %d\n", ipar);
      ab_secrete = ipar;
    }
    else if( strcmp( buf2, "ag_injection_schedule" ) == 0 ){
      strcpy( buf2, buf );
      ag_nr_injections = 0;
      strtok( buf, " " );
      while( (tokbuf = strtok( NULL, " " )) ) ag_nr_injections ++; 
      if( ag_injection_schedule ) delete [] ag_injection_schedule;
      ag_injection_schedule = new int[ag_nr_injections];
      i = 0;
      strtok( buf2, " " );
      while( (tokbuf = strtok( NULL, " " )) ){
		sscanf( tokbuf, "%d", &ipar );
		ag_injection_schedule[i++] = ipar; 
      }
      printf("ag_injection_schedule: ");
      for( i = 0 ; i < ag_nr_injections ; i ++ ){
      	printf( "%d ",ag_injection_schedule[i] );
      }
      printf("\n");
    }
    else{
      printf("unrecognized option : %s\n", buf);
    }
  }
}


