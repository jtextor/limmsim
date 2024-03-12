/* Settings are global variables */ 

#ifndef __SETTINGS

#define __SETTINGS

#include <stdio.h>
#include <string.h>

class Settings{
  public: 

    /** configure limmsim with parameters given on standard input. */ 
    static void configure();

    /** NON-MODEL PARAMS */ 

    static int dump_bmps;
    static int threads;

    /** MODEL PARAMS */ 

    /** number of bits to represent receptors, must be even */ 
    static int nbitstr;

    /** initial population sizes (all per grid point) */
    static int capacity;
    static int init_macrophages;
    static int init_tcells;
    static int init_bcells;
    
    /** antigen injection schedule */
    static int * ag_injection_schedule;
    /** this is deduced from the injection schedule and kept for portability reasons (array size) */ 
    static int ag_nr_injections;

    /** simulation parameters */
    // diameter of grid in x and y direction. must be even! 
    // on smp, best performance for gridsize == 0 mod 3 
    static int gridsize; 
    static int timesteps; 

    /** diffusion speeds */ 
    static double vdiff_ag; // = 0.5;
    static double vdiff_ab; // = 0.5;
    static double vdiff_ic; // = 0.5;
    static double vdiff_t; // = 0.01;
    static double vdiff_b; // = 0.01;
    static double vdiff_ma; // = 0.01;

    /** expected lifetimes */ 
    static double tau_ag;
    static double tau_ab;
    static double tau_ic;
    static double tau_t;
    static double tau_t_eff;
    static double tau_b;
    static double tau_ma; 
    static double tau_plasma;
    // memory b and t cells
    static double tau_memory;
    // probability for b cell to become memory cell after having gone through all prolif cycles
    static double p_b_become_mem;
    // probability for macrophage to become memory cell and keeping the presented antigen
    static double p_ma_become_mem;
    // probability for t effector cell to become memory cell
    static double p_t_eff_become_mem;
    // probability for t effector cell to revert back to a naive cell
    static double p_t_eff_become_naive;

    /** interaction (base) probabilities */ 
    // unspecific phagocytosis of ag by ma
    static double p_ma_ag;
    // unspecific phagocytosis of ic by ma
    static double p_ma_ic;
    // base prop. for specific phagocytosis of ag by b
    static double p_b_ag;
    // base prop. for ag presentation by ma to t 
    static double p_ma_t;
    // base prop. for ag presentation by b to t 
    static double p_b_t;
    // base prop. for ag-ab binding
    static double p_ag_ab;

    // removal of ag peptide from presenting cell
    static double p_remove_ag;
    // per-bit mutation probability of b-cell receptor during mitosis
    static double p_b_rec_mutate;

    /** proliferation speeds */ 
    static double v_prolif_b;
    static double v_prolif_t;
    static double v_prolif_ag;
    static int prolif_cycles_b;
    static int prolif_cycles_t;
    static int ab_secrete;

    /** affinity parameters */ 
    static double afflevel;
    static int minmatch;
    static double affenhance;
    static int mhc;

    /** approximation thresholds */ 
    // threshold for exact diffusion. must be dividable by 6 
    static int th_diffusion;
    // threshold for exact reaction.
    static int th_reaction;

    /** antigen complexity */ 
    static int ag_epitopes;
    static int ag_peptides;

    /** values below are NOT parameters, but convenience functions of other parameters */ 
    // nr of new b and t cells per lattice point per time step
    static double new_bcells;
    static double new_tcells;

    // internal size of lists
    static int initial_listsize;
    // number of reactions
    static int nr_reactions;
};

#endif
