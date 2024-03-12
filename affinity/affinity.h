#ifndef __AFFINITY

#define __AFFINITY

#include "../settings/settings.h"

class Affinity {
  public:

    static inline double mhc_affinity( int u, int v ){
      return mhc_affinity_lut[hamming_distance(u,v,Settings::nbitstr>>1)];
    };
    static inline double affinity( int u, int v ){
      return affinity_lut[hamming_distance(u,v,Settings::nbitstr)];
    };
    static inline double cd4_affinity( int u, int v ){
      return affinity_lut[(Settings::nbitstr>>1)+hamming_distance(u,v,Settings::nbitstr>>1)];
    };
    static void initialize();

    static inline int hamming_distance( int u, int v, int bits ){
      int ret = 0;
      for( int i = 0 ; i < bits ; i ++ ){
        ret += u & 1 ^ v & 1 ;
        u = u >> 1 ; v = v >> 1 ;
      }
      return ret;
    } 

  private:
    static double * affinity_lut;
    static double * mhc_affinity_lut;


};

#endif
