#include "affinity.h"

void Affinity::initialize(){
  int i;

  mhc_affinity_lut = new double[(Settings::nbitstr>>1)+1];
  affinity_lut = new double[Settings::nbitstr+1];

  for( i = 0 ; i < Settings::minmatch; i ++ )
    affinity_lut[i] = 0;
  affinity_lut[i++] = Settings::afflevel;
  for( ; i <= Settings::nbitstr ; i ++ ){
    double next = (double)i/(double)(Settings::nbitstr-i+1)*Settings::affenhance*affinity_lut[i-1]; 
    affinity_lut[i] = next < 1.0 ? next : 1.0;
  } 

  mhc_affinity_lut[Settings::nbitstr >> 1] = 1;
  for( i = (Settings::nbitstr >> 1) - 1 ; i >= 0 ; i -- ){
    mhc_affinity_lut[i] = mhc_affinity_lut[i+1] * 0.5;
  }
}

double * Affinity::affinity_lut = 0;
double * Affinity::mhc_affinity_lut = 0;
