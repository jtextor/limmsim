#ifndef __MACROS

#define __MACROS

#include <stdlib.h>

inline int rand_int(int n){
  return rand_r()%n;
}

inline bool coin(double p){
  return (double)rand()/(double)RAND_MAX < p;
}

#endif
