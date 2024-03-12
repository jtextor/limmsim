#include "../settings/settings.h"
#include "antigen.h"

#include <iostream>

/** For now, we define equality on antigens in a biological fashion: 
  two antigens are equal iff their epitopes are equal */


Antigen::Antigen( int * _epitopes, int * _peptides ){
  epitopes = new int[Settings::ag_epitopes];
  peptides = new int[Settings::ag_peptides];

  for( int i = 0 ; i < Settings::ag_epitopes ; i ++ )
    epitopes[i] = _epitopes[i];

  for( int i = 0 ; i < Settings::ag_peptides ; i ++ )
    peptides[i] = _peptides[i];
}

Antigen::Antigen(){
  epitopes = new int[Settings::ag_epitopes];
  peptides = new int[Settings::ag_peptides];
}

Antigen::~Antigen(){
  delete [] epitopes;
  delete [] peptides;
}

bool Antigen::operator==(const Entity &other) {
  int * _op = ((Antigen *) &other) -> epitopes;
  for( int i = 0 ; i < Settings::ag_epitopes ; i ++ ){
    if( epitopes[i] != _op[i] ) return false;
  }

  /** for now, epitopes are equal iff peptides are equal to speed up. 
    needs to change if ag ever should mutate. */ 
  return true;
  /*_op = ((Antigen *) &other) -> peptides;
  for( int i = 0 ; i < Settings::ag_peptides ; i ++ ){
    if( peptides[i] != _op[i] ) return false;
  }*/
}

Entity * Antigen::clone(){
  return (Entity *) new Antigen( epitopes, peptides );
}


