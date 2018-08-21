#include "static_variables.hpp"
#include<stdio.h>
/**
 *This is the subroutine: Static_Variables_mem_alloc(kpoints,species)
 Static_Variables_mem_alloc(kpoints,species) ![X]
  Subroutine Static_Variables_mem_alloc(Ialloc_k,Ialloc_s)

 */

Static_Variables :: Static_Variables(int _Ialloc_k, int _Ialloc_s, double dk)
  :
  _kpoints(_Ialloc_k),
  _dk(dk),
  sk  (new double[_Ialloc_k * _Ialloc_s * _Ialloc_s]),
  ski (new double[_Ialloc_k * _Ialloc_s * _Ialloc_s]),
  ck  (new double[_Ialloc_k * _Ialloc_s * _Ialloc_s]),
  hk  (new double[_Ialloc_k * _Ialloc_s * _Ialloc_s]),
  _k   (new double[_Ialloc_k])
{
}


Static_Variables:: ~Static_Variables() {
  delete[] sk;
  delete[] ski;
  delete[] ck;
  delete[] hk;
}


