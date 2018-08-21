#include<iostream>

class Static_Variables {
public:

  
  Static_Variables(int _Ialloc_k, int _Ialloc_s, double dk);
  ~Static_Variables();
protected:
  int _kpoints;
  double _dk;
  double *_k;
  //Array 3D with a single pointer
  double  *sk;/**<Guarda el conjunto de datos que conforman al sk de esfera dura.*/
  double *ski;/**<Guarda el conjunto de datos que conforman a otro  sk analogo de esfera dura.*/
  double *ck;
  double *hk;
    
    
  
};
