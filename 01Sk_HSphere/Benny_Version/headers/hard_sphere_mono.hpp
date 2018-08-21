#include "static_variables.hpp"
#include "system_variables.hpp"


#include <math.h>
#include<stdio.h>
/**@brief realiza todas las operaciones para el calculo de Sk de esfera dura
   
*/
class Hard_Sphere_Mono : public Static_Variables, public Sys_Variables {
public:
  Hard_Sphere_Mono();
  Hard_Sphere_Mono(double _vol_fraction);
  void calc_static_K();
  void calc_Sk_hs_py_mono(bool VW_op);
  double* get_K();
  int get_KPoints();
  int get_Species();
  double* get_sk ();
  double get_phi();
};
