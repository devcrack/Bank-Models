#include "hard_sphere_dynamics.hpp"

/**
 *Constructor default
*/
Hard_Sphere_Dynamic :: Hard_Sphere_Dynamic(bool VW_op, bool Gamm_op, bool SK_writting_op)
  :
  _mono_hard_sphere(),
  _file("./data/sk_HSpheere.dat", true),
  _VW_op(VW_op),
  _Gamm_op(Gamm_op),  
  _SK_writting_op(SK_writting_op)

{

}

/**
 * En este constructo se considera agregar la fraccion de volumen como parametro, dicha fraccion de volumen viene
 * dada por el usuario.
 * @param vol_fraction Fraccion de volumen, determinada por el usuario.
*/
Hard_Sphere_Dynamic :: Hard_Sphere_Dynamic(bool VW_op, bool Gamm_op, bool SK_writting_op, double vol_fraction)
  :
  _mono_hard_sphere(vol_fraction),
  _file("./data/sk_HSpheere.dat", true),
  _VW_op(VW_op),
  _Gamm_op(Gamm_op),  
  _SK_writting_op(SK_writting_op)

{
}


void
Hard_Sphere_Dynamic ::start() {
  this->_mono_hard_sphere.calc_static_K();
  this->_mono_hard_sphere.calc_Sk_hs_py_mono(true);
  this->_file.save_sk(this->_mono_hard_sphere.get_sk(),
  		      this->_mono_hard_sphere.get_K(),
  		      this->_mono_hard_sphere.get_KPoints(),
  		      this->_mono_hard_sphere.get_Species());
  
}




Hard_Sphere_Dynamic :: ~Hard_Sphere_Dynamic() {
  //All here   is destroy automatically
}

double
Hard_Sphere_Dynamic :: get_volumen_fraction() {
  return this->_mono_hard_sphere.get_phi();
}
