#include "hard_sphere_mono.hpp"


Hard_Sphere_Mono :: Hard_Sphere_Mono()
  :
  Static_Variables(pow(2,12), 1,1.0e-2),
  Sys_Variables(1)
{
  
}

Hard_Sphere_Mono :: Hard_Sphere_Mono(double _vol_fraction)
  :
  Static_Variables(pow(2,12), 1,1.0e-2),
  Sys_Variables(1, _vol_fraction)
{
  
}


/*
 * Subroutine for the calculation of a static number of equally spaced wave vectors
 * Calc_static_k(dk,kpoints)
 */
void
Hard_Sphere_Mono :: calc_static_K() {
  Sys_Variables::D0M_ini();
  Sys_Variables::rho_ini();
  for (int i = 0; i < Static_Variables :: _kpoints ; i++)     
    Static_Variables :: _k[i] = (i+1) * Static_Variables ::_dk;
}


void 
Hard_Sphere_Mono :: calc_Sk_hs_py_mono( bool VW_op) {
  double *_k_vw;
  double
    eta_vw,
    _calc_1,
    _calc_2,
    _calc_3,
    _calc_4,
    _calc_sin,
    _calc_cos,
    tmp_aux;
  // Param checks
  if(Sys_Variables :: _SDimen != 3 || Sys_Variables :: _Species != 1) {
    printf("Error, cannot compute the system direct Correlaction function");      
    return;
  }

  _k_vw = new double[Static_Variables ::_kpoints];
  if(VW_op) {
    eta_vw = Sys_Variables :: _eta[0] * (1.e0 - Sys_Variables::_eta[0] / 16.e0);
    for (int i = 0 ; i < Static_Variables::_kpoints; i++) 
      _k_vw[i] = Static_Variables :: _k[i] * pow(eta_vw / Sys_Variables :: _eta[0], 1.e0 / 3.e0);
  }
  else {
    eta_vw = Sys_Variables ::_eta[0];
    for (int i = 0; i < Static_Variables :: _kpoints; i++) 
      _k_vw[i] = Static_Variables :: _k[i];	
  }
  _calc_4 = pow(1.0 - eta_vw, 4);
  _calc_1 = -pow(1.0 + 2.0 * eta_vw,2) / _calc_4;
  _calc_2 = ( 6.0 *eta_vw * pow(1.0 + eta_vw / 2.0, 2)  )  / _calc_4;
  _calc_3 = ( -eta_vw * pow(1.0 + 2.0 * eta_vw,2) / 2.0 ) / _calc_4;

  for (int i = 0; i < Static_Variables :: _kpoints; i++) {
    _calc_sin = sin(_k_vw[i]);
    _calc_cos = cos(_k_vw[i]);
    *(Static_Variables::ck + i *_Species) =
      (_calc_1 * (_calc_sin - _k_vw[i] * _calc_cos) / pow(_k_vw[i],2.e0))      
      +
      (_calc_2 *( ( 2.e0 * _k_vw[i]* _calc_sin ) + ( ( -pow(_k_vw[i], 2.e0) + 2.e0) * _calc_cos) - 2.e0) / pow(_k_vw[i],3.e0))
      +
      ( _calc_3
	*
	(
	   (4.e0 * pow(_k_vw[i], 3.e0) - 24.0 * _k_vw[i]) * _calc_sin
	   +
	   ( (-pow(_k_vw[i], 4.e0) + 12.e0 * pow(_k_vw[i], 2.e0) - 24.e0) * _calc_cos)
	   +
	   24.e0
	)
	/ pow(_k_vw[i], 5)
      );
    *(Static_Variables::ck + i *_Species)    = 24.e0 * eta_vw * *(Static_Variables::ck + i *_Species) / _k_vw[i];
    *(Static_Variables :: sk + i * _Species) = 1.e0 / (1.e0 - *(Static_Variables::ck + i * _Species));
    *(Static_Variables :: ski + i * _Species) = 1.e0 - *(Static_Variables::ck + i * _Species);
  }

}


double*
Hard_Sphere_Mono::get_K(){
  return Static_Variables::_k;
}

double*
Hard_Sphere_Mono::get_sk (){
  return Static_Variables::sk;
}

int
Hard_Sphere_Mono::get_KPoints()  {
  return Static_Variables::_kpoints;
}

int
Hard_Sphere_Mono::get_Species() {
  return Sys_Variables::_Species;
}
  
double
Hard_Sphere_Mono::get_phi() {
  return Sys_Variables::_eta[0];
}
