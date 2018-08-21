#include "SS_Potential.hpp"


SS_Potential :: SS_Potential() :
  HS_Potential(0.56),
  dss(0.0),
  tmp_Init(0.8)
{
  
}


SS_Potential::SS_Potential(double _vol_fraction):
  HS_Potential(_vol_fraction),
  dss(0.0),
  tmp_Init(0.0)
{
  
}


SS_Potential::SS_Potential(double _vol_fraction, double _tmp_init):
  HS_Potential(_vol_fraction),
  dss(0.0),
  tmp_Init(_tmp_init)
{
  
}

double
SS_Potential :: blip(double temp){
  double
    eps,
    dr,
    r,
    dss;

  eps = 1/temp;
  dr = 1.0 / 50000;
  dss = 0;

  if(temp < 0.0000001)
    dss = 1;
  else {
    for(r = dr; r <= 1; r+=dr)
      dss += dr * r*r * exp(-eps*(1.0/pow(r,12) - 2.0/pow(r,6) + 1));
    dss = pow(1 - 3*dss,1/3.0);
  }
  return dss;
}
 
/*@brief Funcion que carga la configuracion incial desde un archivo dado.
**/
// void
// SS_Potential::set_up_from_config_file(std::vector<std::string>_lines) {
//   HS_Potential::set_up_from_config_file(_lines);
  
//   this->tmp_Init = atof(_lines[5].c_str());
//   printf("Lectura de TmpIni = %.17f\n", this->tmp_Init);
// }


/*@brief Funcion que carga la configuracion incial desde un archivo dado.
**/
void
SS_Potential::set_up_from_config_file(std::vector<std::string>_lines) {
  HS_Potential::set_up_from_config_file_2(_lines);
  
  this->tmp_Init = atof(_lines[1].c_str());
}


void
SS_Potential::set_up_start_dss() {
  this->dss = this->blip(this->tmp_Init);
}

void
SS_Potential :: adjust_phi() {
  double tmp_phi = HS_Potential :: _phi;
  
  HS_Potential :: _phi = tmp_phi *pow(this->dss,3);
}



void
SS_Potential :: Structure_Factor_SSphere() {
  int index = 1;
 
  HS_Potential::_sk[0] = HS_Potential::KZero(HS_Potential::_phi);

  for (
        double i  = HS_Potential :: _deltaK * this->dss ;
        index     < HS_Potential :: _nKas;
        i+= (HS_Potential::_deltaK * this->dss) ,++index
      ) {
    HS_Potential::_sk[index] = HS_Potential::sK(HS_Potential::_phi, i);
  }

  
}


int
SS_Potential :: get_nKas() {
  return HS_Potential :: _nKas;
}

double
SS_Potential :: get_Phi() {
  return HS_Potential:: _phi;
}


double
SS_Potential::get_deltaK() {
  return HS_Potential::_deltaK;
}

double
SS_Potential::get_dss() {
  return this->dss;
}

double*
SS_Potential::get_SK() {
  return HS_Potential::_sk;
}


