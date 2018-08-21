#include "Yuk_HS_Potential.hpp"


Yuk_HS_Potential :: Yuk_HS_Potential() :
  z_yuk(2.e0), 			// This corresponds to Subroutine z_yuk_mono(dummy,species_number)
  Common_Potential()
{
}



Yuk_HS_Potential :: Yuk_HS_Potential(double _vol_fraction, double _temperature) :
  z_yuk(2.e0),
  temperature(_temperature),
  Common_Potential(_vol_fraction)
{
}

Yuk_HS_Potential::~Yuk_HS_Potential() {}

void
Yuk_HS_Potential::sharmax2_attractive_yukawa(double tmp) {
  double
    dmb_ck,
    uyuk;

  Common_Potential::calc_static_K();
  Common_Potential::calc_Sk_hs_py_mono(true);
  for (int i = 0; i < Common_Potential::nKas_Kpoints; ++i) {
    dmb_ck = 1.e0 / Common_Potential::_sk[i];
    // printf("ck = %.16f\n", dmb_ck);    
    uyuk = 24.e0 * _phi_eta * (_k[i] * cos(_k[i]) + z_yuk * sin(_k[i])) / (tmp *_k[i] * (pow(_k[i], 2.e0) + pow(z_yuk,2.e0)) );
    // printf("uyuk = %.16f\n", uyuk);
    dmb_ck = dmb_ck - uyuk ;
    Common_Potential :: _sk[i] = 1.e0 / dmb_ck;
  }
}


void
Yuk_HS_Potential::sharmax2_attractive_yukawa() {
  double
    dmb_ck,
    uyuk;

  Common_Potential::calc_static_K();
  Common_Potential::calc_Sk_hs_py_mono(true);
  for (int i = 0; i < Common_Potential::nKas_Kpoints; ++i) {
    dmb_ck = 1.e0 / Common_Potential::_sk[i];
    uyuk = 24.e0 * _phi_eta * (_k[i] * cos(_k[i]) + z_yuk * sin(_k[i])) / (this->temperature *_k[i] * (pow(_k[i], 2.e0) + pow(z_yuk,2.e0)) );
    dmb_ck = dmb_ck - uyuk ;
    Common_Potential :: _sk[i] = 1.e0 / dmb_ck;
  }
}






void
Yuk_HS_Potential::set_config(std::vector<std::string>_lines) {
  Common_Potential :: set_up_from_config_file(_lines);
  Common_Potential :: nKas_Kpoints = atoi(_lines[3].c_str());
  Common_Potential :: _SDimen      =  atoi(_lines[4].c_str());
  Common_Potential :: _Species     =  atoi(_lines[5].c_str());
  Common_Potential :: _sk          = new double[Common_Potential::nKas_Kpoints];
  Common_Potential :: ck           = new double[Common_Potential::nKas_Kpoints]; 
  Common_Potential :: _k           = new double[Common_Potential::nKas_Kpoints]; 
}

/**
 *Carga la configuracion leida de archivo.
*/
void
Yuk_HS_Potential::set_config_2(std::vector<std::string>_lines) {
  Common_Potential::_phi_eta = atof(_lines[0].c_str());
  this->temperature = atof(_lines[1].c_str());
}

void
Yuk_HS_Potential::print_optimal_values() {
  std::cout << "Volumen Fraction = " << Common_Potential::_phi_eta <<"\n";
  std::cout << "Temperature = " << this->temperature <<"\n";
}
