/** @file */ 
#include<hard_sphere_dynamics.hpp>
#include<iostream>
#include<stdio.h>
#include <cstdlib>
#include<sstream>
#include<string.h>

/**
 * Funcion principal de donde parte la ejecucion para el calculo del factor de estructura para 
 * para esfera dura.
 */
int main(int argc, const char *argv[]) {
  Hard_Sphere_Dynamic *hrd_sphr_dynamic;
  double eta_phi_volume_fraction = 0.0;
  File_Manage file;
  std::ostringstream str_converter;

  if(argc == 2) { //The user gives the volumen fraction
    if(!strcmp(*(argv + 1), "load")){
      std::vector<std::string> lines;
      lines = file.get_just_numbrs(file.read_file_line_by_line_g("./config/options_Hard_Sphere.txt"));
      eta_phi_volume_fraction = atof(lines[0].c_str());
    }
    else 
      eta_phi_volume_fraction = atof(*(argv+1));
    std::cout << "Volumen Fraction = " << eta_phi_volume_fraction <<"\n";
    hrd_sphr_dynamic = new Hard_Sphere_Dynamic(true, true, true,eta_phi_volume_fraction);
  }
  else
    hrd_sphr_dynamic = new Hard_Sphere_Dynamic(true, true, true);
  str_converter << hrd_sphr_dynamic->get_volumen_fraction();
  file.save_config(str_converter.str());
  hrd_sphr_dynamic->start();

  delete hrd_sphr_dynamic;
  return 0; 
}
