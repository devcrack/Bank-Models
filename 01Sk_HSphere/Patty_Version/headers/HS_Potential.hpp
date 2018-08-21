#ifndef SCGLE_HPP
#define SCGLE_HPP

#include "HPP_headers.hpp"


class HS_Potential {
private:
  double         _kMax;
  double         _deltaK;
  double         _phi;
  int            _nKas;
  double         _SK_max;
  double         _SK_min;
  double         *_sk;
public:
  HS_Potential();
  void Structure_Factor();
  double SK();
  double KZero(double phi);
  void strcuture_factor();
  void set_values_from_Config_File();
  void set_up_from_config_file(std::vector<std::string>_lines);
  double sK(double phi, double k);
  double PZero(double phi, double k);
  double P_One(double phi, double k);
  double P_Two(double phi, double k);
  int get_nKas();
  double get_deltaK();
  double* get_SK();
  double get_Kmax();
};

#endif
