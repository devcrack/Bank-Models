#include "HPP_headers.hpp"


class HS_Potential {
protected:
  double         _kMax;
  double         _deltaK;
  double         _phi;
  int            _nKas;
  double         _SK_max;
  double         _SK_min;
  double         *_sk;
public:
  HS_Potential();
  HS_Potential(double _vol_fraction);
  void Structure_Factor();
  double SK();
  double KZero(double phi);
  void strcuture_factor();
  void set_values_from_Config_File();
  void set_up_from_config_file(std::vector<std::string>_lines);
  void set_up_from_config_file_2(std::vector<std::string>_lines);
  double sK(double phi, double k);
  double PZero(double phi, double k);
  double P_One(double phi, double k);
  double P_Two(double phi, double k);
};

