#include "HPP_headers.hpp"


class Common_Potential {
protected:
  double         _kMax;
  double         deltaK_dk;
  double         _phi_eta;
  int            nKas_Kpoints;
  double         _SK_max;
  double         _SK_min;
  double         *_sk;
  double         *_k;
  int            _SDimen;
  int            _Species;
  double         *ck;
public:
  Common_Potential();
  Common_Potential(double vol_fract);
  ~Common_Potential();
  Common_Potential(int SDimen, int Species);
  void calc_Sk_hs_py_mono( bool VW_op); 
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
  void calc_static_K();
  double* get_SK();
  double* get_K();
  int get_KPoints();
};

