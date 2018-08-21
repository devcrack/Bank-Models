#ifndef SS_POTENTIAL_HPP
#define SS_POTENTIAL_HPP

#include "HS_Potential.hpp"
#include "HPP_headers.hpp"

class SS_Potential : public HS_Potential {
private:
  double dss;
  double tmp_Init;
public:
  SS_Potential();
  double blip(double temp);
  void set_up_start_dss();
  void set_up_from_config_file(std::vector<std::string>_lines);
  void adjust_phi();
  void Structure_Factor_SSphere();
  int get_nKas();
  double get_Phi();
  double get_dss();
  double get_deltaK();
  double* get_SK();
  
};

#endif
