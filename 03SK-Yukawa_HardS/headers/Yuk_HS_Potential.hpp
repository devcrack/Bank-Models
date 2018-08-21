#include "Common_Potential.hpp"

class Yuk_HS_Potential : public Common_Potential {
private:
  double z_yuk;
  double temperature;
public:
  Yuk_HS_Potential();
  Yuk_HS_Potential(double _vol_fraction, double _temperature);
  ~Yuk_HS_Potential();
  void sharmax2_attractive_yukawa();
  void sharmax2_attractive_yukawa(double tmp);
  void set_config(std::vector<std::string>_lines);
  void print_optimal_values();
  void set_config_2(std::vector<std::string>_lines);
};
  
