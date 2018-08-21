#include "HPP_headers.hpp"


int main(int argc, const char *argv[]) {
  File_Manage *file = new File_Manage("../data/sk_HS_Mono_Patty.dat", false); // Delete this in explicit way
  HS_Potential *ss_potential; 
  std::vector<std::string> lines;

  ss_potential = new HS_Potential();
  lines = file->read_file_line_by_line_g("options_HS.txt");
  ss_potential->set_up_from_config_file(lines);
  ss_potential->Structure_Factor();
  file->save_Sk_HSphere(ss_potential->get_Kmax(), ss_potential->get_deltaK(), ss_potential->get_SK());
  
  return 0;  
}
