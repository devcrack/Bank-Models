#include  "HPP_headers.hpp"
#include  "SS_Potential.hpp"

int main(int argc, const char *argv[]) {
  SS_Potential *_ss= new SS_Potential();
  File_Manage  *file = new File_Manage("data/sk_SSphere.dat", false); // Delete this in explicit way
  std::vector<std::string> lines;
    
  lines = file->read_file_line_by_line_g("./config/options_SSphere.txt");
  std::cout << lines.size();
  _ss->set_up_from_config_file(lines);
  _ss->set_up_start_dss();
  _ss->adjust_phi();
  _ss->Structure_Factor_SSphere();
  file->save_Sk_SSPhere(_ss->get_nKas(),_ss->get_Phi(),_ss->get_deltaK(), _ss->get_dss(), _ss->get_SK());
  return 0;
}
  
