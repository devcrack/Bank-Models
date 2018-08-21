#include  "HPP_headers.hpp"
#include  "SS_Potential.hpp"

int main(int argc, const char *argv[]) {
  SS_Potential *_ss;//= new SS_Po2tential();
  File_Manage  *file = new File_Manage("./data/sk_SSphere.dat", false); // Delete this in explicit way  
  //In this case the user just give the volumen fraction
  if(argc ==  2) {
    if(!strcmp(*(argv+1),"load")){ // Here load the configuration from file;
      std::vector<std::string> lines;
      
      _ss= new SS_Potential();      
      lines = file->get_just_numbrs(file->read_file_line_by_line_g("./config/options_SSphere.txt"));
      std::cout << "Configuration loaded\n";
      _ss->set_up_from_config_file(lines);
    }
    else { // Set a given volume fraction and save it;
      std :: cout << "Volume Fraction = "<< atof(*(argv + 1)) << "\n";
      _ss= new SS_Potential(atof(*(argv + 1)));
      file->load_config(*(argv + 1), "0.1");
    }
  }
  if(argc ==  3) {		// In this case the user gives the volumen fraction and the incial temperature
    std :: cout << "Volume Fraction = "     << atof(*(argv + 1)) << "\n";
    std :: cout << "Initial Temperature = " << atof(*(argv + 2)) << "\n";
    _ss= new SS_Potential(atof(*(argv + 1)), atof(*(argv + 2)));
    // std::cout << "Saving set configuration\n";
    file->load_config(*(argv + 1), *(argv + 2));
  }
  if(argc == 1) { //In this case the user doesn't gives anything and the program operates with the best configuration    
    _ss= new SS_Potential();
    // std::cout << "Saving deafault configuration\n";
    file->load_config("0.56", "0.1"); // default values
  }
  

  _ss->set_up_start_dss();
  _ss->adjust_phi();
  _ss->Structure_Factor_SSphere();
  file->save_Sk_SSPhere_V2(_ss->get_nKas(),_ss->get_Phi(),_ss->get_deltaK(), _ss->get_dss(), _ss->get_SK());

  return 0;
}
  
