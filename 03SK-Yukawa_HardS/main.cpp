
#include  "HPP_headers.hpp"
#include "Yuk_HS_Potential.hpp"

int main(int argc, const char *argv[]) {
  Yuk_HS_Potential *yuk_hs;// =  new Yuk_HS_Potential();
  // File_Manage  *file = new File_Manage("../data/sk_MonoYuk.dat", false); // Delete this in explicit way
    File_Manage  *file = new File_Manage("./data/sk_MonoYuk.dat", false); // Delete this in explicit way
    std::vector<std::string> lines,split;
    if(argc==2)
      if(!(strcmp(*(argv+1), "load"))) {		
	yuk_hs =  new Yuk_HS_Potential();
	lines = file->get_just_numbrs(file->read_file_line_by_line_g("./config/Yuk_options.txt"));
	yuk_hs->set_config_2(lines);
      }
      else
	std::cout<<"Arguments Error\n";
  if(argc == 3) {
    yuk_hs =  new Yuk_HS_Potential(atof(*(argv+1)), atof(*(argv+2)));
    std::cout << "Vol_fraction = " << atof(*(argv+1)) << "Temperature = " << atof(*(argv+2)) << "\n";    
    file->load_config(*(argv+1), *(argv+2));	  
  }
  if(argc == 1) {
    yuk_hs =  new Yuk_HS_Potential();
    std::cout<<"Volumen Fraction = 0.5 \n" << "Temperature = 1\n";
    lines.push_back("0.5");
    lines.push_back("1.e0");
    yuk_hs->set_config_2(lines);
    file->load_config("0.5", "1.e0");	  
  }
  yuk_hs->sharmax2_attractive_yukawa();
  file->save_SK_MONO(yuk_hs->get_SK(), yuk_hs->get_K(), yuk_hs->get_KPoints());

  delete yuk_hs;
  delete file;

  return 0;
}
