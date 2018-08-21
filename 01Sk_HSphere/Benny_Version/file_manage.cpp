#include "file_manage.hpp"


File_Manage :: File_Manage() {
}

File_Manage::File_Manage(const std::string& str, bool truncate):
  file_name(str),
  truncate_file(truncate)
{
  
}

void 
File_Manage :: save_sk(double *sk, double* k, int Kpoints, int species) {
  outfile.open(this->file_name.c_str());
  this->outfile << "#k" << "\t\t\t\t\t";
  this->outfile << "0S(k)_ 1_ 1" << "\n";
  
  for (int i = 0; i < Kpoints; i++) {
    // this->outfile << std::fixed << std::setprecision(17) << k[i] << "\t\t\t";    
    this->outfile << std::scientific << k[i] << "\t\t\t";
    for(int j = 0; j < species; j++) 
      for (int k= 0; k < species; k++)
	this->outfile << std::scientific << *(sk + i * species * species + j * species * k) << "\n"; // Original Line Code
	// this->outfile << std::fixed << std::setprecision(17) << *(sk + i * species * species + j * species * k) << "\n"; // Original Line Code

    
  }
  outfile.clear();
  outfile.close();
}



bool
File_Manage :: read_Sk() {
  return false;
}

/**
 * Carga la configuracion al archivo de configuracion. Dicha configuracion solo contiene el factor 
 * de estructura. Los demas parametros del modelo permancen constantes.
*/
void
File_Manage::save_config(const std::string vol_fract) {
  this->outfile.open("./config/options_Hard_Sphere.txt");
  this->outfile << vol_fract << "   ;Fraccion de Volumen\n";
  this->outfile.clear();
  this->outfile.close();
}


std::vector<std::string>
File_Manage::get_just_numbrs(std::vector<std::string>lines) {
  std::vector<std::string>split,nmbrs;

  for (int i = 0; i < lines.size(); ++i) {
    split = this->my_split(lines[i],';');
    if(split.size() > 0)
      nmbrs.push_back(split[0]);
  }

  return nmbrs;
}

std::vector<std::string>
File_Manage:: read_file_line_by_line_g(std::string str) {
  std::vector<std::string>lines;
  std::string line;

  infile.open(str.c_str());
  while(std::getline(this->infile, line)) {
    lines.push_back(line);
  }
  infile.clear();
  infile.close();
  
  return lines;
}

std::vector<std::string>
File_Manage::my_split(const std::string &str, char delim) {
  std::vector<std::string> elems;
  this->split(str, delim, elems);
  
  return elems;
}


std::vector<std::string>&
File_Manage :: split(const std::string &s, char delim,std::vector<std::string> &elems) {
  std::stringstream ss(s);	// This line require #include<sstream>
  std::string item;
  while (std::getline(ss, item, delim)) {
    if (item.length() > 0) {
      elems.push_back(item);  
    }
  }
  return elems;
}

