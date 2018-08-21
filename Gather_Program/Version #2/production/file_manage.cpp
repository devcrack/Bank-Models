#include "HPP_headers.hpp"

/**
 *brief Constructor default de la clase.
*/
File_Manage :: File_Manage():
  truncate_file(false)
{
}

/**
 *Crea un archivo con el nombre especificado, ademas de indicar si este mismo desea truncarlo o no.
 @param str nombre del archivo, este incluye la ruta del mismo.
 @param truncate simple bandera para indicar si desea que este sea truncado o no.
 */
File_Manage::File_Manage(const std::string& str, bool truncate):
  file_name(str),
  truncate_file(truncate)
{
  
}

/**
 *Destructor default de la clase.
*/
File_Manage :: ~File_Manage() { delete [] _directory;}

/**
 *Guarda el factor de estructura para esfera dura.
 @param sk arreglo que contiene los datos del factor de estructura.
 @param k arreglo que contiene datos relativos al factor de estructura.
 @param kpoints numero de puntos que contiene el factor de estructura.
 @param species numero de especies en el sistema.
 */
void 
File_Manage :: save_sk(double *sk, double* k, int Kpoints, int species) {
  std::cout << this->file_name;
  outfile.open(this->file_name.c_str());
  this->outfile << "#k" << "\t\t\t\t\t";
  this->outfile << "0S(k)_ 1_ 1" << "\n";
  
  for (int i = 0; i < Kpoints; i++) {
    this->outfile << std::fixed << std::setprecision(17) << k[i] << "\t\t\t";
    for(int j = 0; j < species; j++) 
      for (int k= 0; k < species; k++)		
	this->outfile << std::fixed << std::setprecision(17) << *(sk + i * species * species + j * species * k) << "\n"; // Original Line Code
  }
  outfile.clear();
  outfile.close();
}

/*
 *Definitive save SK
 */
void
File_Manage :: save_SK_MONO(double *sk, double* k, int Kpoints) {
  std::cout << this->file_name;
  outfile.open(this->file_name.c_str());
  this->outfile << "#k" << "\t\t\t";
  this->outfile << "0S(k)_ 1_ 1" << "\n";
  for (int i = 0; i < Kpoints; ++i) {
    this->outfile << std::scientific << k[i] << "\t\t";
    this->outfile << std::scientific << sk[i] << "\n";
  }
  outfile.clear();
  outfile.close();  
}


// void
// File_Manage::save_Sk_HSphere(double kMax,double deltaK, double *sk) {
//     double k,i;
//     int index;
  
//   outfile.open(this->file_name.c_str());
//   for(i = 0.0, index = 0; i < kMax; i+=deltaK, ++index) {
//     this->outfile << std::fixed << std::setprecision(std::cout.precision()) << i << "\t\t";
//     this->outfile << std::fixed << std::setprecision(18) << sk[index] << "\n"; // Original Line Code
//     printf("Sk[%d] = %.17f\n", index,sk[index]);
//   }
//   outfile.clear();
//   outfile.close();  
// }





bool
File_Manage :: read_Sk() {
  return false;
}

void
File_Manage:: read_file_line_by_line_p() {  
  std::string line;
  int nmb_lne = 1;
  
  infile.open("options.txt");
  while(std::getline(this->infile, line)) {
    std::cout<<nmb_lne << "\t";
    std::cout<<line << "\n";
    ++nmb_lne;
  }
}


/**
 *Read a text file line by line, load this information in a vector and return it.
 */
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
  std::cout<<"Read the Options File\n";

  
  return lines;
}

/**
 * Get  the current date and time and set this in the _date array.
 * this must to read the options file for read some properties as like final 
 * tempeture, incial tempeture between others.
 * THIS IS CAN CHANGE IS A HORRIBLE SNIPET CODE
*/
void
File_Manage::set_date_FILE_from_Options_file() {
  time_t t = time(NULL);
  std::cout<<"Load the date\n";
  struct tm tm = *localtime(&t);
  std::vector<std::string> lines = this->read_file_line_by_line_g("options.txt");
  sprintf(this->_date,"%f-%f %d-%d %d:%d:%d",atof(lines[5].c_str()), atof(lines[6].c_str()), tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  std::cout <<"Current date:" << this->_date <<"\n";
  
}

/**
 *Just create a simple normal directory.
 @param dir_name - Is the name of the directory.
 */
void
File_Manage::create_directory(const char *dir_name[]) {
  mkdir(_date, 0770);
}


/**
 *Add a file to a directory with the name in _date
 */
void
File_Manage::Add_a_File_to_Directory_Date(const char file_name[]) {
  sprintf(this->_directory_name, "%s/%s",this->_date, file_name);
  this->_directory = fopen(this->_directory_name, "w");  
  fclose(this->_directory);
}

char*
File_Manage::get_date() {  
  return this->_date;
}


void
File_Manage::save_Sk_SSPhere(int nKas,double phi,double deltaK,double _dss ,double *sk) {
  double k,i;
  int index;
  
  outfile.open(this->file_name.c_str());
  
  for(
       index = 0, i = 0.0;
       index < nKas;
       i+=(deltaK * _dss), ++index
      ) {
    k = i * pow( (phi - (phi * phi) / 16.0) / phi, 1/3.0);
    // this->outfile << std::fixed << std::setprecision(18) << sk[index] << "\n"; // Original Line Code
    this->outfile << std::fixed << std::setprecision(std::cout.precision()) << index *deltaK << "\t\t";
    this->outfile << std::fixed << std::setprecision(std::cout.precision()) << k << "\t\t";
    this->outfile << std::fixed << std::setprecision(std::cout.precision()) << k / _dss << "\t\t";
    this->outfile << std::fixed << std::setprecision(18) << sk[index] << "\n"; // Original Line Code
    
  }
  outfile.clear();
  outfile.close();
}


/*
 *Source: https://stackoverflow.com/questions/16749069/c-split-string-by-regex
 */
std::vector<std::string>&
File_Manage :: split(const std::string &s, char delim,std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (item.length() > 0) {
            elems.push_back(item);  
        }
    }
    return elems;
}

/*
 *Source: https://stackoverflow.com/questions/16749069/c-split-string-by-regex
 */
std::vector<std::string>
File_Manage::my_split(const std::string &str, char delim) {
  std::vector<std::string> elems;
  this->split(str, delim, elems);
  
  return elems;
}


/*
 *This is method gets the part of vector that corresponds only to numbers.
*/
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


void
File_Manage::print_file(const std::string &str) {
  std::string line;

  infile.open(str.c_str());
  while(std::getline(this->infile, line))
    std::cout << line << "\n"; 
  infile.clear();
  infile.close();
}
