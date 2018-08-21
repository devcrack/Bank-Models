#include <sstream>
#include<fstream> 		// Class for use files
#include<string>
#include<iostream>
#include<iomanip>
#include<limits>
#include <cstdlib>
#include <vector>
/** @brief Clase para el manejo de archivos 
    
    Esta clase permite la manipulacion de archivos de tal manera que se pueden realizar las acciones basicas como son,
    abrir guardar y modificar un archivo. Estos archivos basicamente contienen los datos calulados del factor de estructura.

    @author: delcranck
    @date: 28/08/2017 
*/
class File_Manage {
public:
  /**
   *Constructor default
  */
  File_Manage();
    /**
   *Constructor que puede especificar el archivo y si desea truncarse
  */
  File_Manage(const std::string& str, bool truncate);
    /**
   *Guarda el sk
  */
  void save_sk(double *sk, double* k, int Kpoints, int species);
  /**
   * Lee el sk 
  */
  bool read_Sk();  
  
  void save_config(const std::string vol_fract);

  std::vector<std::string> get_just_numbrs(std::vector<std::string>lines);
  std::vector<std::string> read_file_line_by_line_g(std::string str);
  std::vector<std::string> my_split(const std::string &str, char delim);
  std::vector<std::string> &split(const std::string &s, char delim,std::vector<std::string> &elems);
protected:
  std::string file_name;/**<Nombre del archivo, dicho nombre contiene la ruta del mismo*/
  bool truncate_file;/**<Bandera que determina si el archivo debe de ser truncado o no*/
  std :: ofstream outfile;/**<Permite realizar operaciones de escritura*/
  std :: ifstream infile;/**<Permite realizar operaciones de lectura*/
};
