#include<fstream> 		// Class for use files
#include<string>
#include<iostream>
#include<iomanip>
#include<limits>
#include <vector>
#include<cstdlib>
#include <sys/stat.h>

/**
 *@brief Clase que se encarga de toda las funciones que tienen que ver con la manipulacion de archivos, 
 *ademas de la interaccion del sistema de archivos del host.
*/
class File_Manage {
  
public:
  //Constructor Default
  File_Manage();
  //Constructor que indica la ruta y el nombre del archivo ademas de si este, hay que truncarlo.
  File_Manage(const std::string& str, bool truncate);
  //Destructor Default
  ~File_Manage();
  //Guarda el factor de estructura para esfera dura
  void save_sk(double *sk, double* k, int Kpoints, int species);
  //Lee los factores de estructura.
  bool read_Sk();
  //Lee un archivo linea por linea
  void read_file_line_by_line_p();
  //Le un archivo linea por linea y regresa un vector con dichas lineas.
  std::vector<std::string> read_file_line_by_line_g(std::string);
  //Establece la fecha desde el archivo de configuracion.
  void set_date_FILE_from_Options_file();
  //Crea un directorio
  void create_directory(const char *dir_name[]);
  //Agrega un archivo al directorio
  void Add_a_File_to_Directory_Date(const char dir_name[]);
  //Guarda el factor de estructura de esfera suave.
  void save_Sk_SSPhere(int nKas,double phi,double deltaK,double _dss ,double *sk);
  //Obtiene la fecha actual. 
  char* get_date();
  //Sirve para obtener la separacion de un texto dada por un caracter especificado.
  std::vector<std::string> my_split(const std::string &str, char delim);
  //Metodo complemetario a my_split
  std::vector<std::string> &split(const std::string &s, char delim,std::vector<std::string> &elems);
  //Obtiene solo los numeros de un conjunto de cadenas
  std::vector<std::string> get_just_numbrs(std::vector<std::string>lines);
  //Guarda el factor de estructura para Yukawa 3D
  void save_SK_MONO(double *sk, double* k, int Kpoints);
  //Imprime en pantalla un archivo determinado
  void print_file(const std::string &str);
  //Imprime el contenido de un directorio
  int print_directory();
protected:
  std::string file_name;              /**<Nombre del archivo, este ademas incluye la ruta del mismo*/
  bool truncate_file;                /**<Indica si el archivo se tiene que truncar o no. */
  std :: ofstream outfile;          /**<Para escribir archivos.*/
  std :: ifstream infile;          /**<Para leer archivos.*/
  char _date[50];                 /**<Para guardar la fecha actual.*/
  char _directory_name[100];     /*<Para guardar el nombre del directorio*/
  FILE *_directory;             /*<Apuntador a dicho directorio*/
  
};
