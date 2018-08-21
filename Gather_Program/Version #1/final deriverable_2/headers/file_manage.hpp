#include<fstream> 		// Class for use files
#include<string>
#include<iostream>
#include<iomanip>
#include<limits>
#include <vector>
#include<cstdlib>
#include <sys/stat.h>

class File_Manage {
  
public:
  File_Manage();
  File_Manage(const std::string& str, bool truncate);
  ~File_Manage();
  void save_sk(double *sk, double* k, int Kpoints, int species); 
  bool read_Sk();
  void read_file_line_by_line_p();
  std::vector<std::string> read_file_line_by_line_g(std::string);
  void set_date_FILE_from_Options_file();
  void create_directory(const char *dir_name[]);
  void Add_a_File_to_Directory_Date(const char dir_name[]);
  void save_Sk_SSPhere(int nKas,double phi,double deltaK,double _dss ,double *sk);
  char* get_date();
  std::vector<std::string> my_split(const std::string &str, char delim);
  std::vector<std::string> &split(const std::string &s, char delim,std::vector<std::string> &elems);
  std::vector<std::string> get_just_numbrs(std::vector<std::string>lines);
  void save_SK_MONO(double *sk, double* k, int Kpoints);
  void print_file(const std::string &str);
  int print_directory();
protected:
  std::string file_name;
  bool truncate_file;
  std :: ofstream outfile;
  std :: ifstream infile;
  char _date[50];
  char _directory_name[100];
  FILE *_directory;
  
};
