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
  void save_sk(double *sk, double* k, int Kpoints, int species); 
  bool read_Sk();
  void read_file_line_by_line_p();
  std::vector<std::string> read_file_line_by_line_g(std::string);
  void set_date_FILE_from_Options_file();
  void create_directory(const char *dir_name[]);
  void Add_a_File_to_Directory_Date(const char dir_name[]);
  void save_Sk_SSPhere(int nKas,double phi,double deltaK,double _dss ,double *sk);
  void save_Sk_HSphere(double kMax, double deltaK, double *sk);
  char* get_date();
protected:
  std::string file_name;
  bool truncate_file;
  std :: ofstream outfile;
  std :: ifstream infile;
  char _date[50];
  char _directory_name[100];
  FILE *_directory;
  
};
