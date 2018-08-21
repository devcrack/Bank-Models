
#include "HPP_headers.hpp"


Menu_Class::Menu_Class() {
  
}

void
Menu_Class::print_preface() {
  std::cout<< "Bank Models Version 0.1\n";
  std::cout<< "Contributors:\n";
  std::cout<< "\t -Ezequiel Pedro Ramirez\n";
  std::cout<< "\t -Jesus Benigno Zempeda\n";
  std::cout<< "\t -Magdaleno Medina\n";
  std::cout<< "\t -Patricia Mendoza Mendez\n";
  std::cout<< "Lanimfe Bank Models is a propietary software \n\n";

  std::cout<< "Basic use: \n";
  std::cout<< "./blds [HS1]   [HS2]   [SS]   [Yuk]\n\n";
  
  std::cout<< "For more details, type: ./blds [help]\n";
}


void
Menu_Class::print_main_loop_menu() {
  std::cout<< "Bank Models Version 0.1\n";
  std::cout<< "Contributors:\n";
  std::cout<< "\t -Ezequiel Pedro Ramirez\n";
  std::cout<< "\t -Jesus Benigno Zempeda\n";
  std::cout<< "\t -Magdaleno Medina\n";
  std::cout<< "\t -Patricia Mendoza Mendez\n";
  std::cout<< "Lanimfe Bank Models is a propietary software \n\n";

  std::cout<< "Basic use: \n";
  std::cout<< "./blds [HS1]   [HS2]   [SS]   [Yuk]\n\n";
  
  std::cout<< "For more details, type: ./blds [help]\n";
}




void
Menu_Class::print_Hard_Sphere_options() {
  
}

void
Menu_Class::print_help() {
  system("tput reset");

  std::cout << "About calculates\n";
  
  std::cout << "\t   ./blds [HS1]\n";
  std::cout << "\t\t Executes the process for calculate the Hard Sphere structure factor\n";
  
  std::cout << "\t   ./blds [HS2]\n";
  std::cout << "\t\t Executes the process for calculate the Hard Sphere structure factor but, with a different procedure\n";;
  
  std::cout << "\t   ./blds [SS]\n";
  std::cout << "\t\t Executes the process for calculate the Soft Sphere structure factor\n";
  
  std::cout << "\t   ./blds [Yuk]\n";
  std::cout << "\t\t Executes the process for calculate the Yukawa 3D structure factor\n";

  
  std::cout <<"\n\nAbout printing and plotting\n";
  
  std::cout<< "\t./blds [dir]\n";
  std::cout << "\t\t List all the files avaible of the structures factors\n";

  std::cout<< "\t./blds [plot HS1]\n";
  std::cout << "\t\t Plot the Data of Hard Sphere Structure Factor by the procedure #1\n";

  std::cout<< "\t./blds [plot HS2]\n";
  std::cout << "\t\t Plot the Data of Hard Sphere Structure Factor by the procedure #2\n";

  std::cout<< "\t./blds [plot SS]\n";
  std::cout << "\t\t Plot Data of Soft Sphere Structure factor \n";

  std::cout<< "\t./blds [plot Yuk]\n";
  std::cout << "\t\t Plot Data of Soft Yukawa 3D Structure factor \n";

  std::cout<< "\t./blds [print HS1]\n";
  std::cout << "\t\t Print the data calculated about Structure factor Hard Sphere by the procedure #1\n";

  std::cout<< "\t./blds [print HS2]\n";
  std::cout << "\t\t Print the data calculated about Structure factor Hard Sphere by the procedure #2\n";
  
  std::cout<< "\t./blds [print SS]\n";
  std::cout << "\t\t Print the data calculated about Structure factor Soft Sphere\n";

    std::cout<< "\t./blds [print Yuk]\n";
  std::cout << "\t\t Print the data calculated about Structure factor Yukawa\n";
}


int
Menu_Class :: print_directory(){
  DIR *dir;
  struct dirent *ent;
  std::cout <<"Files of Strucutures Factors:\n\n";
  if ((dir = opendir ("./data")) != NULL) {
  /* print all the files and directories within directory */
    while ((ent = readdir (dir)) != NULL) {
      if(strcmp(ent->d_name,".") && strcmp(ent->d_name,".."))
	printf ("<%s>\n", ent->d_name);
    }
  closedir (dir);
}
  else {
  /* could not open directory */
    std::cout<<"Error to open the data directory\n";
    return EXIT_FAILURE;
  }
}


