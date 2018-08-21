/** @file */ 
#include "HPP_headers.hpp"
/**
 * Esta funcion basicamente es el clasico main.
 * Dicha funcion lo que hace es dirigir le interaccion con el usuario para permitirle 
 * usar los diversas opciones que permite el banco de modelos. En otras palabras, esta funcion
 * lo que hace es presentar un menu para poder operar las distintas funcionalidades que ofrece Bank Models LANIMFE.
 */
int main(int argc, const char *argv[]) {
  int status = 0;
  char option;
  File_Manage * fm = new File_Manage();
  bool flag = false;
  Menu_Class mcls;
  
  std::cout<<"num args ="<<argc<<"\n";
  system("tput reset");
  // This section corresponds to plotting things
  if(argc >= 2){
    if(!(strcmp(*(argv+1), "plot"))) {
      if(argc == 3) {
  	if(!(strcmp(*(argv+2), "HS")))
  	  status = system("xmgrace ./data/sk_HSpheere.dat &");
  	if(!(strcmp(*(argv+2), "SS")))
  	  status = system("xmgrace ./data/sk_SSphere.dat &");
  	if(!(strcmp(*(argv+2), "Yuk"))) 
  	  status = system("xmgrace ./data/sk_MonoYuk.dat &");
      }      
      else	
  	std::cout<<"May be you meant:\n\t./blds plot HS | SS | Yuk\n\n";      
    }
    // This section corresponds to printing
    else {			
      if(!(strcmp(*(argv+1), "print"))) {
  	if(argc == 3) {
  	  if(!(strcmp(*(argv+2), "HS")))
  	    status = system("nano ./data/sk_HSpheere.dat");
  	  if(!(strcmp(*(argv+2), "SS")))
  	    status = system("nano ./data/sk_SSphere.dat");
  	  if(!(strcmp(*(argv+2), "Yuk")))
  	    status = system("nano ./data/sk_MonoYuk.dat");
  	}
  	else
  	  std::cout<<"May be you meant:\n\t./blds print HS1 | HS2 | SS | Yuk\n\n";
      }
    }
    if(!(strcmp(*(argv+1), "help")))
      mcls.print_help();
    if(!(strcmp(*(argv+1), "HS"))) {
      system("./models/01Hard_Spheere");
    }
    if(!(strcmp(*(argv+1), "SS")))  {
      system("./models/02SSphere");      
    }
    if(!(strcmp(*(argv+1), "Yuk"))) {
      system("./models/03Yuk");
    }
    if(!(strcmp(*(argv+1), "dir")))
      mcls.print_directory();
    std::cout <<"\n\nfor more details type:\n\t\t\t./blds -help\n";
  }
  // Start the loop of main menu
  else {
    bool flag_loop = true;
    while(flag_loop){
      mcls.print_preface();
      std::cin.get(option);
      switch (option) {
      case '1':
	flag_loop = loop_simple_menu(1);
	getchar();
	break;
      case '2': 
	flag_loop = loop_simple_menu(2);
	getchar();
	break;
      case '3':
	flag_loop = loop_simple_menu(3);
	getchar();
	break;
      case '4':
	return 0;
	break;
      }
    }
  }
}



/**
 * Presenta un sencillo menu para las cosas relacionadas al factor de estructura de Esfera Dura y,  por obvias razones ademas permite
 * la operatividad de este.
 * Es importante notar que aqui para la ruta de los archivos que manejan  los distintos ejecutables, es decir el ejecutable para Esfera dura(01Hard_Spheere..) esfera suave y demas, toman 
 * la ruta del proceso padre es decir toman como referencia para indexarse en el arbol de directorios la ruta en la que se esta corriendo "blds":
 @param option_sk indica que factor de estructura realizara sus operaciones.
 */
bool loop_simple_menu(int option_sk) {
  char option;
  Menu_Class mcls;

  getchar();
  while(true) {
    switch (option_sk) {	// Acerca de la impresion del menu
    case 1:
      mcls.print_options("Hard Spheere");
      break;
    case 2:
      mcls.print_options("Soft Spheere");
      break;
    case 3:
      mcls.print_options("Yukawa 3D");
      break;
    }
    
    std::cin.get(option);
    switch (option) {
    case '1':			// Acerca de la ejecucion de los factores de estrcutura
      switch (option_sk) {
      case 1:
	system("./models/01Hard_Spheere");
	break;
      case 2:
	system("./models/02SSphere");      
	break;
      case 3:
	system("./models/03Yuk");
	break;
      }
      getchar();
      break;
    case '2':			// Acerca de la impresion de los datos
      switch (option_sk) {
      case 1:
	system("nano ./data/sk_HSpheere.dat");
	break;
      case 2:
	system("nano ./data/sk_SSphere.dat");
	break;
      case 3:
	system("nano ./data/sk_MonoYuk.dat");
	break;
      }
      getchar();
      break;
    case '3':			// Acerca de la graficacion de los datos
      switch (option_sk) {
      case 1:
	system("xmgrace ./data/sk_HSpheere.dat &");
	break;
      case 2:
	system("xmgrace ./data/sk_SSphere.dat &");
	break;
      case 3:
	system("xmgrace ./data/sk_MonoYuk.dat&");
	break;
      }
      getchar();
      break;
    case '4':
      return true;
    case '5':
      return false;;
    }
  }
  return false;
}
	 
