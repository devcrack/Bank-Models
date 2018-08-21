


/** @file */ 
#include "HPP_headers.hpp"
/**
 * Esta funcion basicamente es el clasico main.
 * Dicha funcion lo que hace es dirigir le interaccion con el usuario para permitirle 
 * usar los diversas opciones que permite el banco de modelos. En otras palabras, esta funcion
 * lo que hace es presentar un menu para poder operar las distintas funcionalidades que ofrece Bank Models LANIMFE.
 */
int main(int argc, char *argv[]) {
  Menu_Class mcls;
  
  if(argc > 1)
    if(!(strcmp(*(argv + 1),"--help")))
      mcls.print_help();
    else
      process_args(argc, argv);
  else {
        
    bool flag_loop = true;
    char option;
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
  std::string temperature_ini, vol_frac, exec;
  char dmb_char;
  
  getchar();
  while(true) {
    switch (option_sk) {	// Acerca de la impresion del menu
    case 1:			// Imprime opciones para esfera dura
      mcls.print_options("Hard Spheere");
      break;
    case 2:			// Imprime opciones para esfera suave
      mcls.print_options("Soft Spheere");
      break;
    case 3:			// Imprime opciones para Yukawa
      mcls.print_options("Yukawa 3D");
      break;
    }
    
    std::cin.get(option);
    switch (option) {
    case '1':			// Acerca de la ejecucion de los factores de estrcutura
      switch (option_sk) {
      case 1:			// Calculos de Esfera Dura
	{
	  getchar();
	  vol_frac = get_value("Type the Fraction Volume :");	
	  if(vol_frac.empty()) {
	    std::cout<< "Default value for the fraction volumen have stablished \n";
	    exec = "./models/01Hard_Spheere";
	  }
	  else {
	    std::cout<< "Volumen Fraction = " << vol_frac <<"\n";
	    exec = "./models/01Hard_Spheere " + vol_frac;
	  }
	  system(exec.c_str());
	  std::cout << "\n\nPress Enter to Continue\n";
	  break;
	}
      case 2:			// Calculos de Esfera Suave
	{
	  getchar();
	  vol_frac = get_value("Type the Fraction Volume: ");	
	  temperature_ini = get_value("Type the Intial temperature: ");	
	  if(vol_frac.empty() && temperature_ini.empty()) {
	    std::cout<< "Default values for the fraction volumen and Initial Temperature have stablished \n";
	    exec = "./models/02SSphere";
	  }
	  else {
	    if(vol_frac.empty())
	      vol_frac = "0.56";
	    if(temperature_ini.empty())
	      temperature_ini =  "0.1";
	    exec = "./models/02SSphere " + vol_frac + " " + temperature_ini;
	  }
	  system(exec.c_str());
	  std::cout << "\n\nPress Enter to Continue\n";
	}
	break;
      case 3:			// Calculos de Yukawa 3D
	getchar();
	vol_frac = get_value("Type the Fraction Volume: ");	
	temperature_ini = get_value("Type the Intial temperature: ");	
	if(vol_frac.empty() && temperature_ini.empty()) {
	  std::cout<< "Default values for the fraction volumen and Initial Temperature have stablished \n";
	  exec = "./models/03Yuk";	  
	}
	else {
	  if(vol_frac.empty())
	    vol_frac = "0.5";
	  if(temperature_ini.empty())
	    temperature_ini =  "1.e0";
	  exec = "./models/03Yuk "  + vol_frac + " " + temperature_ini;
	}	  	
	system(exec.c_str());
	std::cout << "\n\nPress Enter to Continue\n";
	break;
      }
      getchar();
      break;
    case '2':			// Acerca de la impresion de los datos
      switch (option_sk) {
      case 1:			// Imprime datos esfera Dura
	system("nano ./data/sk_HSpheere.dat");
	break;
      case 2:			// Imprime datos Esfera Suave
	system("nano ./data/sk_SSphere.dat");
	break;
      case 3:			// Imprime datos Yukawa
	system("nano ./data/sk_MonoYuk.dat");
	break;
      }
      getchar();
      break;
    case '3':			// Acerca de la graficacion de los datos
      switch (option_sk) {
      case 1:			// Grafica Datos Esfera dura 
	system("xmgrace ./data/sk_HSpheere.dat &");
	break;
      case 2:
	system("xmgrace ./data/sk_SSphere.dat &");
	break;
      case 3:
	system("xmgrace ./data/sk_MonoYuk.dat &");
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
	 

std::string get_value(const std::string msg) {
  char dmb_char;
  std::string vol_frac;
  
  std::cout << msg;
  dmb_char = getchar();
  
  if(dmb_char == '\n') 
    vol_frac="\0";
  else
    std::getline(std::cin, vol_frac,'\n');
  return vol_frac;
}


void process_args(int argc, char* argv[]) {

  std::string temperature_ini = "\0";
  std::string vol_frac        = "\0";
  bool exe = false;
  
  // const char* const short_opts = "hdsypt:vo";
  const char* const short_opts = "hdfypt:v:ox";
  const option long_opts[] = {
    {"help",  no_argument,             0, 'l'},
    {"hs",    no_argument,             0, 'h'},
    {"ss",    optional_argument,       0, 'f'},
    {"yuk",   optional_argument,       0, 'y'},
    {"vf",   required_argument,        0, 'v'},
    {"ti",   required_argument,        0, 't'},
    {"print", no_argument,             0, 'p'},
    {"plot",  no_argument,             0, 'o'},
    {"exe",  no_argument,             0, 'x'},
    {0, 0, 0, 0}  
  };
  int option;
  std::list<int> actions;
  std::list<int> list_acts;
  std::list<int>::iterator it;
  
  while(( option = getopt_long(argc, argv, short_opts, long_opts, 0 )) // Si ninguna opcion fue ingresada por obvias razones no entra al loop de evaluacion de argumentos
	!= -1 )
    {
      switch(option)
	{
	  // Parametros inciales para los calculos
	case'v' : 		// Se carga la fraccion de Volumen
	  if(optarg)
	    vol_frac = optarg;
	  else 
	    std::cout << "You must correct the option" << std::endl;
	  // std::cout << "volume fraction typed " << vol_frac << std::endl;
	  break;	  
	case 't': 		// Se carga la temperatura Inicial
	  if(optarg)
	    temperature_ini = optarg;
	  else
	    std::cout << "You must correct the option" << std::endl;
	  // std::cout << "intial temperature typed " << temperature_ini << std::endl;
	  break;
	
	case 'h':		// Potencial Esfera dura
	  actions.push_back(1);
	  // potential = "./models/01Hard_Spheere";
	  // std::cout << "You hit for Hard Sphere Potential" << std::endl;
	  break;
	case 'f':		// Potencial Esfera Suave 
	  actions.push_back(2);	  
	  break;
	case 'y':		// Potencial de Yukawa
	  actions.push_back(3);
	  break;	
	case 'o':		// Graficar
	  actions.push_back(4);
	  std::cout << "Plot Option" << std::endl;	  
	  break;
	case 'l':		// Imprimir Ayuda
	  std::cout << "You hit help" << std::endl;
	  actions.push_back(6);
	  break;
	case 'p':		// Imprime la ayuda
	  std::cout << "Print Option" << std::endl;
	  actions.push_back(5);
	  break;
	case 'x':
	    exe = true;
	  break;
	}
    }

  // Ejecutar comando
  for (it = actions.begin(); it != actions.end(); ++it) // Ordena lista de acciones
    inserts_ordered(*it,list_acts);
  exeute_command(list_acts, exe, vol_frac, temperature_ini);
}


// Inserta ordenado para la lista de acciones
void inserts_ordered(int value, std::list<int> &lst_actns) {
  std::list<int>::iterator it;
  std::list<int>::iterator pos_insert;;
  bool stop_loop = false;
    
  for (it=lst_actns.begin(); it!=lst_actns.end() && stop_loop == false; ++it) {
    if(value < *it) {
      stop_loop = true;
    }
  }
  if(!stop_loop) 
    lst_actns.push_back(value);
  else
    lst_actns.insert(--it,value);
}


void exeute_command(std::list<int> &lst_actns, bool exe, std::string &vf, std::string &ti) {
  std::string exec;
  std::list<int>::iterator it;
  bool break_loop = false;
  
  for (it=lst_actns.begin(); it!=lst_actns.end() && break_loop == false; ++it) {
    switch (*it) {
    case 1:
      {				// Ejecutar calculo de esfera dura
	if(exe)
	  {
	    // std::cout << "Ejecutando Calculo Esfera Dura\n";
	    exec = "./models/01Hard_Spheere" ;
	    if(vf.empty()) 
	      std::cout<< "Default value for the fraction volumen have stablished \n";
	    else
	      exec = exec  + " " + vf;
	    system(exec.c_str());
	    exec = "\0";      	   
	  }
	break;
      }
    case 2:
      {				// Ejecutar calculo de Esfera Suave
	if(exe)
	  {
	    // std::cout << "Ejecutando Calculo Esfera Dura\n";
	    exec = "./models/02SSphere";
	    if(vf.empty() && ti.empty())
	      std::cout<< "Default values for the fraction volumen and Initial Temperature have stablished \n";
	    else {
	      if(vf.empty())
		vf = "0.56";
	      if(ti.empty())
		ti =  "0.1";
	      exec = exec + " " + vf + " " + ti;
	    }
	  }
	break;
      }
    case 3:
      {				// Ejecutar Calculo de Yukawa
	if(exe)
	  {
	    // std::cout << "Ejecutando Calculo Yukawa\n";
	    exec = "./models/03Yuk";
	    if(vf.empty() && ti.empty())
	      std::cout<< "Default values for the fraction volumen and Initial Temperature have stablished \n";
	    else {
	      if(vf.empty())
		vf = "0.5";
	      if(ti.empty())
		ti = "1.e0";
	      exec = exec + " " +  vf + " " + ti;
	    }
	  }
	break;
      }
    case 4:			// GRAFICACION
      {
	switch(lst_actns.front())
	  {
	  case 1:			// Grafica Datos esfera Dura
	    exec =  "xmgrace ./data/sk_HSpheere.dat &";
	    break;
	  case 2:			// Grafica Datos Esfera Suave
	    exec = "xmgrace ./data/sk_SSphere.dat &";
	    break;
	  case 3:			// Grafica Datos Yukawa
	    exec = "xmgrace ./data/sk_MonoYuk.dat &";
	    break;
	  }
	break_loop = true;
	break;
      }
    case 5:			// IMPRESION
      {
	switch(lst_actns.front()) {
	case 1:		// Imprime datos esfera dura
	  exec = "nano ./data/sk_HSpheere.dat";
	  break;
	case 2:			// Imprime datos Esfera suave
	  exec = "nano ./data/sk_SSphere.dat";
	  break;
	case 3:			// Imprime datos Yukawa
	  exec = "nano ./data/sk_MonoYuk.dat";
	  break;
	}
	break_loop = true;
	break;
      }
    }
    system(exec.c_str());
    exec = "\0";
  }
}




// Snippet code from : https://stackoverflow.com/questions/8809196/how-do-i-use-getopt-long-to-parse-multiple-arguments
// for (int i = optind; i < argc; i++) { // Resuelve el hecho de como manejar  mas de un argumento por opcion;
// 	std::cout << "non-option arg: " << argv[i] << std::endl;
// }  
