#include "HPP_headers.hpp"

int main(int argc, const char *argv[]) {
  int status = 0;
  char option;
  File_Manage * fm = new File_Manage();
  bool flag = false;
  Menu_Class mcls;
  
  // std::cout<<"num args ="<<argc<<"\n";
  system("tput reset");
  if(argc >= 2){
    if(!(strcmp(*(argv+1), "plot"))) {
      if(argc == 3) {
	if(!(strcmp(*(argv+2), "HS1")))
	  status = system("xmgrace ./data/sk_HSMono_Benny.dat &");
	if(!(strcmp(*(argv+2), "HS2")))
	  status = system("xmgrace ./data/sk_HS_Mono_Patty.dat &");
	if(!(strcmp(*(argv+2), "SS")))
	  status = system("xmgrace ./data/sk_SSphere.dat &");
	if(!(strcmp(*(argv+2), "Yuk")))
	  status = system("xmgrace ./data/sk_MonoYuk.dat &");
	if(status !=0)
	  std::cout << "Plase install grace for this execute the command: \"apt-get install grace\""  << "\n";	
      }
      else 
	std::cout<<"May be you meant:\n\t./blds -plot HS1 | HS2 | SS | Yuk\n\n";
    }
    else {
      if(!(strcmp(*(argv+1), "print"))) {
	if(argc == 3) {
	  if(!(strcmp(*(argv+2), "HS1")))
	    status = system("nano ./data/sk_HSMono_Benny.dat");
	  if(!(strcmp(*(argv+2), "HS2")))
	    status = system("nano ./data/sk_HS_Mono_Patty.dat");
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
    if(!(strcmp(*(argv+1), "HS1"))) {
      std::cout <<"Executing Hard Spheere version #1 Calculates\n";
      system("./models/00hardS_Benny");
    }
    if(!(strcmp(*(argv+1), "HS2"))) {
      std::cout <<"Executing Hard Spheere version #2 Calculates\n";
      system("./models/01hardS_Patty");
    }
    if(!(strcmp(*(argv+1), "SS")))  {
      std::cout <<"Executing Soft Spheere Calculates\n";
      system("./models/02SSphere");      
    }
    if(!(strcmp(*(argv+1), "Yuk"))) {
      std::cout <<"Executing Yukawa 3D Calculates\n";
      system("./models/03Yuk");
    }
    if(!(strcmp(*(argv+1), "dir")))
      mcls.print_directory();
    std::cout <<"\n\nfor more details type:\n\t\t\t./blds -help\n";
  }     
  else
    mcls.print_preface();
  
  // switch(argv[]) {
  // }

  // usleep(1800000);
  // system("tput reset");
  // print_menu();    
  
  // while(true) {
  //   std::cin.get(option);
  //   switch(option) {
  //   case 48:
  //     system("tput reset");
  //     std::cout<<"\n\n\t\tSee you next time...\n\n";
  //     usleep(2000000);
  //     system("tput reset");
  //     return 0;
  //   case 49:      
  //     status = system("./models/00hardS_Benny");
  //     std::cout << "\nStatus..."  << status << "\n";
  //     flag = true;
  //     getchar();
  //     break;
  //   case 50:
  //     status = system("./models/01hardS_Patty");
  //     std::cout << "\nStatus..."  << status << "\n";
  //     flag = true;
  //     getchar();
  //     break;
  //   case 51:
  //     status = system("./models/02SSphere");
  //     std::cout << "\nStatus..."  << status << "\n";
  //     flag = true;
  //     getchar();
  //     break;
  //   case 52:
  //     status = system("./models/03Yuk");
  //     std::cout << "\nStatus..."  << status << "\n";
  //     flag = true;
  //     getchar();
  //     break;    
  //   case 97:
  //   case 67: 			// Aa
  //     fm->print_file("./data/sk_HSMono_Benny.dat");
  //     std::cout << "\ndata sk HSMono_Bnny"  << "\n";
  //     flag = true;       
  //     getchar();
  //     break;
  //   case 98:
  //   case 68:			// Bb
  //     fm->print_file("./data/sk_HS_Mono_Patty.dat");
  //     std::cout << "\ndata sk HSMono_Ptty"  << "\n";
  //     flag = true;
  //     getchar();
  //     break;
  //   case 99:
  //   case 69:      		// Cc
  //     fm->print_file("./data/sk_SSphere.dat");
  //     std::cout << "\ndata sk SSPhere"  << "\n";
  //     flag = true;
  //     getchar();
  //     break;
  //   case 100:
  //   case 70:      		// Dd
  //     fm->print_file("./data/sk_MonoYuk.dat");
  //     std::cout << "\ndata sk Yukawa"  << "\n";
  //     flag = true;
  //     getchar();      
  //     break;
  //   case 101:
  //   case 71:      		// Ee
  //     status = system("xmgrace ./data/sk_HSMono_Benny.dat");
  //     if(status !=0)
  // 	std::cout << "Plase install grace for this execute the command: \"apt-get install grace\""  << "\n";
  //     std::cout << "\nStatus..."  << status << "\n";
  //     flag = true;
  //     getchar();      
  //     break;
  //   case 102:
  //   case 72:      		// Ff
  //     status = system("xmgrace ./data/sk_HS_Mono_Patty.dat");
  //     if(status !=0)
  // 	std::cout << "Plase install grace for this execute the command: \"apt-get install grace\""  << "\n";
  //     std::cout << "\nStatus..."  << status << "\n";
  //     flag = true;
  //     getchar();      
  //     break;
  //   case 103:
  //   case 73:      		// Gg
  //     status = system("xmgrace ./data/sk_SSphere.dat");
  //     if(status !=0)
  // 	std::cout << "Plase install grace for this execute the command: \"apt-get install grace\""  << "\n";
  //     std::cout << "\nStatus..."  << status << "\n";
  //     flag = true;
  //     getchar();      
  //     break;
  //   case 104:
  //   case 74:      		// Hh 
  //     status = system("xmgrace ./data/sk_MonoYuk.dat");
  //     if(status !=0)
  // 	std::cout << "Plase install grace for this execute the command: \"apt-get install grace\""  << "\n";
  //     std::cout << "\nStatus..."  << status << "\n";
  //     flag = true;
  //     getchar();      
  //     break;
  //   }    
  //   if(flag) {
  //     std::cout << "\nPress Enter to Continue..."  << "\n";
  //     getchar();
  //     flag = false;
  //   }
  //   system("tput reset");
  //   print_menu();
  // }
  
}


void print_menu(){
  std::cout<<"\t Select an option please\n\n";  
  std::cout<<"\t 1).-SK_Hard_Spheere_Bnny\n";
  std::cout<<"\t a).-List data SK_Hard_Spheere_Bnny\n";
  std::cout<<"\t e).-Plot SK_Hard_Sphere_Bnny\n\n";
  
  std::cout<<"\t 2).-SK_Hard_Spheere_Ptty\n";
  std::cout<<"\t b).-List data SK_Hard_Spheere_Ptty\n";
  std::cout<<"\t f).-Plot SK_Hard_Sphere_Ptty\n\n";
  
  std::cout<<"\t 3).-SK_Soft_Spheere\n";
  std::cout<<"\t c).-List data SK_Soft_Spheere\n";
  std::cout<<"\t g).-Plot SK_Soft_Spheere\n\n";
  
  std::cout<<"\t 4).-SK_Yukawa_Spheere\n";
  std::cout<<"\t d).-List data SK_Yukawa\n";
  std::cout<<"\t h).-Plot SK_Yukawa\n\n";
  
  std::cout<<"\t 0).-EXIT\n"; 
}
  
		 

