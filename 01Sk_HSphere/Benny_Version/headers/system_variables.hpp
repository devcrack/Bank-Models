#ifndef SYS_VAR
#define SYS_VAR

#include <math.h>
#include<iostream>
#include<stdio.h>
#include<iomanip> 		// Include this header for that you can use the std::limits :P
const double pi =  4.0 * atan(1.0);;
const double  distu = 1.e0;

class Sys_Variables {
  /**@brief clase para manejar las variables del sistema.
     
     @author Bennigno zempeda & delcranck
     @date 28/08/2017
  */
  
 public :
  /**
   *Inicializacion de todas las variables del sistema.
  */  
  explicit Sys_Variables (int _Ialloc);

  /**
   *Este constructor incluye la fraccion de Volumen.
   */
  Sys_Variables(int _Ialloc, double _vol_fraction);
  
  
  /**
   *Destructor libera la memoria de las variables dinamicas usadas.
   */
  ~Sys_Variables ();

  /**
   *Initialization of D0M variable
   * Here do not have any argument cause in the fortran program this is : call D0M_ini(species,sigma),
   * and hence we can see that species and sigma are member variable of this class.
   */
  void D0M_ini();
  void rho_ini ();

  void print_rho();
 protected:

  
  int _SDimen;                         /**<Dimension del espacio */
  int _Species;                       /**<NUmero de especies usadas en el sistema */
  double *_eta;                      /**<Densidad de la dimension de cada especia, en otras palabras la fraciones de volumen*/
  double *_sigma;                   /**<Dimensioless characteristc distance of each species, Es la dimension de distanciao unidades de distancia,*/
  double *m_rho;
  double *m_rho_i;
  double *d0m;
  
};  
  

#endif
