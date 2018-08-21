#ifndef HARD_SHPERE_DYN
#define HARD_SHPERE_DYN

#include "hard_sphere_mono.hpp"
#include<iostream>
#include<file_manage.hpp>

class Hard_Sphere_Dynamic {
  /**@brief Variables del sistema

     Permite tener todas las variables del sistema encapsuladas en este modula, notar que esto no es optimo
     para los fines de este programa.
     @Author's: benigno zempeda & delcranck
     @Date: 28/08/2017
  */
public :  
  /**
   *Se inicializan todas las varibles del sistema
  */
  Hard_Sphere_Dynamic(bool VW_op, bool Gamm_op, bool SK_writting_op);
  Hard_Sphere_Dynamic(bool VW_op, bool Gamm_op, bool SK_writting_op, double vol_fraction);
  ~Hard_Sphere_Dynamic();

  void print_kpionts_value();

  const int* get_ik_test();
  void print_ik_test();
  void start();
  double get_volumen_fraction();
  
  // const S* getArrayPointer() const { return &myArray; }
protected:
  bool _VW_op;                           /**<Determina si esto tiene que escribirse en archivo, esto tiene que reescribirse, no es optimo */
  bool _Gamm_op;                 
  bool _SK_writting_op;        
  Hard_Sphere_Mono _mono_hard_sphere;  /**<Se encapsulan todas las operaciones para el factor de estructura*/
  File_Manage _file;                  /**<Permite tener un acceso rapido a las operaciones con archivos*/
};

#endif
