#include "system_variables.hpp"


/**
 *Initializes all variables used in this class
 * This is the subroutine Sys_Variables_Alloc in System_Variables_Module file variables_modules.f90
 */
Sys_Variables :: Sys_Variables(int _Ialloc)
  :
  _Species(_Ialloc),
  _eta(new double[_Ialloc]),
  _sigma(new double[_Ialloc]),
  _SDimen(3),
  m_rho(new double[_Ialloc * _Ialloc]),
  m_rho_i(new double[_Ialloc * _Ialloc]),
  d0m(new double[_Ialloc * _Ialloc])
{
  this->_sigma[0]    = distu;
  this->_eta[0]      = 0.3e0;	// esto es phi
}


Sys_Variables :: Sys_Variables(int _Ialloc, double _vol_fraction)
  :
  _Species(_Ialloc),
  _eta(new double[_Ialloc]),
  _sigma(new double[_Ialloc]),
  _SDimen(3),
  m_rho(new double[_Ialloc * _Ialloc]),
  m_rho_i(new double[_Ialloc * _Ialloc]),
  d0m(new double[_Ialloc * _Ialloc])
{
  this->_sigma[0]    = distu;
  this->_eta[0]      = _vol_fraction;	// esto es phi
}


/**
 * Here do not have any argument cause in the fortran program this is : 
 *                                                                      call D0M_ini(species,sigma)
 * and hence we can see that species and sigma are member variable of this class.
 */
void
Sys_Variables :: D0M_ini() {
  for (int i1 = 0; i1 < this->_Species; i1++)
    for (int i2 = 0; i2 < _Species; i2++) {
      if(i1 == i2) 
	*(d0m + i1 *_Species + i2)   = 1.0 / this->_sigma[i1];    
      else  
	*(d0m + i1 *_Species + i2)   = 0.0;	
    }
}

void
Sys_Variables :: rho_ini() {  
  for (int i1 = 0; i1 < this->_Species; i1++) 
    for (int i2 = 0; i2 < this->_Species; i2++) {
      if(i1 == i2) {
	if(_SDimen == 3) 
	  *(m_rho+i1*_Species+i2) = 6.0 * this->_eta[i1] / (pi * ( pow(this->_sigma[i1],3)));
	else if(_SDimen == 2) 
	  *(m_rho+i1*_Species+i2) = 4.0 * this->_eta[i1] / (pi * ( pow(this->_sigma[i1],2)));
	*(m_rho   + i1 *_Species + i2)   = sqrt(*(m_rho+i1*_Species+i2));
	*(m_rho_i + i1 *_Species + i2)   = 1.0 / *(m_rho_i+i1*_Species+i2);
      }
      else {
	*(m_rho+i1*_Species+i2)   = 0.0;
	*(m_rho_i+i1*_Species+i2) = 0.0;
      }
    }
}

void
Sys_Variables :: print_rho(){
  std::cout<<"Imprimiendo valores de m_rho\n";
  for (int i1 = 0; i1 < this->_Species; i1++) 
    for (int i2 = 0; i2 < this->_Species; i2++) 
      std::cout << std::fixed << std::setprecision(17) << *(m_rho+i1*this->_Species+i2) <<"\n";
}
      

Sys_Variables ::~Sys_Variables () {
  delete[] m_rho;
  delete[] m_rho_i;
  delete[] d0m;
}

//Teoria de los liquidos Simples
//Theory of simple liquids
//MccDonald 
//Sergay 
