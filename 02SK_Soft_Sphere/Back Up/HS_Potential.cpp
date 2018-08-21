#include "HS_Potential.hpp"
/*
 *Simple Constructor
 */
HS_Potential::HS_Potential():
  _kMax(0.0),  
  _deltaK(0.0),
  _phi(0.0),
  _nKas(0),
  _SK_max(0.0),
  _SK_min(0.0)
{

}


/*
 *Loads the setup of this structure factor from the options.txt file
 */
void
HS_Potential::set_up_from_config_file(std::vector<std::string>_lines) {
  this->_kMax   = atof(_lines[0].c_str());
  this->_deltaK = atof(_lines[1].c_str());
  this->_phi    = atof(_lines[2].c_str());
  this->_nKas = (int)(this->_kMax / this->_deltaK);
  this->_sk = new double[this->_nKas];
}
  

void
HS_Potential::Structure_Factor() {
  double skMax;
  int index = 1;

  this->_sk[0] = skMax = this->KZero(this->_phi);
  for(double i = this->_deltaK; i <= this->_kMax; i+= this->_deltaK, ++index) {
    this->_sk[index] = this->sK(this->_phi, i);
    if(this->_sk[index] > skMax)
      skMax = this->_sk[index];
  }
}



double
HS_Potential::sK(double phi, double k){//Calculo de S(K)
  double
    sk,ck;

  k = k * pow((phi - (phi*phi)/16.0)/phi,1/3.0);
  phi = phi - (phi*phi)/16.0;
  ck = (-24.0 * phi / (pow(k,6) * pow(1-phi,4))) * (PZero(phi,k) + P_One(phi,k)*k*sin(k) + P_Two(phi,k)*cos(k));
  sk = 1/(1-ck);
  
  return sk;
}


double
HS_Potential::PZero(double phi, double k) {//Calculo de P0
  double p = 0;
  
  p = 3*phi*(pow(k,2)*pow(2+phi,2) + 4*pow(1+2*phi,2));
  
  return p;
}

double
HS_Potential::P_One(double phi, double k) {//Calculo de P1
  double p = 0;
  
  p = -12*phi*pow(1+2*phi,2) +k*k*(1-6*phi + 5*pow(phi,3));
  
  return p;
}

double
HS_Potential::P_Two(double phi, double k) {//Calculo de P2
  double p = 0;
  
  p = -12*phi*pow(1 + 2*phi,2) + 3*phi*k*k*(-2+4*phi+7*phi*phi) - (pow(k,4)/2.0)*(2-3*phi+pow(phi,3));
  
  return p;
}


double
HS_Potential:: KZero(double phi) {
  double ck, sk;
  
  phi = phi - (phi*phi)/16.0;
  ck = phi*(phi*phi*phi - 4*phi*phi + 2*phi -8)/pow(1-phi,4);
  sk = 1/(1.0-ck);

  return sk;
}

