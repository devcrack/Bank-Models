#include "Common_Potential.hpp"
/*
 *Simple Constructor
 */
Common_Potential::Common_Potential():
  _kMax(0.0),  
  deltaK_dk(1.0e-1),
  _phi_eta(0.5e0),
  nKas_Kpoints(199),
  _SK_max(0.0),
  _SK_min(0.0),
  _Species(1),
  _SDimen(3),
  _sk (new double[199]),
  ck  (new double[199]),
  _k   (new double[199])
{
}


Common_Potential::Common_Potential(double vol_fract):
  _kMax(0.0),  
  deltaK_dk(0.0),
  _phi_eta(vol_fract),
  nKas_Kpoints(0),
  _SK_max(0.0),
  _SK_min(0.0)
{
  this->_kMax = 0;
  this->deltaK_dk  = 1.0e-1;
  this->nKas_Kpoints = 199;
  this->_Species = 1;
  this->_SDimen = 3;
  this->_sk = new double[199];
  this->ck  = new double[199];
  this->_k =   new double[199];  
}

Common_Potential::Common_Potential(int SDimen, int Species):
  _kMax(0.0),  
  deltaK_dk(0.0),
  _phi_eta(0.0),
  nKas_Kpoints(0),
  _SK_max(0.0),
  _SK_min(0.0),
  _SDimen(SDimen),
  _Species(Species)
{

}




Common_Potential::~Common_Potential() {
  delete [] _sk;
  delete [] _k;
  delete [] ck; 
}


/*
 *Loads the setup of this structure factor from the options.txt file
 */
void
Common_Potential :: set_up_from_config_file(std::vector<std::string>_lines) {
  this->_kMax   = atof(_lines[0].c_str());
  this->deltaK_dk = atof(_lines[1].c_str());
  this->_phi_eta    = atof(_lines[2].c_str());
  this->nKas_Kpoints = (int)(this->_kMax / this->deltaK_dk);
  this->_sk = new double[this->nKas_Kpoints];
  this->ck =  new double[this->nKas_Kpoints];
}



void
Common_Potential :: Structure_Factor() {
  double skMax;
  int index = 1;

  this->_sk[0] = skMax = this->KZero(this->_phi_eta);
  for(double i = this->deltaK_dk; i <= this->_kMax; i+= this->deltaK_dk, ++index) {
    this->_sk[index] = this->sK(this->_phi_eta, i);
    if(this->_sk[index] > skMax)
      skMax = this->_sk[index];
  }
}



double
Common_Potential::sK(double phi, double k){//Calculo de S(K)
  double
    sk,ck;

  k = k * pow((phi - (phi*phi)/16.0)/phi,1/3.0);
  phi = phi - (phi*phi)/16.0;
  ck = (-24.0 * phi / (pow(k,6) * pow(1-phi,4))) * (PZero(phi,k) + P_One(phi,k)*k*sin(k) + P_Two(phi,k)*cos(k));
  sk = 1/(1-ck);
  
  return sk;
}


double
Common_Potential::PZero(double phi, double k) {//Calculo de P0
  double p = 0;
  
  p = 3*phi*(pow(k,2)*pow(2+phi,2) + 4*pow(1+2*phi,2));
  
  return p;
}

double
Common_Potential :: P_One(double phi, double k) {//Calculo de P1
  double p = 0;
  
  p = -12*phi*pow(1+2*phi,2) +k*k*(1-6*phi + 5*pow(phi,3));
  
  return p;
}


double
Common_Potential :: P_Two(double phi, double k) {//Calculo de P2
  double p = 0;
  
  p = -12*phi*pow(1 + 2*phi,2) + 3*phi*k*k*(-2+4*phi+7*phi*phi) - (pow(k,4)/2.0)*(2-3*phi+pow(phi,3));
  
  return p;
}


double
Common_Potential:: KZero(double phi) {
  double ck, sk;
  
  phi = phi - (phi*phi)/16.0;
  ck = phi*(phi*phi*phi - 4*phi*phi + 2*phi -8)/pow(1-phi,4);
  sk = 1/(1.0-ck);

  return sk;
}

void
Common_Potential :: calc_static_K() {
  for (int i = 0; i <  this->nKas_Kpoints ; i++) {
    this->_k[i] = (i+1) * this->deltaK_dk;
  }
}




void 
Common_Potential :: calc_Sk_hs_py_mono( bool VW_op) {
  double *_k_vw;
  double
    eta_vw,
    _calc_1,
    _calc_2,
    _calc_3,
    _calc_4,
    _calc_sin,
    _calc_cos,
    tmp_aux;
  // Param checks
  if(this->_SDimen != 3 || this->_Species != 1) {
    printf("Error, cannot compute the system direct Correlaction function");      

    return;
  }

  _k_vw = new double[this->nKas_Kpoints];
  if(VW_op) {
    eta_vw = this->_phi_eta * (1.e0 - this->_phi_eta / 16.e0);
    for (int i = 0 ; i < this->nKas_Kpoints; i++) 
      _k_vw[i] = this->_k[i] * pow(eta_vw / this->_phi_eta, 1.e0 / 3.e0);
  }
  else {
    eta_vw = this->_phi_eta;
    for (int i = 0; i < this->nKas_Kpoints; i++) 
      _k_vw[i] = this->_k[i];	
  }
  _calc_4 = pow(1.0 - eta_vw, 4);
  _calc_1 = -pow(1.0 + 2.0 * eta_vw,2) / _calc_4;
  _calc_2 = ( 6.0 *eta_vw * pow(1.0 + eta_vw / 2.0, 2)  )  / _calc_4;
  _calc_3 = ( -eta_vw * pow(1.0 + 2.0 * eta_vw,2) / 2.0 ) / _calc_4;

  for (int i = 0; i < this->nKas_Kpoints; i++) {
    _calc_sin = sin(_k_vw[i]);
    _calc_cos = cos(_k_vw[i]);
    ck[i] =  (_calc_1 * (_calc_sin - _k_vw[i] * _calc_cos) / pow(_k_vw[i],2.e0))      
      +
      (_calc_2 *( ( 2.e0 * _k_vw[i]* _calc_sin ) + ( ( -pow(_k_vw[i], 2.e0) + 2.e0) * _calc_cos) - 2.e0) / pow(_k_vw[i],3.e0))
      +
      ( _calc_3
	*
	(
	 (4.e0 * pow(_k_vw[i], 3.e0) - 24.0 * _k_vw[i]) * _calc_sin
	 +
	 ( (-pow(_k_vw[i], 4.e0) + 12.e0 * pow(_k_vw[i], 2.e0) - 24.e0) * _calc_cos)
	 +
	 24.e0
	 )
	/ pow(_k_vw[i], 5)
       );
    ck[i] = 24.e0 * eta_vw * ck[i] / _k_vw[i];
    _sk[i] = 1.e0 / (1.e0 - ck[i]);
  }
}


double*
Common_Potential :: get_SK() { return this->_sk; }

double*
Common_Potential :: get_K() { return this->_k; }

int
Common_Potential::get_KPoints() { return this->nKas_Kpoints;}
  


