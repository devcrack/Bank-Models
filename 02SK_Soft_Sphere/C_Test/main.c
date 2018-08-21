#include "scgle.h"



int main() {
  sumTiempo = 0;
  inputOptions = fopen("options.txt","r");
  data = nuevaSCGLE();
  
  allocateMemorry(&data,inputOptions); /* Yep */
  t = data->temI;		       /* Yep */
  phi = data->phi;		       /* Some thing like that */  
  deltaTemp = (data->temI-data->temF)/data->ts; /* Discarded */
  data->dss = blip(t);				/* Yepi */  
  data->phi = phi*pow(data->dss,3); /* Aqui desmadran Phi */
  /* printf("Dss: %f, phi:%f\n",data->dss,data->phi); */
  data->ski = factorDeEstructura(data,"factor.dat");





  /* Estos desvergues despues */
  /* printf("Factor de Estructura iteraciones = %d\n", data->kas); */
  /* for (int iter = 0; iter < data->kas; iter++) { */
  /*   printf("SK = %.17f\n",*(data->ski+iter));  */
  /* } */

  /* sprintf(nom,"%s/skmax%d.dat",date,0); */
  /* skmax = fopen(nom, "w"); */
  /* sprintf(nom,"%s/grmax%d.dat",date,0); */
  /* grmax = fopen(nom, "w"); */
  /* t-=deltaTemp; */
  /* data->dss = blip(t);//\**** */
  /* data->phi = phi*pow(data->dss,3); */
  /* printf("Dss: %f, phi:%f\n",data->dss,data->phi); */
  /* sprintf(nom,"factor%d.dat",0); */
  /* data->skf = factorDeEstructura(data,"factor.dat"); */
  /* data->sk = (double*)malloc(sizeof(double)*data->kas); */
  /* data->phi = phi; */
  return 0;  
}
