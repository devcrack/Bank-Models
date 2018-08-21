#include "scgle.h"
#include <math.h>

/* Constructor de scgle*/
scgle *nuevaSCGLE(){
  scgle *nva;
  
  nva = (scgle*)malloc(sizeof(scgle));
  nva->kMax = 0;
  nva->deltaK = 0;
  nva->phi = 0;
  nva->deltaT = 0.0000001;
  nva->skMax = 0;
  nva->skMin = 0;
  nva->mu = 2;
  nva->alfa = 1.305;
  nva->d = 1;
  nva->b = 0;
  nva->deltaZ = (double*)malloc(sizeof(double)*(TMMAXT));
  nva->uMax = 2;
  nva->arrestado = false;

  return nva;
}

/* Constructor de scgle */
void allocateMemorry(scgle **data, FILE *inputOptions){
  scgle *aux = *data;
  char *input,directory[100], ops[500];
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  FILE *options;
  
  input = (char*)malloc(sizeof(char)*200);
  fscanf(inputOptions, "%s", input);
  aux->kMax = atof(input);	/* L1 */
  /* printf("Kmax:%f\n", aux->kMax); */
  fscanf(inputOptions, "%s", input);
  aux->deltaK = atof(input);	/* L2 */
  /* printf("deltaK:%.17f\n", aux->deltaK); */
  fscanf(inputOptions, "%s", input);
  aux->phi = atof(input);	/* L3 CHIDA*/
  fscanf(inputOptions, "%s", input);
  aux->deltaT = atof(input);	/* L4 */
  fscanf(inputOptions, "%s", input);
  aux->us = atoi(input);	/* L5 */
  fscanf(inputOptions, "%s", input);
  aux->temI = atof(input);	/* L6 CHIDA*/
  fscanf(inputOptions, "%s", input);
  aux->temF = atof(input);	/* L7 CHIDA*/
  fscanf(inputOptions, "%s", input);
  aux->ts = atoi(input);	/* L8 */
  fscanf(inputOptions, "%s", input);
  aux->tiempoMax = atof(input);	/* L9 */
  fscanf(inputOptions, "%s", input);
  aux->device = atoi(input);	/* L10 */
  fscanf(inputOptions, "%s", input);
  aux->decimaciones = atoi(input); /* ;11 */


  aux->kas = (int)(aux->kMax/aux->deltaK);
  /* printf("kas:%d\n", aux->kas); */
  aux->sk = (double*)malloc(sizeof(double)*aux->kas);
  aux->gr = (double*)malloc(sizeof(double)*aux->kas);
  aux->Fcol = (double*)malloc(sizeof(double)*(TMMAXT)*aux->kas);
  aux->Fself = (double*)malloc(sizeof(double)*(TMMAXT)*aux->kas);
  //Creating the Directory
  sprintf(date,"%f-%f %d-%d %d:%d:%d",aux->temI, aux->temF, tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  mkdir(date, 0770);
  /* End of create Directory */
  /* Creating files in directory */
  /* directory : Directory  is the String where date going to be stored */
  sprintf(directory, "%s/efeself.dat",date); 
  efeself = fopen(directory, "w");
  /* printf("Que putas tiene directory = %s\n", directory); */
  /* sprintf(directory, "%s/datos.txt",date); */
  /* options = fopen(directory, "w"); */
  /* sprintf(ops, "\tdeltaK:%f\n\ttemI:%f\n\ttemF:%f\n\tus:%f\n\tts:%f\n\tphi:%f",aux->deltaK,aux->temI,aux->temF,aux->us,aux->ts,aux->phi); */
  /* fprintf(options, "%s" , ops); */
  fclose(efeself);
  /* printf("Tiempo maximo:%f\n",aux->tiempoMax); */
}






double blip(double temp){
  double eps, dr,r,dss;
  
  eps = 1/temp;
  dr = 1.0 / 50000;
  dss = 0;

  if(temp < 0.0000001){
    dss = 1;
  } else {
    for(r = dr; r <= 1; r+=dr){
      dss += dr * r*r * exp(-eps*(1.0/pow(r,12) - 2.0/pow(r,6) + 1));
    }
    dss = pow(1 - 3*dss,1/3.0);
  }

  return dss;
}


double* factorDeEstructura(scgle *data, char ruta_arch[]){
  double sk,skMax,skMin, i,k;
  double *arr = (double*)malloc(sizeof(double)*data->kas);
  int indice = 1;
  FILE *arch;
  char directory[100];
  /* printf("%s\n",ruta_arch); */
  skMax = skMin = sk = kCero(data->phi);
  printf("\n\ndeltaK %.17f\ndss =  %.17f\n",data->deltaK, data->dss );
  for(i = data->deltaK*data->dss; i <= data->kMax; i+=(data->deltaK*data->dss)){
    sk > skMax?skMax=sk:sk;
    sk < skMin?skMin=sk:sk;
  }
  data->skMax = skMax;
  data->skMin = skMin;
  arr[0] = sk = kCero(data->phi);
  printf("arr[0] = %.17f\n",arr[0] );
  for(i = data->deltaK*data->dss; indice < data->kas; i+=(data->deltaK*data->dss), indice++){
    arr[indice] = sK(data->phi,i);
    printf("arr[%d] = %.17f\n",indice,arr[indice] );
  }
  if(ruta_arch != NULL){
    sprintf(directory, "%s/%s", date,ruta_arch);
    arch = fopen(directory,"w");
    if(arch != NULL){
      for(i = 0, indice = 0; indice < data->kas; i+=(data->deltaK*data->dss), indice++){
	k = i * pow((data->phi - (data->phi*data->phi)/16.0)/data->phi,1/3.0);
	/* fprintf(arch,"%f %f %f %20.18f \n", indice*data->deltaK, k, k/data->dss, arr[indice]); */
	fprintf(arch,"%f %20.18f \n", k, arr[indice]);
      }
      fclose(arch);
    }
  }

  return arr;
}


double sK(double phi, double k){//Calculo de S(K)
  double sk,ck;

  k = k * pow((phi - (phi*phi)/16.0)/phi,1/3.0);
  phi = phi - (phi*phi)/16.0;
  ck = (-24.0 * phi / (pow(k,6) * pow(1-phi,4))) * (PCero(phi,k) + PUno(phi,k)*k*sin(k) + PDos(phi,k)*cos(k));
  sk = 1/(1-ck);
  return sk;
}


double kCero(double phi){
  double ck, sk;
        
  /* printf("#######Phi = %f\n",phi ); */
  phi = phi - (phi*phi)/16.0;
  ck = phi*(phi*phi*phi - 4*phi*phi + 2*phi -8)/pow(1-phi,4);
  sk = 1/(1.0-ck);

  return sk;
}



double PCero(double phi, double k){//Calculo de P0
  double p = 0;
  p = 3*phi*(pow(k,2)*pow(2+phi,2) + 4*pow(1+2*phi,2));
  return p;
}

double PUno(double phi, double k){//Calculo de P1
  double p = 0;
  p = -12*phi*pow(1+2*phi,2) +k*k*(1-6*phi + 5*pow(phi,3));
  return p;
}

double PDos(double phi, double k){//Calculo de P2
  double p = 0;
  p = -12*phi*pow(1 + 2*phi,2) + 3*phi*k*k*(-2+4*phi+7*phi*phi) - (pow(k,4)/2.0)*(2-3*phi+pow(phi,3));
  return p;
}

int  minmaxSk(scgle **data){
  int i,j;
  scgle *aux = *data;

  aux->skMax = aux->sk[0];
  for(i = 0; i < aux->kas; i++){
    if(aux->sk[i] > aux->skMax){
      aux->skMax = aux->sk[i];
      aux->skK = i;
      j = i;
    }
  }
  aux->skMin = aux->skMax;
  for(i = j; i < aux->kas; i++){
    if(aux->sk[i] < aux->skMin)
      aux->skMin = aux->sk[i];
  }

  return j;
}
