/*
 * scgle.h
 *
 *  Created on: 01/06/2016
 *  Author: Daniel Varela Varela
 *  @: Instituto de Fisica Manuel Sandoval Vallarta UASLP
 */
#ifndef SCGLE_H
#define SCGLE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "scgleCUDA.h" 		/*!!!!!!!!!!!!!!!!!!!!! */
#include "termalProps.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#define TCMAXT 10
#define TMMAXT 1000

char date[50];
FILE *efeself;
FILE *skmax;
FILE *grmax;
double sumTiempo = 0;

scgle *nuevaSCGLE();
double* loadSKFile(char *nombreArchivo, int kas, int kMax);
double PCero(double phi, double k);
double PUno(double phi, double k);
double PDos(double phi, double k);
double sK(double phi, double k);
double* factorDeEstructura(scgle *data, char* ruta_arch);
void calculoTiemposCortos(scgle *data, char* archZ, char* archF);
void guardaDatosTiemposCortos(scgle* data, char* rutaArchZ, char* rutaArchF);
void calculoTiemposMedios(scgle *data, char* archTiMed, double ntInicial, double ntFinal);
double sumaF(scgle *data, int nt, int nk, int i);
double sumaFself(scgle *data, int nt, int nk, int i);
double fcolTM(int nt, int nk, double sumFcol, scgle *data, double sk, double dzN);
double fseflTM(int nt, int nk, double sumFself, scgle *data, double dzN);
void decimaciones(scgle *data, int numDec);
void dinamico(scgle *data);
int minmaxSk(scgle **data);
double blip(double temp);
double uDeArresto(scgle *data);
void encuentraFself();
double quench(scgle *data, double uIni, double uFin, int num);


/* Parte fundamental del constructor */
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

void allocateMemorry(scgle **data, FILE *inputOptions){
  scgle *aux = *data;
  char *input,directory[100], ops[500];
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  FILE *options;

  input = (char*)malloc(sizeof(char)*200);
  fscanf(inputOptions, "%s", input);
  aux->kMax = atof(input);	/* L1 */
  printf("Kmax:%f\n", aux->kMax);
  fscanf(inputOptions, "%s", input);
  aux->deltaK = atof(input);	/* l2 */
  printf("deltaK:%f\n", aux->deltaK);
  fscanf(inputOptions, "%s", input);
  aux->phi = atof(input);	/* L3 CHIDA*/
  printf("phi:%f\n", aux->phi);
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
  aux->decimaciones = atoi(input); /* L11 */


  aux->kas = (int)(aux->kMax/aux->deltaK);
  printf("kas:%d\n", aux->kas);
  aux->sk = (double*)malloc(sizeof(double)*aux->kas);
  aux->gr = (double*)malloc(sizeof(double)*aux->kas);
  aux->Fcol = (double*)malloc(sizeof(double)*(TMMAXT)*aux->kas);
  aux->Fself = (double*)malloc(sizeof(double)*(TMMAXT)*aux->kas);
  sprintf(date,"%f-%f %d-%d %d:%d:%d",aux->temI, aux->temF, tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  mkdir(date, 0770);
  sprintf(directory, "%s/efeself.dat",date);
  efeself = fopen(directory, "w");
  sprintf(directory, "%s/datos.txt",date);
  options = fopen(directory, "w");
  sprintf(ops, "\tdeltaK:%f\n\ttemI:%f\n\ttemF:%f\n\tus:%f\n\tts:%f\n\tphi:%f",aux->deltaK,aux->temI,aux->temF,aux->us,aux->ts,aux->phi);
  fprintf(options, "%s" , ops);
  fclose(options);
  printf("Tiempo maximo:%f\n",aux->tiempoMax);
}

double kCero(double phi){
  double ck, sk;
        
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

double sK(double phi, double k){//Calculo de S(K)
  double sk,ck;

  k = k * pow((phi - (phi*phi)/16.0)/phi,1/3.0);
  phi = phi - (phi*phi)/16.0;
  ck = (-24.0 * phi / (pow(k,6) * pow(1-phi,4))) * (PCero(phi,k) + PUno(phi,k)*k*sin(k) + PDos(phi,k)*cos(k));
  sk = 1/(1-ck);
  return sk;
}

double* factorDeEstructura(scgle *data, char ruta_arch[]){
  double sk,skMax,skMin, i,k;
  double *arr = (double*)malloc(sizeof(double)*data->kas);
  int indice = 1;
  FILE *arch;
  char directory[100];

  skMax = skMin = sk = kCero(data->phi);
  for(i = data->deltaK*data->dss; i <= data->kMax; i+=(data->deltaK*data->dss)){
    sk > skMax?skMax=sk:sk;
    sk < skMin?skMin=sk:sk;
  }
  data->skMax = skMax;
  data->skMin = skMin;
  arr[0] = sk = kCero(data->phi);
  for(i = data->deltaK*data->dss; indice < data->kas; i+=(data->deltaK*data->dss), indice++){
    arr[indice] = sK(data->phi,i);
  }
  if(ruta_arch != NULL){
    sprintf(directory, "%s/%s", date,ruta_arch);
    arch = fopen(directory,"w");
    if(arch != NULL){
      for(i = 0, indice = 0; indice < data->kas; i+=(data->deltaK*data->dss), indice++){
	k = i * pow((data->phi - (data->phi*data->phi)/16.0)/data->phi,1/3.0);
	fprintf(arch,"%f %f %f %20.18f \n", indice*data->deltaK, k, k/data->dss, arr[indice]);
      }
      fclose(arch);
    }
  }

  return arr;
}

void calculoTiemposCortos(scgle *data, char* archZ, char* archF){
  int  ti, nk;
  double suma,k,i;
  double *Fself, *Fcol;

  Fself = data->Fself;
  Fcol = data->Fcol;

  for(i = data->deltaT, ti = 0; ti < TCMAXT; i+=data->deltaT, ti++){
    suma = 0;
    for(nk = 0, k = 0; nk < data->kas; nk++, k += data->deltaK){
      Fself[nk + (ti * data->kas)] = exp(-data->d*k*k*i);
      Fcol[nk + (ti * data->kas)] = exp((-data->d*k*k*i)/data->sk[nk]) * data->sk[nk];
      suma += k*k*k*k * Fself[nk + (ti * data->kas)] * Fcol[nk + (ti * data->kas)] * ((data->sk[nk]-1.0) / data->sk[nk] ) * ((data->sk[nk]-1.0) / data->sk[nk] );
    }
    data->deltaZ[ti] = (suma * data->deltaK * data->d) /(36*PI*(data->phi));
  }
}

void calculoTiemposMedios(scgle *data, char* archTiMed, double ntInicial, double ntFinal){
  int nt,i,j;//Host Memory
  double  *fself, *fcol, *sumDZ, sum, error, dzaprox, dzN, ERROR;//Host Memory
  double *d_sumDZ,*d_fself, *d_fcol, *d_deltaZ, *d_sk, *d_arrFself, *d_arrFcol; //Device Memory

  ERROR = 0.0001;
  fcol = data->Fcol;
  fself = data->Fself;
  sumDZ = (double*)malloc(sizeof(double)*data->kas);
  cudaMalloc((void**)&d_deltaZ, sizeof(double)*(TMMAXT));
  cudaMalloc((void**)&d_sk, sizeof(double)*data->kas);
  cudaMalloc((void**)&d_fself, sizeof(double)*(TMMAXT)*data->kas);
  cudaMalloc((void**)&d_fcol, sizeof(double)*(TMMAXT)*data->kas);
  cudaMalloc((void**)&d_sumDZ, sizeof(double)*(TMMAXT)*data->kas);
  cudaMemcpy(d_fself, fself, sizeof(double)*(TMMAXT)*data->kas, cudaMemcpyHostToDevice);
  cudaMemcpy(d_fcol, fcol, sizeof(double)*(TMMAXT)*data->kas, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sk, data->sk, sizeof(double)*data->kas, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_arrFself, sizeof(double)*(TMMAXT)*data->kas);
  cudaMalloc((void**)&d_arrFcol, sizeof(double)*(TMMAXT)*data->kas);
  cudaMemcpy(d_deltaZ, data->deltaZ, sizeof(double)*(TMMAXT), cudaMemcpyHostToDevice);
  for(nt = ntInicial; nt < ntFinal; nt++){
    i = (int)(nt/2 - 1);
    dzN = data->deltaZ[nt - 1];
    do{
      sum = 0;
      sumTM<<<data->kas,i>>>(nt, d_deltaZ, d_fself, d_fcol, d_sk, *data, d_sumDZ, d_arrFself, d_arrFcol);
      cudaMemcpy(sumDZ, d_sumDZ, sizeof(double)*data->kas, cudaMemcpyDeviceToHost);
      for(j = 0; j < data->kas; j++){
	sum += sumDZ[j];
      }      
      dzaprox = (data->d * sum * data->deltaK) / (36*(data->phi)*PI);
      error = fabs( ( dzaprox - dzN ) )/dzaprox;
      dzN = dzaprox;
    }while(error > ERROR);
    data->deltaZ[nt] = dzN;
    cudaMemcpy(d_deltaZ+nt, &data->deltaZ[nt], sizeof(double), cudaMemcpyHostToDevice);
  }
  cudaMemcpy(fself,d_fself, sizeof(double)*(TMMAXT)*data->kas, cudaMemcpyDeviceToHost);
  cudaMemcpy(fcol, d_fcol, sizeof(double)*(TMMAXT)*data->kas, cudaMemcpyDeviceToHost);
  cudaFree(d_deltaZ);
  cudaFree(d_sk);
  cudaFree(d_fcol);
  cudaFree(d_fself);
  cudaFree(d_sumDZ);
  cudaFree(d_arrFself);
  cudaFree(d_arrFcol);
  free(sumDZ);
}

void decimaciones(scgle *data, int numDec){
  int i, nt, j, k;
  FILE *fself;
  double t = 0.0000001;
  char nom[50];
  //arch = fopen("deltaZ01.dat","w");
  sprintf(nom, "%s/fself.dat",date);
  fself = fopen(nom,"w");

  for(j = 0; j < (TMMAXT); j++, t += data->deltaT){
    data->b += data->deltaZ[j] * data->deltaT;
    //fprintf(arch, "%20.18f\n", *(data->deltaZ + j));
    fprintf(fself, "%20.18f %20.18f\n", t, data->Fself[35 + j * data->kas]);
  }

  for(i = 0; i < numDec; i++){
    data->deltaT = data->deltaT*2;
    for(nt = 0; nt < (int)(TMMAXT)/2; nt++){
      for(k = 0; k < data->kas; k++){
	data->Fcol[k + (nt*data->kas)] = data->Fcol[k + (nt*2+1) * data->kas];
	data->Fself[k + (nt*data->kas)] = data->Fself[k + (nt*2+1) * data->kas];
      }
      data->deltaZ[nt] = data->deltaZ[nt*2+1];
    }
    calculoTiemposMedios(data, NULL, (int)(TMMAXT)/2, TMMAXT);
    for(j = (int)TMMAXT/2; j < TMMAXT; j++, t = t + data->deltaT){
      data->b += data->deltaZ[j] * data->deltaT;
      //	        fprintf(arch, "%20.18f\n", *(data->deltaZ + j));
      fprintf(fself, "%20.18f %20.18f\n", t, data->Fself[35 + j * data->kas]);
    }
  }
  data->b = 1/(1 + data->b);
  printf("Coeficiente de difucion: %20.18f\n",data->b);
  //   fclose(arch);
  fclose(fself);
  encuentraFself();
}

void dinamico(scgle *data){
  data->b = 0;
  data->deltaT = 0.0000001;
  calculoTiemposCortos(data, "deltaZ.dat", "efeData.dat");
  calculoTiemposMedios(data,NULL,TCMAXT,TMMAXT);
  decimaciones(data,30);
}

bool criterioDeArresto(scgle *data){
  bool band = false;
  double *d_sk; //variables del dispositivo
  double gamadsol,error,gamad,*d_gamadsol; //variables del host
  int i = 0;

  gamad = 0.000001;
  cudaMalloc((void**)&d_sk, sizeof(double)*data->kas);
  cudaMalloc((void**)&d_gamadsol, sizeof(double));
  cudaMemcpy(d_sk, data->sk, sizeof(double)*data->kas, cudaMemcpyHostToDevice);
  do{
    sumCriterioDeAresto<<<1, data->kas, sizeof(double)*data->kas>>>(*data, d_sk, d_gamadsol, gamad);
    cudaMemcpy(&gamadsol, d_gamadsol, sizeof(double), cudaMemcpyDeviceToHost);
    error = abs((gamadsol- gamad)/gamadsol);
    gamad = gamadsol;
    if(gamadsol > 1000){
      break;
    }
    if(error < 0.0001){
      band = true;
      break;
    }
    i++;
  }while(i < 1000);

  cudaFree(d_gamadsol);
  cudaFree(d_sk);

  return band;
}

double uDeArresto(scgle *data){
  double u,uAux,adk, alfadk, deltaU;
  int k;
  uAux = 0;
  deltaU = 0.1;

  do{
    for(u = uAux; u <= data->uMax; u+=deltaU){
      for(k = 0; k < data->kas; k++){
	adk =  (k * data->deltaK)/data->dss;
	alfadk = (2*adk*adk)/data->skf[k];
	data->sk[k] = data->ski[k]*exp(-u*alfadk) + data->skf[k]*(1-exp(-u*alfadk));
      }
      data->arrestado = criterioDeArresto(data);
      if(data->arrestado){
	break;
      }
    }
    uAux = u - deltaU;
    deltaU = deltaU/10.0;
  }while(deltaU >= 0.0000001);
  printf("U de arresto: %f\n", u);
  return u;
}

double quench(scgle *data, double uIni, double uFin, int num){
  double u,adk, alfadk,t, du, grMax;
  int i,k,aux;
  char bdtfile[100];
  sprintf(bdtfile, "%s/bdt%d.dat",date,num);
  FILE *bdt = fopen(bdtfile, "w");
  t = 0;
  printf("UFinal[%d]:%f\n", num, uFin);
  t = 0;
  for(u = uIni, i = 0; (u < uFin && t < data->tiempoMax) || data->arrestado; u+=data->deltaU, i++){
    for(k = 0, adk = 0; k < data->kas; k++, adk+=data->deltaK){
      alfadk = (2*adk*adk)/data->skf[k];
      data->sk[k] = data->ski[k]*exp(-u*alfadk) + data->skf[k]*(1-exp(-u*alfadk));
    }
    grMax = fft(data, data->sk);
    aux = minmaxSk(&data);
    dinamico(data);
    fprintf(bdt, "%20.18f\t %20.18f\n", data->b, t+sumTiempo);
    fprintf(skmax, "%20.18f\t %20.18f\t %20.18f\t %d\n", data->u + u, t+sumTiempo,  data->skMax, aux);
    fprintf(grmax, "%20.18f\t %20.18f\t %20.18f\n", data->u + u, t+sumTiempo,  grMax);
    t += data->deltaU/data->b;
  }

  du = ((data->tiempoMax - t)/10) * data->b;
  while(t < data->tiempoMax){
    fprintf(skmax, "%20.18f\t %20.18f\t %20.18f\n", data->u + u, t+sumTiempo,  data->skMax);
    fprintf(grmax, "%20.18f\t %20.18f\t %20.18f\n", data->u + u, t+sumTiempo,  grMax);
    u += du;
    t += du/data->b;
  }
  sumTiempo = (num+1) * data->tiempoMax;
  data->u += u;
  fclose(bdt);

  return t;
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

void encuentraFself(){
  char t[50], f[50], nom[50];
  double n;
  FILE* arch;

  sprintf(nom, "%s/fself.dat",date);
  arch = fopen(nom,"w");
  if(arch){
    while(!feof(arch)) {
      fscanf(arch, "%s", t);
      fscanf(arch, "%s", f);
      n = atof(f);
      if(n <= 0.36787944117){
	fprintf(efeself, "%s\t %s\n",t,f);
	break;
      }
    }
    fclose(arch);
  } else {
    printf("Imposible abrir archivo\n");
  }
}

void guardaSK(scgle *data, char dir[], int n){
  FILE *ski, *skf;
  char nom[50];
  int i;
  float k;

  sprintf(nom, "%s/ski%d.dat", dir, n);
  ski = fopen(nom, "w");
  sprintf(nom, "%s/skf%d.dat", dir, n);
  skf = fopen(nom, "w");

  for(i = 0, k = 0; i < data->kas; i++, k+=data->deltaK){
    fprintf(ski, "%f\t%f\n", k, data->ski[i]);
    fprintf(skf, "%f\t%f\n", k, data->skf[i]);
  }
  fclose(ski);
  fclose(skf);
}

void guardaGR(scgle *data,char dir[], int n){
  FILE *gr;
  char nom[50];
  int i;
  float k;

  sprintf(nom, "%s/g(r)%d.dat", dir, n);
  gr = fopen(nom, "w");

  for(i = 0, k = 0; i < data->kas; i++, k+=data->deltaK){
    fprintf(gr, "%f\t%f\n", k, data->gr[i]);
  }
  fclose(gr);
}

#endif







