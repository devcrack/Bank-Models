ifndef SCGLE_H
#define SCGLE_H

#include <stdio.h>
#include <stdlib.h>
#include "scgleCUDA.h"

#define TCMAXT 10 //estatica
#define TMMAXT 1000//estarica

scgle *nuevaSCGLE();//constructor chafa
int  allocateMemorry(scgle **data, char *arch_sk);////constructor chafa

//HS 
double PCero(double phi, double k);
double PUno(double phi, double k);
double PDos(double phi, double k);
double sK(double phi, double k);
void factorDeEstructura(scgle *data, char* ruta_arch);
//HS

//Modulo dínamico
void leeFactorDeEstructura(scgle *data, char* ruta_arch);
void calculoTiemposCortos(scgle *data, char* archZ, char* archF);
void guardaDatosTiemposCortos(scgle* data, char* rutaArchZ, char* rutaArchF);
void calculoTiemposMedios(scgle *data, char* archTiMed, double ntInicial, double ntFinal, double tiempoAcumlado);
void decimaciones(scgle *data, int numDec);
//Modulo dínamico

scgle *nuevaSCGLE(){
  scgle *nva;

        
  nva = (scgle*)malloc(sizeof(scgle));
  nva->kMax = 0;
  nva->deltaK = 0;
  nva->phi = 0;
  nva->deltaT = 0;
  nva->skMax = 0;
  nva->skMin = 0;
  nva->mu = 2;
  nva->alfa = 1.305;
  nva->d = 1;
  nva->kas = 0;
	
  nva->deltaZ = (double*)malloc(sizeof(double)*(TMMAXT+TCMAXT));
        
  return nva;
}

int  allocateMemorry(scgle **data, char *arch_sk){
  scgle *aux = *data;
    	
  if(arch_sk != NULL) {
    leeFactorDeEstructura(aux, arch_sk);
    aux->Fcol = (double*)malloc(sizeof(double)*(TMMAXT+TCMAXT)*aux->kas);
    aux->Fself = (double*)malloc(sizeof(double)*(TMMAXT+TCMAXT)*aux->kas);
    return 1;
  } else {
    aux->kas = (int)(aux->kMax/aux->deltaK);
    aux->sk = (double*)malloc(sizeof(double)*aux->kas);
    aux->Fcol = (double*)malloc(sizeof(double)*(TMMAXT+TCMAXT)*aux->kas);
    aux->Fself = (double*)malloc(sizeof(double)*(TMMAXT+TCMAXT)*aux->kas);
  }	
  return 0;
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


/* 
 *From here towards pCero, PUno, PDos
 */
double sK(double phi, double k){//Calculo de S(K)
  double sk,ck;
        
  k = k * pow((phi - (phi*phi)/16.0)/phi,1/3.0);
  phi = phi - (phi*phi)/16.0;
  ck = (-24.0 * phi / (pow(k,6) * pow(1-phi,4))) * (PCero(phi,k) + PUno(phi,k)*k*sin(k) + PDos(phi,k)*cos(k));
  sk = 1/(1-ck);
  return sk;
}


/* 
 * From here towards sk
 */
void factorDeEstructura(scgle *data, char* ruta_arch){
  double sk,skMax,skMin, i;
  int indice = 1;
  FILE *arch;

  skMax = skMin = sk = kCero(data->phi);
  for(i = data->deltaK; i <= data->kMax; i+=data->deltaK){
    sk > skMax?skMax=sk:sk;
    sk < skMin?skMin=sk:sk;
  }
  data->skMax = skMax;
  data->skMin = skMin;
  skMax = data->sk[0] = sk = kCero(data->phi);
  for(i = data->deltaK; i <= data->kMax; i+=data->deltaK, indice++){
    data->sk[indice] = sK(data->phi,i);
    if(data->sk[indice] > skMax){
      data->skMaxPos = indice;
      skMax = data->sk[indice];
    }
  }
  if(ruta_arch != NULL){
    arch = fopen(ruta_arch,"w");
    if(arch != NULL){
      for(i = 0, indice = 0; i < data->kMax; i+= data->deltaK, indice++){
	fprintf(arch,"%f %20.18f \n",i,data->sk[indice]);
      }
      fclose(arch);
    }
  }
}

void leeFactorDeEstructura(scgle *data, char* ruta_arch){
  FILE *arch;
  double k1,k2,sk;
  int i = 0;
  printf("hola\n");
  arch = fopen(ruta_arch, "r");
  if(arch){
    fscanf(arch, "%f", &k1);
    fscanf(arch, "%f", &sk);
    fscanf(arch, "%f", &k2);
    data->deltaK = k2 - k1;
    fseek(arch, 0, SEEK_SET);
    printf("aaaa %f\n", data->deltaK);
    do{
      fscanf(arch, "%f", &k1);
      fscanf(arch, "%f", &sk);
    }while(fgetc(arch) != EOF);	
    data->kMax = k1;
    data->kas = k1/data->deltaK;
    data->sk = (double*)malloc(sizeof(double)*data->kas);
    data->Fcol = (double*)malloc(sizeof(double)*(TMMAXT+TCMAXT)*data->kas);
    data->Fself = (double*)malloc(sizeof(double)*(TMMAXT+TCMAXT)*data->kas);
    fseek(arch, 0, SEEK_SET);
    do{
      fscanf(arch, "%f", &k1);
      fscanf(arch, "%f", &sk);
      data->sk[i] = sk;
      i++;
      sk > data->skMax?data->skMax=sk:sk;
      sk < data->skMin?data->skMin=sk:sk;
    }while(fgetc(arch) != EOF);	
  }
}

void calculoTiemposCortos(scgle *data, char* archZ, char* archF){
  int  ti, nk;
  double suma,k,i;
  double *Fself, *Fcol;

  Fself = data->Fself;
  Fcol = data->Fcol;

  for(i = data->deltaT, ti = 0; ti < TCMAXT; i+=data->deltaT, ti++){
    suma = 0;
    for(k = 0, nk = 0; k < data->kMax; k += data->deltaK, nk++){
      Fself[nk + (ti * data->kas)] = exp(-data->d*k*k*i);
      Fcol[nk + (ti * data->kas)] = exp((-data->d*k*k*i)/data->sk[nk]) * data->sk[nk];
      suma += k*k*k*k * Fself[nk + (ti * data->kas)] * Fcol[nk + (ti * data->kas)] * ((data->sk[nk]-1.0) / data->sk[nk] ) * ((data->sk[nk]-1.0) / data->sk[nk] );
    }
    data->deltaZ[ti] = (suma * data->deltaK * data->d) /(36*PI*(data->phi));
  }
}

void calculoTiemposMedios(scgle *data, char* archTiMed, double ntInicial, double ntFinal, double tiempoAcumlado){
  int nt,i,j;//Host Memorry
  double  *fself, *fcol, *sumDZ, sum, dzaprox, dzN, *etaT, sumEta;//Host Memorry
  double *d_sumDZ,*d_fself, *d_fcol, *d_deltaZ, *d_sk, *d_etaT,*d_arrFself, *d_arrFcol; //Device Memory
  double error,ERROR;

  cudaSetDevice(0);
  ERROR = 0.0001;
  fcol = data->Fcol;
  fself = data->Fself;
  sumDZ = (double*)malloc(sizeof(double)*data->kas);
  etaT = (double*)malloc(sizeof(double)*data->kas);
  cudaMalloc((void**)&d_deltaZ, sizeof(double)*(TMMAXT+TCMAXT));
  cudaMalloc((void**)&d_sk, sizeof(double)*data->kas);
  cudaMalloc((void**)&d_fself, sizeof(double)*(TMMAXT+TCMAXT)*data->kas);
  cudaMalloc((void**)&d_fcol, sizeof(double)*(TMMAXT+TCMAXT)*data->kas);
  cudaMalloc((void**)&d_sumDZ, sizeof(double)*data->kas);
  cudaMalloc((void**)&d_etaT, sizeof(double)*data->kas);
  cudaMemcpy(d_fself, fself, sizeof(double)*(TMMAXT+TCMAXT)*data->kas, cudaMemcpyHostToDevice);
  cudaMemcpy(d_fcol, fcol, sizeof(double)*(TMMAXT+TCMAXT)*data->kas, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sk, data->sk, sizeof(double)*data->kas, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_arrFself, sizeof(double)*(TMMAXT+TCMAXT)*data->kas);
  cudaMalloc((void**)&d_arrFcol, sizeof(double)*(TMMAXT+TCMAXT)*data->kas);
  cudaMemcpy(d_deltaZ, data->deltaZ, sizeof(double)*(TMMAXT+TCMAXT), cudaMemcpyHostToDevice);

  for(nt = ntInicial; nt < ntFinal; nt++){
    i = (int)(nt/2 - 1);
    dzN = data->deltaZ[nt - 1];
    do{
      sum = 0;
      sumEta = 0;
      sumTM<<<data->kas,i>>>(nt, d_deltaZ, d_fself, d_fcol, d_sk, *data, d_sumDZ, d_arrFself, d_arrFcol, d_etaT);
      cudaMemcpy(sumDZ, d_sumDZ, sizeof(double)*data->kas, cudaMemcpyDeviceToHost);    
      cudaMemcpy(etaT, d_etaT, sizeof(double)*data->kas, cudaMemcpyDeviceToHost);    
      for(j = 0; j < data->kas; j++){
	sum += sumDZ[j];
	sumEta += etaT[j];
      }      
      dzaprox = (data->d * sum * data->deltaK) / (36*(data->phi)*PI);
      error = fabs( ( dzaprox - dzN ) )/dzaprox;
      dzN = dzaprox;
    }while(error > ERROR);
    data->deltaZ[nt] = dzN;
    cudaMemcpy(d_deltaZ+nt, &data->deltaZ[nt], sizeof(double), cudaMemcpyHostToDevice);
  }
  cudaMemcpy(fself,d_fself, sizeof(double)*(TMMAXT+TCMAXT)*data->kas, cudaMemcpyDeviceToHost);
  cudaMemcpy(fcol, d_fcol, sizeof(double)*(TMMAXT+TCMAXT)*data->kas, cudaMemcpyDeviceToHost);
  cudaFree(d_deltaZ);
  cudaFree(d_sk);
  cudaFree(d_fcol);
  cudaFree(d_fself);
  cudaFree(d_sumDZ);
  cudaFree(d_arrFself);
  cudaFree(d_arrFcol);
  cudaFree(d_etaT);
  free(sumDZ);
}

void decimaciones(scgle *data, int numDec){
  int i, nt, j, k, x;
  FILE *fself,*arch,*coef;
  double t = 0.0000001;
  arch = fopen("deltaZ01.dat","w");
  fself = fopen("Efes.dat","w");
  coef = fopen("Coeficiente.dat","w");
	
  x = (int)(2.7182818/data->deltaK);
  for(j = 0; j < (TMMAXT); j++, t += data->deltaT){
    data->b += data->deltaZ[j] * data->deltaT;
    fprintf(arch, "%20.18f\n", *(data->deltaZ + j));
    fprintf(fself, "%20.18f %20.18f %20.18f\n", t, data->Fcol[x + j * data->kas],data->Fself[x + j * data->kas]);
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
    calculoTiemposMedios(data, NULL, (int)(TMMAXT)/2, TMMAXT, t);
    for(j = (int)TMMAXT/2; j < TMMAXT; j++, t = t + data->deltaT){
      data->b += data->deltaZ[j] * data->deltaT;
      fprintf(arch, "%20.18f\n", data->deltaZ[j]);
      fprintf(fself, "%20.18f %20.18f %20.18f\n", t, data->Fcol[x + j * data->kas],data->Fself[x + j * data->kas]);
    }
  }
  data->b = 1/(1 + data->b);
  fprintf(coef, "%20.18f\n", data->b);
  fclose(arch);
  fclose(fself);
  fclose(coef);
}

#endif
