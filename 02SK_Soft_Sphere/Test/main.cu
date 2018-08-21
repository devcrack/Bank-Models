/*
 * main.cu
 *
 * Created on: 01/06/2016
 * Author: Daniel Varela Varela
 * @: Instituto de Fisica Manuel Sandoval Vallarta UASLP
 * nvcc -arch=sm_35 -rdc=true -lcudadevrt main.cu -o aging
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "scgle.h"
#include <omp.h>

int main(){
  scgle *data;
  FILE *inputOptions;
  double phi,deltaTemp,t, uArresto;
  inputOptions = fopen("options.txt","r");
  data = nuevaSCGLE();
  int i;
  char nom[50];

  if(inputOptions){
    allocateMemorry(&data,inputOptions);
    cudaSetDevice(data->device); /* ._. */
    t = data->temI;
    phi = data->phi;
    deltaTemp = (data->temI-data->temF)/data->ts;
    data->dss = blip(t);
    data->phi = phi*pow(data->dss,3);
    printf("Dss: %f, phi:%f\n",data->dss,data->phi);
    data->ski = factorDeEstructura(data,"factor.dat"); /* EL CHIDO */
    sprintf(nom,"%s/skmax%d.dat",date,0);
    skmax = fopen(nom, "w");
    sprintf(nom,"%s/grmax%d.dat",date,0);
    grmax = fopen(nom, "w");
    t-=deltaTemp;
    data->dss = blip(t);//****
    data->phi = phi*pow(data->dss,3);//***
    printf("Dss: %f, phi:%f\n",data->dss,data->phi);
    sprintf(nom,"factor%d.dat",0);
    data->skf = factorDeEstructura(data,nom);
    data->sk = (double*)malloc(sizeof(double)*data->kas);
    data->phi = phi;
    minmaxSk(&data);		/*  */
    uArresto = uDeArresto(data);
    data->deltaU = uArresto/data->us;
    quench(data, 0, 2,0);
    memcpy(data->ski, data->sk, sizeof(double)*data->kas);
    free(data->skf);
    fclose(skmax);
    fclose(grmax);
    fft(data, data->sk);
    for(i = 1; i < data->ts && !data->arrestado; i++, t -= deltaTemp){
      sprintf(nom,"%s/skmax%d.dat",date,i);
      skmax = fopen(nom, "w");
      sprintf(nom,"%s/grmax%d.dat",date,0);
      grmax = fopen(nom, "w");
      data->dss = blip(t-deltaTemp);//*****
      data->phi = phi*pow(data->dss,3);//*****
      printf("Dss: %f, phi:%f\n",data->dss,data->phi);
      sprintf(nom,"factor%d.dat",i);
      data->skf = factorDeEstructura(data,nom);
      data->phi = phi;
      uArresto = uDeArresto(data);
      data->deltaU = uArresto/data->us;
      quench(data, 0, uArresto,i);
      free(data->ski);
      data->ski = data->sk;
      data->sk = (double*)malloc(sizeof(double)*data->kas);
      free(data->skf);
      fclose(skmax);
      fclose(grmax);
    }
    fclose(efeself);
    printf("Fin del Programa\n");
  }

  return 0;
}

