#ifndef SCGLE_H
#define SCGLE_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<time.h>
#include<unistd.h>
#include<math.h>
#include<omp.h>
#include<stdbool.h>

#define PI 3.14159265359

#define TCMAXT 10
#define TMMAXT 1000

/* typedef int bool; */
/* #define true 1 */
/* #define false 0 */

typedef struct{
  double kMax;
  double deltaK;
  double phi;
  double deltaT;
  double skMax;
  double skMin;
  double mu;
  double alfa;
  double d;
  double b;
  double u;
  double us;
  double uMax;
  double dss;
  double deltaU;
  int skK;
  int kas;
  int device;
  double *sk;
  double *gr;
  double *ski;
  double *skf;
  double *Fself;
  double *Fcol;
  double *deltaZ;
  double temI;
  double temF;
  double deltaTemp;
  double ts;
  double tiempoMax;
  int decimaciones;
  bool arrestado;
}scgle;


char date[50];
FILE *efeself;
FILE *skmax;
FILE *grmax;
FILE *inputOptions;
double
sumTiempo,
  phi,
  deltaTemp,
  t,
  uArresto;
scgle *data;
char nom[50];

scgle *nuevaSCGLE();
void allocateMemorry(scgle **data, FILE *inputOptions);
double blip(double temp);
double* factorDeEstructura(scgle *data, char ruta_arch[]);
double sK(double phi, double k);
double kCero(double phi);
double PCero(double phi, double k);
double PUno(double phi, double k);
double PDos(double phi, double k);
int  minmaxSk(scgle **data);

#endif
