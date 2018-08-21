#ifndef SCGLE_H
#define SCGLE_H

#include <stdio.h>
#include <stdlib.h>
#include<math.h>

#define TCMAXT 10 //estatica
#define TMMAXT 1000//estarica
#define PI 3.14159265359
/* #include "cuda.h" */
#define TCMAXT 10
#define TMMAXT 1000

typedef struct{
	int kMax;//estatica
	double deltaK;//estatica
	double phi;//estatica
	double deltaT;//estatica
	double skMax; //Dynamico
	int skMaxPos;//Dynamico
	double skMin;//Dynamico
	double mu;//estatica
	double alfa;//estatica
	double d;//estatica
	double b;//estatica
	int skK; //Dynamico
	int kas;//estatica
	double *sk;//Dynamico
	double *Fself;//Dynamico
	double *Fcol;//Dynamico
	double *deltaZ;//Dynamico
}scgle;



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
/* void calculoTiemposCortos(scgle *data, char* archZ, char* archF); */
/* void guardaDatosTiemposCortos(scgle* data, char* rutaArchZ, char* rutaArchF); */
/* void calculoTiemposMedios(scgle *data, char* archTiMed, double ntInicial, double ntFinal, double tiempoAcumlado); */
void decimaciones(scgle *data, int numDec);
int allocate_some_Memory(scgle *data, FILE * input_Ops);

#endif
