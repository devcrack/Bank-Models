/*
 *  scgleCUDA.h
 *
 *  Created on: 01/06/2016
 *  Author: Daniel Varela Varela
 *  @: Instituto de Fisica Manuel Sandoval Vallarta UASLP
 */
#ifndef SCGLE_CUDA_H
#define SCGLE_CUDA_H

#include "cuda.h"
#define PI 3.14159265359
#define TCMAXT 10
#define TMMAXT 1000

/* Clase : p */
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

__global__ void sumTM(int nt, double *deltaZ, double *Fself, double *Fcol, double *sk, scgle data, double *sumDZ, double* arrFself, double* arrFcol){
        int i = threadIdx.x + 1;
        int threadCount = (int)(nt)/2 -1;

        arrFself[TMMAXT*blockIdx.x + threadIdx.x] = (deltaZ[nt-i] - deltaZ[nt-i-1]) * Fself[i*data.kas + blockIdx.x] + deltaZ[i] * (Fself[(nt-i) * data.kas + blockIdx.x] - Fself[(nt-i-1)*data.kas+blockIdx.x]);

		arrFcol[TMMAXT*blockIdx.x + threadIdx.x] = (deltaZ[nt-i] - deltaZ[nt-i-1]) * Fcol[i*data.kas + blockIdx.x] + deltaZ[i] * (Fcol[(nt-i) * data.kas + blockIdx.x] - Fcol[(nt-i-1)*data.kas+blockIdx.x]);

        __syncthreads();
        if(threadIdx.x == 0){
			int j; 
			double sumFself = 0, sumFcol = 0, dzN = deltaZ[nt-1], lambda, Kcol, Kself, nk = blockIdx.x*data.deltaK;;
			for(j = 0; j < threadCount; j++){
				sumFself += arrFself[TMMAXT*blockIdx.x + j];
				sumFcol += arrFcol[TMMAXT*blockIdx.x + j];
			}
			lambda = nk/(2*PI*data.alfa);
		    lambda = 1.0/(1 + pow(lambda, data.mu));
		    Kcol = (1.0/data.deltaT) + (data.d*pow(nk,2))/sk[blockIdx.x] + lambda * deltaZ[0];
		    Kself = (1.0/data.deltaT) + data.d*pow(nk,2) + lambda * deltaZ[0];

		    Fself[nt * data.kas + blockIdx.x] =  (1.0/Kself)
		                    * ((lambda * dzN) * (1 - Fself[blockIdx.x])
		                    - lambda * deltaZ[(int)nt/2] * Fself[(int)((nt/2) * data.kas + blockIdx.x)]
		                    + lambda * deltaZ[nt - 1] * Fself[blockIdx.x]
		                    + (1.0/data.deltaT + lambda * deltaZ[0]) * Fself[(int)((nt-1) * data.kas + blockIdx.x)]
		                    - lambda * sumFself);

		    Fcol[nt * data.kas + blockIdx.x] =   (1.0/Kcol)
		                    * ((lambda * dzN) * (sk[blockIdx.x] - Fcol[blockIdx.x])
		                    - lambda * deltaZ[(int)nt/2] * Fcol[(int)((nt/2) * data.kas + blockIdx.x)]
		                    + lambda * deltaZ[nt - 1] * Fcol[blockIdx.x]
		                    + (1.0/data.deltaT + lambda * deltaZ[0]) * Fcol[(int)((nt-1) * data.kas + blockIdx.x)]
		                    - lambda * sumFcol);

			sumDZ[blockIdx.x] = nk*nk*nk*nk * Fself[nt * data.kas + blockIdx.x] * Fcol[nt * data.kas + blockIdx.x] * ((sk[blockIdx.x]-1.0) / sk[blockIdx.x] ) * ((sk[blockIdx.x] - 1.0) / sk[blockIdx.x] );
        }
}

__global__ void sumCriterioDeAresto(scgle data, double *sk, double *gamadsol, double gamad){
	double nk, lambda;
	extern __shared__ double num[];

	nk = threadIdx.x*data.deltaK;
	lambda = nk/(2*PI*data.alfa);
	lambda = 1.0/(1 + pow(lambda, data.mu));
	num[threadIdx.x] = (pow(sk[threadIdx.x]-1,2)*lambda*lambda)/
					   ((lambda*sk[threadIdx.x]+gamad*nk*nk)*(lambda+gamad*nk*nk));
	num[threadIdx.x] = num[threadIdx.x]*pow(nk,4);

	__syncthreads();
	if(threadIdx.x == 0){
		double sum = 0;
		for(int i = 0; i < data.kas; i++){
			sum += num[i];
		}
		sum = (sum*data.deltaK)/(36*data.phi*PI);
		*gamadsol = 1/sum;
	}
}

#endif
