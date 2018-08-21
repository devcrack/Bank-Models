#ifndef SCGLE_CUDA_H
#define SCGLE_CUDA_H

#define PI 3.14159265359
#include "cuda.h"
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

__global__ void sumTM(int nt, double *deltaZ, double *Fself, double *Fcol, double *sk, scgle data, double *sumDZ, double* arrFself, double* arrFcol, double* etaT){
        int i = threadIdx.x + 1;
        int threadCount = (int)(nt)/2 -1;

        arrFself[TMMAXT*blockIdx.x + threadIdx.x] = (deltaZ[nt-i] - deltaZ[nt-i-1]) * Fself[i*data.kas + blockIdx.x] + deltaZ[i] * (Fself[(nt-i) * data.kas + blockIdx.x] - Fself[(nt-i-1)*data.kas+blockIdx.x]);

	arrFcol[TMMAXT*blockIdx.x + threadIdx.x] = (deltaZ[nt-i] - deltaZ[nt-i-1]) * Fcol[i*data.kas + blockIdx.x] + deltaZ[i] * (Fcol[(nt-i) * data.kas + blockIdx.x] - Fcol[(nt-i-1)*data.kas+blockIdx.x]);

        __syncthreads();
        if(threadIdx.x == 0){
		int j; 
		double nk = data.deltaK * blockIdx.x, dzN = deltaZ[nt-1], lambda, Kcol, Kself;
		double sumFself = 0, sumFcol = 0;
		for(j = 0; j < threadCount; j++){
			sumFself += arrFself[TMMAXT*blockIdx.x + j];
			sumFcol += arrFcol[TMMAXT*blockIdx.x + j];
		}
		lambda = 1.0/(1 + powf(nk/(2*PI*data.alfa), data.mu));
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
		if(blockIdx.x < data.kas-1)
			etaT[blockIdx.x] = nk*nk*nk*nk*(sk[blockIdx.x+1]/sk[blockIdx.x] -1)*(sk[blockIdx.x+1]/sk[blockIdx.x] -1)* (Fcol[nt * data.kas + blockIdx.x]/sk[blockIdx.x]);
		
		sumDZ[blockIdx.x] = nk*nk*nk*nk * Fself[nt * data.kas + blockIdx.x] * Fcol[nt * data.kas + blockIdx.x] * ((sk[blockIdx.x]-1.0) / sk[blockIdx.x] ) * ((sk[blockIdx.x] - 1.0) / sk[blockIdx.x] );
        }
}

__global__ void sumCriterioDeAresto(scgle data, double *sk, double *gamadsol, double gamad){
	double nk, lambda;
	extern __shared__ double num[];

	nk = data.deltaK * threadIdx.x;
	lambda = 1.0/(1 + powf(nk/(2*PI*data.alfa), data.mu));
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
