/*
 *	fft.h
 *
 *	Created on: 21/03/2017
 *	Author: est1
 */

#ifndef FFT_H
#define FFT_H

#include "scgleCUDA.h"

double fft(scgle *data);

double fft(scgle *data, double *sk) {
	int i,j;
	double vr, rmax, dr, sumak, rho, vk, gMax;
	double *hk,*hr;

	hk = new double[data->kas];
	hr = new double[data->kas];

	gMax = 0;
	rmax = 10;
	dr = rmax/(data->kas*1.0);
	rho = 6.0*data->phi/PI;

	for(i = 0; i < data->kas; i++)
		hk[i] = sk[i] - 1;

	for(i = 1, vr = 0; i < data->kas; i++, vr+=dr) {
		sumak = 0;
		for(j = 0; j < data->kas-1; j++) {
			vk = (j-1)*data->deltaK;
         	sumak += data->deltaK*vk*hk[j]*sin(vk*vr);
		}
		hr[i] = sumak/(2*vr*rho*PI*PI);
		data->gr[i] = hr[i] + 1;
		if(data->gr[i] > gMax)
			gMax = data->gr[i];
	}
	delete(hk);
	delete(hr);

	return gMax;
}

#endif
