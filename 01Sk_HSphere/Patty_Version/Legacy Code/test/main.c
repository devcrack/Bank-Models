#include <stdio.h>
#include <string.h>
#include "scgle.h"

int main(){	
  char *ruta_arch, *input;
  scgle *data;
  FILE *inputOptions;
	
  /* Part of the constructor */
  input = (char*)malloc(sizeof(char)*100);
  inputOptions = fopen("options.txt","r");
  data = nuevaSCGLE();
  /* End of constructor */
  printf("Inician operaciones \n");
  if(allocate_some_Memory(data, inputOptions)) {
    factorDeEstructura(data, NULL);
  }
  /* if(inputOptions){ */
  /*   fscanf(inputOptions, "%s", input); */
  /*   data->kMax = atoi(input); */
  /*   fscanf(inputOptions, "%s", input); */
  /*   data->deltaK = atof(input); */
  /*   fscanf(inputOptions, "%s", input); */
  /*   data->phi = atof(input); */
  /*   fscanf(inputOptions, "%s", input); */
  /*   data->deltaT = atof(input); */
  /*   fscanf(inputOptions, "%s", input); */
  /*   ruta_arch = input; */
  /*   printf("Qsdfsdfdsue pedo\n"); */
  /* factorDeEstructura(data, NULL);  */
  /* printf("Kmax:%d \ndK:%f \nphi:%f \nkas:%d \n", data->kMax, data->deltaK, data->phi,data->kas); */
  /* printf("skMax: %f, %d \n", data->skMax, data->skMaxPos);	 */
  /* calculoTiemposCortos(data, "deltaZ.dat", "efeData.dat");	 */
  /* printf("Fin tiempos cortos, iniciando tiempos medios\n"); */
  /* calculoTiemposMedios(data,NULL,TCMAXT,(TMMAXT+TCMAXT),0); */
  /* printf("Inician decimaciones\n"); */
  /* decimaciones(data,30); */
  printf("Fin del programa\n");
return 0;
}

