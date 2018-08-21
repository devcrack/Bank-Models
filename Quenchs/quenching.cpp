/**
* Enfriamientos 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dk 0.05
#define MAX_K 15/dk


bool asignate_memory(double **arr, int n);
bool asignate_memory_2D(double ***arr_bidi, int nRen, int nCol);
void free_memory_2D(double **arr, int nRen);
void blip(double a_temp, double *a_dss);
void calculate_sk(double *sk, double a_phi);
bool criterio_arresto(double *sk, double a_phi);
double calculate_u_arresto(double a_dss, double a_phi, double *ski, double *skf);
void calcula(double a_temp_initial, double a_temp_final, double a_phi, double *u_arresto);
void save_sks_in_file(double **arr, int nRen, int nCol);
  
int main(int argc, const char *argv[]){
  double u, tmp_1, tmp_2, phi;
    if(argc == 4) {
      tmp_1 = atof(*(argv + 1));
      tmp_2 = atof(*(argv + 2));
      phi = atof(*(argv + 3));
      
    }
    else { 
      tmp_1 = 6.5431891;
      tmp_2 = 0.0000045;
      phi = 0.537;
    }
    calcula(tmp_1, tmp_2, phi, &u);
    return 1;
}

bool asignate_memory(double **arr, int n){
    *arr = (double *) malloc(sizeof(double) * n);
    if(*arr)
        return true;
    else
        return false;
}

/*Función Asigna Memoría Arreglo Bidimensional Dinámico*/
bool asignate_memory_2D(double ***arr_bidi, int nRen, int nCol)
{
    bool res = false;
    int i, k;

    *arr_bidi = (double **)malloc(sizeof(double *) * nRen);
    if(*arr_bidi){
         res = true;
         for(k = 0; k < nRen && res ;k++){
              *(*arr_bidi + k) = (double *)malloc(sizeof(double) * nCol);
              if(!*(*arr_bidi + k)){
                  while(--k >= 0)
                      free(*(*arr_bidi + k));
                  free(*arr_bidi);
                  res = false;
              }
         }
    }
}

/*Liberar Memoria Arreglo Bidimensional Dinámico*/
void free_memory_2D(double **arr, int nRen){
    for(int i = 0; i < nRen; i++)
        free(*(arr + i));  
    free(arr); 
}

void blip(double a_temp, double *a_dss){
    double eps, dr;

    eps = 1 / a_temp;
    dr = 1.00 / 50000;
    (*a_dss) = 0;
    if(a_temp <= 0.0000001){
        (*a_dss) = 1;
    }else{
        for(double r = dr; r <= 1 ; r += dr){
            (*a_dss) = (*a_dss) + dr * r * r * exp(-eps * (1.0/pow(r, 12) - 2.0 / pow(r, 6) + 1));  
        }
        (*a_dss) = pow(1 - 3 * (*a_dss), 1/3.0);
    }
}

void calculate_sk(double *sk, double a_phi){
    double sk_aux, ck, p0, p1, p2, k;
    
    sk[0] = 1/(1.0 - a_phi * (a_phi * a_phi * a_phi - 4 * a_phi * a_phi + 2*a_phi -8) /pow(1 - a_phi, 4));
    printf("Sk[0] = %lf", sk[0]);
    for(int i = 1; i < MAX_K ; i++){
        k = i *dk;
        p0 = 3 * a_phi * ( pow(k, 2) * pow(2 + a_phi, 2) + 4 *pow(1+2*a_phi,2));
	// printf("p0 = %lf\n", p0);
        p1 = -12*a_phi*pow(1+2*a_phi,2)+k*k*(1-6*a_phi + 5*pow(a_phi,3));
	// printf("p1 = %lf\n", p1);
        p2 = -12 *a_phi*pow(1 + 2*a_phi,2) + 3*a_phi*k*k*(-2+4*a_phi+7*a_phi*a_phi) - (pow(k,4)/2.0)*(2-3*a_phi+pow(a_phi,3));
	// printf("p2 = %lf\n", p2);
        ck = (-24.0 * a_phi / (pow(k,6) * pow(1-a_phi,4))) * (p0 + p1*k*sin(k) + p2*cos(k));
	// printf("ck = %lf\n", ck);
        sk_aux = 1/(1.0 - ck);
        sk[i] = sk_aux;
	printf("sk[%d]  = %lf\n", i, sk[i]);
    }
}

bool criterio_arresto(double *sk, double a_phi){
    double gamma_sol, error, gamma, lambda, sum;
    bool band = false;
    int i = 0;
    int j = 0;
    
    gamma = 0.000001;
    do{
        for(float nk = 0; j < MAX_K; nk += dk, j++){
            lambda = nk / (2 * M_PI * 1.0305);
            lambda = 1.0 / (1 + pow(lambda, 2));
            sum += ((pow(sk[j]-1, 2) * lambda * lambda) / ((lambda * sk[j] + gamma * nk * nk) * (lambda + gamma * nk * nk))) * pow(nk, 4);
        }
        sum = (sum * dk) / (36 * a_phi * M_PI);
        gamma_sol = 1 / sum;
        error = abs((gamma_sol - gamma) / gamma_sol);
        gamma = gamma_sol;
        if(gamma_sol > 1000)
            break;
        if(error < 0.0001){
            band = true;
            break;
        }
        i++;
    }while(i < 1000);
    return band;
}

double calculate_u_arresto(double a_dss, double a_phi, double *ski, double *skf){
    double u, u_aux, adk, alfadk, delta_u, k;
    double *sk;

    u_aux = 0;
    delta_u = 0.01;
    do{
        for(u = u_aux; u <= 0.2; u += delta_u){
            if(asignate_memory(&sk, MAX_K)){
                for(int k = 0 ;k < MAX_K; k++){
                    adk = (k * dk) / a_dss;
                    alfadk = (2 * adk * adk) / skf[k];
                    sk[k] = ski[k] * exp(-u * alfadk) + skf[k] * (1 - exp(-u * alfadk));
                }
                if(criterio_arresto(sk, a_phi)){
                    free(sk);
                    break;
                }                    
                free(sk);
            }
        }
        u_aux = u - delta_u;
		delta_u = delta_u / 10.0;
    }while(delta_u >= 0.0000001);   
    return u;
}

void calcula(double a_temp_initial, double a_temp_final, double a_phi, double *u_arresto){
    double alfa_dk, delta_u, a_dss;
    double phi_initial, phi_final;
    double *sk_initial;
    double *sk_final;
    double **sks;
    int i;
    double u, adk;

    if(asignate_memory(&sk_initial, MAX_K) && asignate_memory(&sk_final, MAX_K) && asignate_memory_2D(&sks, 11, MAX_K)){
        blip(a_temp_initial, &a_dss);
        phi_initial = a_phi * pow(a_dss, 3);
        calculate_sk(sk_initial, phi_initial);    
        blip(a_temp_final, &a_dss);
        phi_final = a_phi * pow(a_dss, 3);
        calculate_sk(sk_final, phi_final);
        printf("Valor de phi = %lf\n", phi_final);
        if (phi_final > 0.582) {
            printf("ARRESTADO\n");
            (*u_arresto) = fabs(calculate_u_arresto(a_dss, a_phi, sk_initial, sk_final));
            delta_u = (*u_arresto) / 10;
        } else {
            printf("EQUILIBRIO\n");
            (*u_arresto)  = 0.2;
            delta_u = 0.02;
        }  
        for(u = 0, i = 0; i <= 10; u += delta_u, i++){
            // printf("valor de U = %lf\n", u);
            adk = 0;
            for(int k = 0; k < MAX_K; k++, adk += dk){
                alfa_dk = (2.0 * adk * adk) / sk_final[k];
                // printf("alfa_dk = %lf\n", alfa_dk);
                sks[i][k] = sk_initial[k] * exp(-u * alfa_dk) + sk_final[k] * (1.0 - exp(-u * alfa_dk));
            }
            // printf("%d\n", i);
        }       
        free(sk_initial);
        free(sk_final);
        save_sks_in_file(sks, 11, MAX_K);
        free_memory_2D(sks, 11);
    } 
}

void save_sks_in_file(double **arr, int nRen, int nCol){
    FILE *f;

    f = fopen("sks.dat", "w");
    if(f){    
        for(int i = 0; i < nRen; i++){
            for(int j = 0; j < nCol; j++){
                fprintf(f, "%lf\n", arr[i][j]);
            }
        }
        fclose(f);
    }
}

