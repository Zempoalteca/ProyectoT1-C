//
//  red.h
//  
//
//  Created by Angel Gonzalez Mendez on 11/05/13.
//
//

#ifndef _pred_h
#define _pred_h

//Constantes
#define ARRIBA    0
#define ABAJO     1
#define IZQUIERDA 2
#define DERECHA   3
#define FRENTE    4
#define FONDO  5

#define long_max  2147483648.0
#define dos_long_max 4294967296.0
#define Pi 3.1416  //6.2832
#define Nknudt 256

//Numeros aleatorios
typedef struct {
    int knudt[Nknudt];
    int v_255[Nknudt], v_122[Nknudt], general;
} rand_t;

rand_t global_rand;

//Tipos
struct poro {
    double radio;
    double radio_lleno;
    double esplib;
};
typedef struct poro radios;

struct elemento {
    radios sitio;
    radios xb; //arriba
    radios yb; //izuiqerda
    radios zb; //frente
};

typedef struct elemento elem_i;

typedef elem_i *** RED;

typedef char cadena[256];

//Funciones
void Red_Translape(double xms, double xmb, double sigma, double *Translape, double *Ri);
void Red_guarda(RED red, int L, cadena r_red);
void Gen_Rand(rand_t *myrand);

RED Inicializa_Red_RG_CV(int L, double xmb_i, double xms, double sigma, double alturaB, double alturaS, double f0, int randi[55], int * j_ran);
void randomizar(int *randi, int *j_ran);
double NBatidos_T2_RG_MC_Seq(RED red, const int L, unsigned long int *nBatidos_Real, int randi[55], int *j_ran);

//Validaciones
int Num_violacSitios(RED red, int L, int *sit_viol);
int Num_violacGeom_A(RED red, int L, double A, int *sit_viol);
int Num_viola_tipo2(RED red, int L, int *enl_vacios);
double val_enlace(int i, int j, int k, int dir, RED red, int L);

#endif