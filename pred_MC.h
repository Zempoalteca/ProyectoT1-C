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

/*Origen*/
int origen[3];

//Numeros aleatorios

typedef struct {
    int knudt[Nknudt];
    int v_255[Nknudt], v_122[Nknudt], general;
} rand_t;

rand_t global_rand;

typedef struct {
    int xi, xf;
    int yi, yf;
    int zi, zf;
    int nx, ny, nz;
    int psize;
} partition_t;

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

/*Estructura para la creacion de la Red*/
typedef struct {
    //Manaejo de Threads
    int tid;
    int nthreads;
    //Areas de trabajo
    int xi, xf;
    int yi, yf;
    int zi, zf;
    int psize;
    int nx, ny, nz;
    int *gdistribution;
    int *gsize_distribution;
    //Manejo de las redes
    int *gsite_count;
    int *gsite_missing;
    int *lsite_count;
    //funcion
    int randi[55];
    int j_ran;
    //Numeros aleatorios
    rand_t trand;
} data_t;

/*Estructura para batidos de la Red*/
typedef struct {
    //Manaejo de Threads
    int tid;
    int nthreads;
    //Particiones y distribucion
    int *gdistribution;
    int *gsize_distribution;
    int npartitions;
    int *partitions;
    int partitions_per_thread;
    partition_t *ps;
    //funcion
    int randi[55];
    int j_ran;
    //Numeros aleatorios
    rand_t trand;
} databat_t;

//Funciones
void Red_Translape(double xms, double xmb, double sigma, double *Translape, double *Ri);
void Red_guarda(RED red, int L, cadena r_red);
void Gen_Rand(rand_t *myrand);

RED Inicializa_Red_RG_CV(int L, double xmb_i, double xms, double sigma, double alturaB, double alturaS, double f0, int randi[55], int * j_ran);
void randomizar(int *randi, int *j_ran);
double NBatidos_T2_RG_MC_Seq(RED red, const int L, double nBATIDOS, unsigned long int *nBatidos_Real, int randi[55], int *j_ran);

//Validaciones
int Num_violacSitios(RED red, int L, int *sit_viol);
int Num_violacGeom_A(RED red, int L, double A, int *sit_viol);
int Num_viola_tipo2(RED red, int L, int *enl_vacios);
double val_enlace(int i, int j, int k, int dir, RED red, int L);

#endif
