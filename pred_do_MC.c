//
//  red.c
//  
//
//  Created by Angel Gonzalez Mendez on 11/05/13.
//	Creacion de Redes Porosas, con sembrado y batido en un origen dinamico
//
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "pred_MC.h"

/****** elige revisado ***********/
void Gen_Rand(rand_t *myrand) {
    int i, jj = 0;
    myrand->general = 0;

    for (i = 0; i < Nknudt; i++) {
        myrand->v_255[i] = (i + 1) % Nknudt;
        myrand->knudt[i] = rand();
        myrand->v_122[i] = (i + 120) % Nknudt;
    }

    for (i = 0; i < 10000; i++) {
        myrand->knudt[jj] = myrand->knudt[jj] + myrand->knudt[myrand->v_122[jj]];
        jj = myrand->v_255[jj];
    }
}

/***************  Metodo del trapecio dada una Funcion*/

/*Funcion a integrar*/
double f(double Ri, double xrb, double sigma) {
    //printf("pi=%f\n",Pi);
    return (1 / (sqrt(2 * Pi) * sigma))*(exp((-1.0 * (pow((Ri - xrb), 2))) / (2 * pow(sigma, 2))));
}

void Trapecio(double n, double *Area, double *Translape, double *Ri, double xrb, double sigma) {
    double x = 0;
    double sup = xrb + (3 * sigma), inf = *Ri;
    double h = 0;
    double So = 0.0;
    int signo = 1;

    if (inf >= sup) {
        inf = inf + 20;
        sup = sup - 10;
        signo = -1;
    }
    h = (sup - inf) / n;
    for (x = inf + h; x < sup; x += h)
        So += f(x, xrb, sigma);
    (*Area) = (signo * h * (f(inf, xrb, sigma) + 2 * So + f(sup, xrb, sigma)) / 2.0);
    (*Translape) = (*Area) * 2;
}

/*Punto de Interseccion*/
double punt_inters(double xrs, double xrb) {
    double ri;
    double Rs = xrs, Rb = xrb;

    ri = ((pow(Rs, 2) - pow(Rb, 2)) / (2 * (Rs - Rb)));
    return (ri);
}

void Red_Translape(double xms, double xmb, double sigma, double *Translape, double *Ri) {

    double n = 999;
    double xrs = xms, xrb = xmb;
    double Area;

    *Ri = punt_inters(xrs, xrb);
    printf("\n ================================================");
    printf("\n Punto de Interseccion = %f", *Ri);
    printf("\n ================================================");
    Trapecio(n, &Area, Translape, Ri, xrb, sigma);
    printf("\n ================================================");
    printf("\n Translape = %5.6f", *Translape);
    printf("\n ================================================\n");
}

/*Funciones Creacion RED*/
void randomizar(int *randi, int *j_ran) {
    int i_ran, j_ran2;
    int num_vueltas = 1000000;

    *j_ran = 0;

    for (i_ran = 0; i_ran < 55; i_ran++)
        randi[i_ran] = (int) rand();

    for (j_ran2 = 0; j_ran2 < num_vueltas; j_ran2++) {
        (*j_ran) = ((*j_ran) + 1) % 55; // mas1_55[(*j_ran)];
        randi[(*j_ran)] = randi[(*j_ran)] + randi[((*j_ran) + 31) % 55]; // randi[mas31_55[(*j_ran)]];
    }
}

RED crea_matriz_poros(int m) {
    int i, j;
    RED arreglo;
    arreglo = (RED) calloc(m, sizeof (elem_i**));
    for (i = 0; i < m; i++) {
        arreglo[i] = (elem_i **) calloc(m, sizeof (elem_i*));
        for (j = 0; j < m; j++)
            arreglo[i][j] = (elem_i *) calloc(m, sizeof (elem_i));
    }
    return (arreglo);
}

// Funciones para calcular limites

double ran01(int randi[55], int *j_ran) {
    (*j_ran) = ((*j_ran) + 1) % 55;
    randi[(*j_ran)] = randi[(*j_ran)] + randi[((*j_ran) + 31) % 55];
    return (((double) randi[(*j_ran)] + long_max) / dos_long_max);
}

//Elige al azar una X y Y de la distribucion gaussiana y evalua la X.
//sale de la funcion cuando el valor evaluado es mayor que el Y al azar

float Sortedi(double delta, double lim_inf, double altura, double xm, double sigma, int randi[55], int *j_ran) {
    double xsorteo, ysorteo, yvalor;
    int resultado;
    resultado = 0;
    while (resultado == 0) {
        xsorteo = ran01(randi, j_ran) * delta + lim_inf;
        ysorteo = ran01(randi, j_ran) * altura;
        yvalor = altura * exp(-(xsorteo - xm)*(xsorteo - xm) / (2.0 * sigma * sigma)); /* Distrib. Gaussiana */
        if (yvalor > ysorteo)
            resultado = 1;
    }
    return ((float) (xsorteo));
}

double calcula_lim_inf_B(double sigma, double xmb) {
    return (xmb - 3 * sigma);
}

double calcula_lim_sup_B(double sigma, double xmb) {
    return (xmb + 3 * sigma);
}

double calcula_deltaB(double lim_sup, double lim_inf) {
    return (lim_sup - lim_inf);
}

double calcula_lim_inf_S(double sigma, double xms) {
    return (xms - 3 * sigma);
}

double calcula_lim_sup_S(double sigma, double xms) {
    return (xms + 3 * sigma);
}

double calcula_deltaS(double lim_sup, double lim_inf) {
    return (lim_sup - lim_inf);
}

/* Genera la red y coloca los sitios de mayor a menor*/
void Ini_red_Sitios_y_Enlaces(RED red, int L, double xmb, double xms, double sigma, double alturaB, double alturaS, double f0, int randi[55], int *j_ran) {

    double lim_inf_B, lim_inf_S, deltaB, deltaS;
    int i, j, k = 0;
    double rxb, ryb, rzb;
    double tini, tfin;
    int cont0 = 0;


    //Lista y red sitios
    tini = time(NULL);
    /**/

    /*paralelización*/


    lim_inf_B = calcula_lim_inf_B(sigma, xmb);
    lim_inf_S = calcula_lim_inf_S(sigma, xms);

    deltaB = calcula_deltaB(calcula_lim_sup_B(sigma, xmb), lim_inf_B);
    deltaS = calcula_deltaS(calcula_lim_sup_S(sigma, xms), lim_inf_S);

    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++)
            for (k = 0; k < L; k++) {
                red[i][j][k].sitio.radio = Sortedi(deltaS, lim_inf_S, alturaS, xms, sigma, randi, j_ran);

                if (ran01(randi, j_ran) > f0)
                    rxb = Sortedi(deltaB, lim_inf_B, alturaB, xmb, sigma, randi, j_ran);
                else {
                    cont0++;
                    rxb = 0.0;
                }
                if (ran01(randi, j_ran) > f0)
                    ryb = Sortedi(deltaB, lim_inf_B, alturaB, xmb, sigma, randi, j_ran);
                else {
                    cont0++;
                    ryb = 0.0;
                }
                if (ran01(randi, j_ran) > f0)
                    rzb = Sortedi(deltaB, lim_inf_B, alturaB, xmb, sigma, randi, j_ran);
                else {
                    cont0++;
                    rzb = 0.0;
                }


                red[i][j][k].xb.radio = rxb;
                red[i][j][k].yb.radio = ryb;
                red[i][j][k].zb.radio = rzb;

                red[i][j][k].sitio.radio_lleno = 0.0;
                red[i][j][k].xb.radio_lleno = 0.0;
                red[i][j][k].yb.radio_lleno = 0.0;
                red[i][j][k].zb.radio_lleno = 0.0;

            }

    tfin = time(NULL);
    printf("Tiempo de inicialización de red:%f\n", tfin - tini);
}

double val_enlace(int i, int j, int k, int dir, RED red, int L) {
    double result = 0.0;

    switch (dir) {
        case IZQUIERDA: result = red[i][j][k].xb.radio;
            break;
        case ABAJO:result = red[i][j][k].yb.radio;
            break;
        case FONDO: result = red[i][j][k].zb.radio;
            break;
        case DERECHA: result = red[(i + 1) % L][j][k].xb.radio;
            break;
        case ARRIBA: result = red[i][(j + 1) % L][k].yb.radio;
            break;
        case FRENTE: result = red[i][j][(k + 1) % L].zb.radio;
            break;
    }
    return result;
}

RED Inicializa_Red_RG_CV(int L, double xmb_i, double xms, double sigma, double alturaB, double alturaS, double f0, int randi[55], int * j_ran) {

    long timeini, timefin;
    RED red;

    red = crea_matriz_poros(L);
    Ini_red_Sitios_y_Enlaces(red, L, xmb_i, xms, sigma, alturaB, alturaS, f0, randi, j_ran);

    return red;
}

//*****************Red_guarda************************************

void Save_Red(RED red, cadena outfile_red, cadena ext, int L) {
    FILE *filered;
    int i, j, k;
    cadena archivo_red;

    strcpy(archivo_red, outfile_red);
    strcat(archivo_red, ext);
    if ((filered = fopen(archivo_red, "wb")) == NULL) {
        printf("Error al guardar la red de poros... %s\n", archivo_red);
        exit(0);
    }
    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++)
            for (k = 0; k < L; k++) {
                fwrite(&red[i][j][k].sitio.radio, sizeof (double), 1, filered);
                fwrite(&red[i][j][k].xb.radio, sizeof (double), 1, filered);
                fwrite(&red[i][j][k].yb.radio, sizeof (double), 1, filered);
                fwrite(&red[i][j][k].zb.radio, sizeof (double), 1, filered);
            }
    fclose(filered);
    printf("Archivo de red guardado: %s \n", archivo_red);
}

void Red_guarda(RED red, int L, cadena r_red) {
    Save_Red(red, r_red, ".bin", L);
}


//**********  funciones de validacion de redes
//************************************************
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//************************************************	Num_viola_tipo2()

//      regresa el numero de violaciones tipo 1 en cada unión enlace-sitio

int Num_violacSitios(RED red, int L, int *sit_viol) {
    int i, j, k, nviolT = 0, flag = 0;
    double enlace_izquierdo, enlace_abajo, enlace_frente, enlace_derecha, enlace_arriba, enlace_atras, sitio;

    *sit_viol = 0;
    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++)
            for (k = 0; k < L; k++) {
                flag = 0;
                enlace_izquierdo = val_enlace(i, j, k, IZQUIERDA, red, L);
                enlace_arriba = val_enlace(i, j, k, ARRIBA, red, L);
                enlace_frente = val_enlace(i, j, k, FRENTE, red, L);
                enlace_derecha = val_enlace(i, j, k, DERECHA, red, L);
                enlace_abajo = val_enlace(i, j, k, ABAJO, red, L);
                enlace_atras = val_enlace(i, j, k, FONDO, red, L);
                sitio = red[i][j][k].sitio.radio;
                if (enlace_izquierdo > sitio) { //printf("ERROR ei:%f > s:%f\n",enlace_izquierdo,sitio);
                    nviolT++;
                    flag = 1;
                }
                if (enlace_arriba > sitio) { //printf("ERROR earr:%f > s:%f\n",enlace_arriba,sitio);
                    nviolT++;
                    flag = 1;
                }
                if (enlace_frente > sitio) { //printf("ERROR ef:%f > s:%f\n",enlace_frente,sitio);
                    nviolT++;
                    flag = 1;
                }
                if (enlace_derecha > sitio) { //printf("ERROR ed:%f > s:%f\n",enlace_derecha,sitio);
                    nviolT++;
                    flag = 1;
                }
                if (enlace_abajo > sitio) { //printf("ERROR eab:%f > s:%f\n",enlace_abajo,sitio);
                    nviolT++;
                    flag = 1;
                }
                if (enlace_atras > sitio) { //printf("ERROR eat:%f > s:%f\n",enlace_atras,sitio);
                    nviolT++;
                    flag = 1;
                }
                if (flag == 1) (*sit_viol)++;
            }
    return nviolT;
}

int errores_geom_A(double en_derecha, double en_abajo, double en_izquierdo, double en_arriba, double en_frente, double en_fondo, double sitio, double A) {
    int viol_geom = 0;

    if ((A * (pow(en_derecha, 2.0) + pow(en_frente, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_abajo, 2.0) + pow(en_frente, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_izquierdo, 2.0) + pow(en_frente, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_arriba, 2.0) + pow(en_frente, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_derecha, 2.0) + pow(en_fondo, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_abajo, 2.0) + pow(en_fondo, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_izquierdo, 2.0) + pow(en_fondo, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_arriba, 2.0) + pow(en_fondo, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_arriba, 2.0) + pow(en_derecha, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_abajo, 2.0) + pow(en_derecha, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_arriba, 2.0) + pow(en_izquierdo, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_abajo, 2.0) + pow(en_izquierdo, 2.0))) > pow(sitio, 2.0)) viol_geom++;

    if (en_derecha > sitio) viol_geom++;
    if (en_abajo > sitio) viol_geom++;
    if (en_izquierdo > sitio) viol_geom++;
    if (en_arriba > sitio) viol_geom++;
    if (en_frente > sitio) viol_geom++;
    if (en_fondo > sitio) viol_geom++;

    return viol_geom;
}

int Num_violacGeom_A(RED red, int L, double A, int *sit_viol) {
    int i, j, k, nviolgeom = 0, nviolG = 0, nsviol = 0;
    double enlace_izquierdo, enlace_abajo, enlace_frente, enlace_derecha, enlace_arriba, enlace_atras, sitio;

    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++)
            for (k = 0; k < L; k++) {
                enlace_izquierdo = val_enlace(i, j, k, IZQUIERDA, red, L);
                enlace_arriba = val_enlace(i, j, k, ARRIBA, red, L);
                enlace_frente = val_enlace(i, j, k, FRENTE, red, L);
                enlace_derecha = val_enlace(i, j, k, DERECHA, red, L);
                enlace_abajo = val_enlace(i, j, k, ABAJO, red, L);
                enlace_atras = val_enlace(i, j, k, FONDO, red, L);
                sitio = red[i][j][k].sitio.radio;
                nviolgeom = errores_geom_A(enlace_derecha, enlace_abajo, enlace_izquierdo, enlace_arriba, enlace_frente, enlace_atras, sitio, A);
                if (nviolgeom > 0)
                    nsviol++;
                nviolG += nviolgeom;
            }
    *sit_viol = nsviol;
    return nviolG; //cuenta doble los enlaces
}

int Num_viola_tipo2(RED red, int L, int *enl_vacios) {
    int i, j, k;
    double sitio_cuadrado;
    double en1c, en2c, en3c, en4c, en5c, en6c;
    int viola = 0;

    *enl_vacios = 0;
    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++)
            for (k = 0; k < L; k++) {
                sitio_cuadrado = pow(red[i][j][k].sitio.radio, 2.0);
                en1c = pow(val_enlace(i, j, k, ARRIBA, red, L), 2.0);
                en2c = pow(val_enlace(i, j, k, ABAJO, red, L), 2.0);
                en3c = pow(val_enlace(i, j, k, IZQUIERDA, red, L), 2.0);
                en4c = pow(val_enlace(i, j, k, DERECHA, red, L), 2.0);
                en5c = pow(val_enlace(i, j, k, FRENTE, red, L), 2.0);
                en6c = pow(val_enlace(i, j, k, FONDO, red, L), 2.0);
                if (((en1c + en3c) > sitio_cuadrado) ||
                        ((en1c + en4c) > sitio_cuadrado) ||
                        ((en1c + en5c) > sitio_cuadrado) ||
                        ((en1c + en6c) > sitio_cuadrado) ||
                        ((en2c + en3c) > sitio_cuadrado) ||
                        ((en2c + en4c) > sitio_cuadrado) ||
                        ((en2c + en5c) > sitio_cuadrado) ||
                        ((en2c + en6c) > sitio_cuadrado) ||
                        ((en3c + en1c) > sitio_cuadrado) ||
                        ((en3c + en2c) > sitio_cuadrado) ||
                        ((en3c + en5c) > sitio_cuadrado) ||
                        ((en3c + en6c) > sitio_cuadrado) ||
                        ((en4c + en1c) > sitio_cuadrado) ||
                        ((en4c + en2c) > sitio_cuadrado) ||
                        ((en4c + en5c) > sitio_cuadrado) ||
                        ((en4c + en6c) > sitio_cuadrado) ||
                        ((en5c + en1c) > sitio_cuadrado) ||
                        ((en5c + en2c) > sitio_cuadrado) ||
                        ((en5c + en3c) > sitio_cuadrado) ||
                        ((en5c + en4c) > sitio_cuadrado) ||
                        ((en6c + en1c) > sitio_cuadrado) ||
                        ((en6c + en2c) > sitio_cuadrado) ||
                        ((en6c + en3c) > sitio_cuadrado) ||
                        ((en6c + en4c) > sitio_cuadrado))
                    viola++;

                if (val_enlace(i, j, k, ABAJO, red, L) == 0.0) (*enl_vacios)++;
                if (val_enlace(i, j, k, IZQUIERDA, red, L) == 0.0) (*enl_vacios)++;
                if (val_enlace(i, j, k, FONDO, red, L) == 0.0) (*enl_vacios)++;
            }
    return (viola);
}

/***************************************************************************************/
/*************************** Funciones para Batidos ************************************/

/***************************************************************************************/

double sitizq(RED red, int X, int Y, int Z, int L) {
    return (red[(X - 1 + L) % L][Y][Z].sitio.radio);
}

double sitabajo(RED red, int X, int Y, int Z, int L) {
    return (red[X][(Y - 1 + L) % L][Z].sitio.radio);
}

double sitatras(RED red, int X, int Y, int Z, int L) {
    return (red[X][Y][(Z - 1 + L) % L].sitio.radio);
}

//Regresa el radio del sitio actual

double sitioac(RED red, int i, int j, int k) {
    return (red[i][j][k].sitio.radio);
};

/*cuenta violaciones geometricas en sitios*/
int error_geom_sitios_A(RED red, int L, int X, int Y, int Z, double A) {
    int nviolgeom = 0;
    double enlace_izquierdo, enlace_abajo, enlace_frente, enlace_derecha, enlace_arriba, enlace_atras, sitio;
    enlace_izquierdo = val_enlace(X, Y, Z, IZQUIERDA, red, L);
    enlace_arriba = val_enlace(X, Y, Z, ARRIBA, red, L);
    enlace_frente = val_enlace(X, Y, Z, FRENTE, red, L);
    enlace_derecha = val_enlace(X, Y, Z, DERECHA, red, L);
    enlace_abajo = val_enlace(X, Y, Z, ABAJO, red, L);
    enlace_atras = val_enlace(X, Y, Z, FONDO, red, L);
    sitio = red[X][Y][Z].sitio.radio;
    nviolgeom = errores_geom_A(enlace_derecha, enlace_abajo, enlace_izquierdo, enlace_arriba, enlace_frente, enlace_atras, sitio, A);
    return nviolgeom;
}

int viol_sitio_sitizq_A(double en_frente, double en_abajo, double en_atras, double en_arriba, double en_frente_sitizq, double en_abajo_sitizq, double en_atras_sitizq, double en_arriba_sitizq, double en_izquierdo, double sitio, double sitio_izq, double A) {
    int viol_geom = 0;
    if ((A * (pow(en_frente, 2.0) + pow(en_izquierdo, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_abajo, 2.0) + pow(en_izquierdo, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_atras, 2.0) + pow(en_izquierdo, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_arriba, 2.0) + pow(en_izquierdo, 2.0))) > pow(sitio, 2.0)) viol_geom++;
    if ((A * (pow(en_frente_sitizq, 2.0) + pow(en_izquierdo, 2.0))) > pow(sitio_izq, 2.0)) viol_geom++;
    if ((A * (pow(en_abajo_sitizq, 2.0) + pow(en_izquierdo, 2.0))) > pow(sitio_izq, 2.0)) viol_geom++;
    if ((A * (pow(en_atras_sitizq, 2.0) + pow(en_izquierdo, 2.0))) > pow(sitio_izq, 2.0)) viol_geom++;
    if ((A * (pow(en_arriba_sitizq, 2.0) + pow(en_izquierdo, 2.0))) > pow(sitio_izq, 2.0)) viol_geom++;
    if (en_izquierdo > sitio) viol_geom++;
    if (en_izquierdo > sitio_izq) viol_geom++;

    return viol_geom;
}

/*junto con viol_sitio_sitizq regresan el numero de errores geometricos que involucran a un enlace Bx IZQUIERDA*/
int cuenta_viogeomBx_A(RED red, int L, int X, int Y, int Z, double A) {
    int violaciones = 0;
    double enlace_izquierdo, enlace_arriba, enlace_frente, enlace_abajo, enlace_atras, sitio;
    double enlace_arriba_sitizq, enlace_frente_sitizq, enlace_abajo_sitizq, enlace_atras_sitizq, sitio_izquierdo;

    enlace_izquierdo = val_enlace(X, Y, Z, IZQUIERDA, red, L);
    enlace_arriba = val_enlace(X, Y, Z, ARRIBA, red, L);
    enlace_frente = val_enlace(X, Y, Z, FRENTE, red, L);
    enlace_abajo = val_enlace(X, Y, Z, ABAJO, red, L);
    enlace_atras = val_enlace(X, Y, Z, FONDO, red, L);

    enlace_arriba_sitizq = val_enlace((X - 1 + L) % L, Y, Z, ARRIBA, red, L);
    enlace_frente_sitizq = val_enlace((X - 1 + L) % L, Y, Z, FRENTE, red, L);
    enlace_abajo_sitizq = val_enlace((X - 1 + L) % L, Y, Z, ABAJO, red, L);
    enlace_atras_sitizq = val_enlace((X - 1 + L) % L, Y, Z, FONDO, red, L);
    sitio = sitioac(red, X, Y, Z);
    sitio_izquierdo = sitizq(red, X, Y, Z, L);

    violaciones = viol_sitio_sitizq_A(enlace_frente, enlace_abajo, enlace_atras, enlace_arriba, enlace_frente_sitizq, enlace_abajo_sitizq, enlace_atras_sitizq,
            enlace_arriba_sitizq, enlace_izquierdo, sitio, sitio_izquierdo, A);
    return violaciones;
}

int viol_sitio_sitabajo_A(double en_frente, double en_izquierdo, double en_atras, double en_derecha, double en_frente_sitabajo, double en_izquierdo_sitabajo, double en_atras_sitabajo, double en_derecha_sitabajo, double en_abajo, double sitio, double sitio_abajo, double A) {
    int viol_geom = 0;
    if ((A * (pow(en_frente, 2.0) + pow(en_abajo, 2.0)) > pow(sitio, 2.0))) viol_geom++;
    if ((A * (pow(en_izquierdo, 2.0) + pow(en_abajo, 2.0)) > pow(sitio, 2.0))) viol_geom++;
    if ((A * (pow(en_atras, 2.0) + pow(en_abajo, 2.0)) > pow(sitio, 2.0))) viol_geom++;
    if ((A * (pow(en_derecha, 2.0) + pow(en_abajo, 2.0)) > pow(sitio, 2.0))) viol_geom++;
    if ((A * (pow(en_frente_sitabajo, 2.0) + pow(en_abajo, 2.0))) > pow(sitio_abajo, 2.0)) viol_geom++;
    if ((A * (pow(en_izquierdo_sitabajo, 2.0) + pow(en_abajo, 2.0))) > pow(sitio_abajo, 2.0)) viol_geom++;
    if ((A * (pow(en_atras_sitabajo, 2.0) + pow(en_abajo, 2.0))) > pow(sitio_abajo, 2.0)) viol_geom++;
    if ((A * (pow(en_derecha_sitabajo, 2.0) + pow(en_abajo, 2.0))) > pow(sitio_abajo, 2.0)) viol_geom++;
    if (en_abajo > sitio) viol_geom++;
    if (en_abajo > sitio_abajo) viol_geom++;


    return viol_geom;
}

/*junto con viol_sitio_sitabajo regresan el numero de errores geometricos que involucran a un enlace By ABAJO*/
int cuenta_viogeomBy_A(RED red, int L, int X, int Y, int Z, double A) {
    int violaciones = 0;
    double enlace_abajo, enlace_izquierdo, enlace_frente, enlace_derecha, enlace_atras, sitio, sitio_abajo;
    double enlace_izquierdo_sitabajo, enlace_frente_sitabajo, enlace_derecha_sitabajo, enlace_atras_sitabajo;

    enlace_izquierdo = val_enlace(X, Y, Z, IZQUIERDA, red, L);
    enlace_frente = val_enlace(X, Y, Z, FRENTE, red, L);
    enlace_derecha = val_enlace(X, Y, Z, DERECHA, red, L);
    enlace_atras = val_enlace(X, Y, Z, FONDO, red, L);
    enlace_abajo = val_enlace(X, Y, Z, ABAJO, red, L);

    enlace_izquierdo_sitabajo = val_enlace(X, (Y - 1 + L) % L, Z, IZQUIERDA, red, L);
    enlace_frente_sitabajo = val_enlace(X, (Y - 1 + L) % L, Z, FRENTE, red, L);
    enlace_derecha_sitabajo = val_enlace(X, (Y - 1 + L) % L, Z, DERECHA, red, L);
    enlace_atras_sitabajo = val_enlace(X, (Y - 1 + L) % L, Z, FONDO, red, L);
    sitio = sitioac(red, X, Y, Z);
    sitio_abajo = sitabajo(red, X, Y, Z, L);

    violaciones = viol_sitio_sitabajo_A(enlace_frente, enlace_izquierdo, enlace_atras, enlace_derecha, enlace_frente_sitabajo, enlace_izquierdo_sitabajo,
            enlace_atras_sitabajo, enlace_derecha_sitabajo, enlace_abajo, sitio, sitio_abajo, A);
    return violaciones;
}

int viol_sitio_sitatras_A(double en_derecha, double en_abajo, double en_izquierdo, double en_arriba, double en_derecha_sitatras, double en_abajo_sitatras, double en_izquierdo_sitatras, double en_arriba_sitatras, double en_atras, double sitio, double sitio_atras, double A) {
    int viol_geom = 0;
    if ((A * (pow(en_derecha, 2.0) + pow(en_atras, 2.0)) > pow(sitio, 2.0))) viol_geom++;
    if ((A * (pow(en_abajo, 2.0) + pow(en_atras, 2.0)) > pow(sitio, 2.0))) viol_geom++;
    if ((A * (pow(en_izquierdo, 2.0) + pow(en_atras, 2.0)) > pow(sitio, 2.0))) viol_geom++;
    if ((A * (pow(en_arriba, 2.0) + pow(en_atras, 2.0)) > pow(sitio, 2.0))) viol_geom++;
    if ((A * (pow(en_derecha_sitatras, 2.0) + pow(en_atras, 2.0))) > pow(sitio_atras, 2.0)) viol_geom++;
    if ((A * (pow(en_abajo_sitatras, 2.0) + pow(en_atras, 2.0))) > pow(sitio_atras, 2.0)) viol_geom++;
    if ((A * (pow(en_izquierdo_sitatras, 2.0) + pow(en_atras, 2.0))) > pow(sitio_atras, 2.0)) viol_geom++;
    if ((A * (pow(en_arriba_sitatras, 2.0) + pow(en_atras, 2.0))) > pow(sitio_atras, 2.0)) viol_geom++;
    if (en_atras > sitio) viol_geom++;
    if (en_atras > sitio_atras) viol_geom++;

    return viol_geom;
}

/*junto con viol_sitio_sitatras regresan el numero de errores geometricos que involucran a un enlace Bz*/
int cuenta_viogeomBz_A(RED red, int L, int X, int Y, int Z, double A) {
    int violaciones = 0;
    double enlace_atras, enlace_izquierdo, enlace_arriba, enlace_derecha, enlace_abajo, sitio, sitio_atras;
    double enlace_derecha_sitatras, enlace_abajo_sitatras, enlace_izquierdo_sitatras, enlace_arriba_sitatras;

    enlace_izquierdo = val_enlace(X, Y, Z, IZQUIERDA, red, L);
    enlace_arriba = val_enlace(X, Y, Z, ARRIBA, red, L);
    enlace_derecha = val_enlace(X, Y, Z, DERECHA, red, L);
    enlace_abajo = val_enlace(X, Y, Z, ABAJO, red, L);
    enlace_atras = val_enlace(X, Y, Z, FONDO, red, L);

    enlace_izquierdo_sitatras = val_enlace(X, Y, (Z - 1 + L) % L, IZQUIERDA, red, L);
    enlace_arriba_sitatras = val_enlace(X, Y, (Z - 1 + L) % L, ARRIBA, red, L);
    enlace_derecha_sitatras = val_enlace(X, Y, (Z - 1 + L) % L, DERECHA, red, L);
    enlace_abajo_sitatras = val_enlace(X, Y, (Z - 1 + L) % L, ABAJO, red, L);
    sitio = sitioac(red, X, Y, Z);
    sitio_atras = sitatras(red, X, Y, Z, L);
    violaciones = viol_sitio_sitatras_A(enlace_derecha, enlace_abajo, enlace_izquierdo, enlace_arriba, enlace_derecha_sitatras, enlace_abajo_sitatras,
            enlace_izquierdo_sitatras, enlace_arriba_sitatras, enlace_atras, sitio, sitio_atras, A);
    return violaciones;
}


//Intercambia dos sitios o enlaces

void swap(double *X1, double *X2) {
    double aux;
    aux = *X1;
    *X1 = *X2;
    *X2 = aux;
}

//Bate una red considerando restricciones geométricas Tipo 2

void batidotEq_MC_Seq(RED red, int L, int randi[55], int *j_ran) {
    int i1, j1, k1, i2, j2, k2, e1t, e2t, control, iLE;
    int p1, p2, nswap;
    //posicion p;

    for (nswap = 0; nswap < 3 * L * L * L; nswap++) {
        //Genera un numero aleatorio entre 0 y 11
        // Para cambiar sitios y enlaces con la misma prob /
        control = (int) (12.0 * ran01(randi, j_ran));

        i1 = (int) (L * ran01(randi, j_ran));
        j1 = (int) (L * ran01(randi, j_ran));
        k1 = (int) (L * ran01(randi, j_ran));
        i2 = (int) (L * ran01(randi, j_ran));
        j2 = (int) (L * ran01(randi, j_ran));
        k2 = (int) (L * ran01(randi, j_ran));


        //printf("c: %d ",control);
        switch (control) {
            case 0: //sitio - sitio
            {
                e1t = error_geom_sitios_A(red, L, i1, j1, k1, 1.0) + error_geom_sitios_A(red, L, i2, j2, k2, 1.0);
                ;
                swap(&(red[i1][j1][k1].sitio.radio), &(red[i2][j2][k2].sitio.radio)); // Cambio los sitios
                e2t = error_geom_sitios_A(red, L, i1, j1, k1, 1.0) + error_geom_sitios_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&(red[i1][j1][k1].sitio.radio), &(red[i2][j2][k2].sitio.radio)); // Re-Cambio los sitios
                break;
            }
            case 1: // enlace_x - enlace_x'	IZQ-IZQ
            {
                e1t = cuenta_viogeomBx_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBx_A(red, L, i2, j2, k2, 1.0);
                swap(&(red[i1][j1][k1].xb.radio), &(red[i2][j2][k2].xb.radio)); // Cambio los enlace
                e2t = cuenta_viogeomBx_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBx_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&red[i1][j1][k1].xb.radio, &red[i2][j2][k2].xb.radio); //Re-Cambio los Bond
                break;
            }
            case 2: // enlace_x - enlace_y'	IZQ-ABA
            {
                e1t = cuenta_viogeomBx_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBy_A(red, L, i2, j2, k2, 1.0);
                swap(&red[i1][j1][k1].xb.radio, &red[i2][j2][k2].yb.radio); // Cambio los enlace
                e2t = cuenta_viogeomBx_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBy_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&red[i1][j1][k1].xb.radio, &red[i2][j2][k2].yb.radio); //Re-Cambio los Bond
                break;
            }
            case 3: // enlace_x - enlace_z'	IZQ-FON
            {
                e1t = cuenta_viogeomBx_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBz_A(red, L, i2, j2, k2, 1.0);
                swap(&red[i1][j1][k1].xb.radio, &red[i2][j2][k2].zb.radio); // Cambio los enlace
                e2t = cuenta_viogeomBx_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBz_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&red[i1][j1][k1].xb.radio, &red[i2][j2][k2].zb.radio); //Re-Cambio los Bond
                break;
            }
            case 4: // enlace_y - enlace_x'	ABA-IZQ
            {
                e1t = cuenta_viogeomBy_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBx_A(red, L, i2, j2, k2, 1.0);
                swap(&red[i1][j1][k1].yb.radio, &red[i2][j2][k2].xb.radio); // Cambio los enlace
                e2t = cuenta_viogeomBy_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBx_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&red[i1][j1][k1].yb.radio, &red[i2][j2][k2].xb.radio); //Re-Cambio los Bond

                break;
            }
            case 5: // enlace_y - enlace_y'   ABA-ABA
            {
                e1t = cuenta_viogeomBy_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBy_A(red, L, i2, j2, k2, 1.0);
                swap(&red[i1][j1][k1].yb.radio, &red[i2][j2][k2].yb.radio); // Cambio los enlace
                e2t = cuenta_viogeomBy_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBy_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&red[i1][j1][k1].yb.radio, &red[i2][j2][k2].yb.radio); //Re-Cambio los Bond_x
                break;
            }
                //		case 6:  // enlace_y - enlace_z'   ABA-FON
                //		{
                //			e1t=cuenta_viogeomBy_A(red,L,i1,j1,k1,1.0)+cuenta_viogeomBz_A(red,L,i2,j2,k2,1.0);
                //			swap(&red[i1][j1][k1].yb.radio,&red[i2][j2][k2].zb.radio); // Cambio los enlace
                //			e2t=cuenta_viogeomBy_A(red,L,i1,j1,k1,1.0)+cuenta_viogeomBz_A(red,L,i2,j2,k2,1.0);
                //			if (e1t < e2t)
                //				swap(&red[i1][j1][k1].yb.radio,&red[i2][j2][k2].zb.radio); //Re-Cambio
                //			break;
                //		}
            case 7: // enlace_z - enlace_z'	FON-FON
            {
                e1t = cuenta_viogeomBz_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBz_A(red, L, i2, j2, k2, 1.0);
                swap(&red[i1][j1][k1].zb.radio, &red[i2][j2][k2].zb.radio); // Cambio los enlaces
                e2t = cuenta_viogeomBz_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBz_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&red[i1][j1][k1].zb.radio, &red[i2][j2][k2].zb.radio); //Re-Cambio los Bond_x

                break;
            }
            case 8: // enlace_z - enlace_y'	FON-ABA
            {
                e1t = cuenta_viogeomBz_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBy_A(red, L, i2, j2, k2, 1.0);
                swap(&red[i1][j1][k1].zb.radio, &red[i2][j2][k2].yb.radio); // Cambio los enlace
                e2t = cuenta_viogeomBz_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBy_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&red[i1][j1][k1].zb.radio, &red[i2][j2][k2].yb.radio); //Re-Cambio los Bond_x
                break;
            }
            case 9: // enlace_z - enlace_x'	FON-IZQ
            {
                e1t = cuenta_viogeomBz_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBx_A(red, L, i2, j2, k2, 1.0);
                swap(&red[i1][j1][k1].zb.radio, &red[i2][j2][k2].xb.radio); // Cambio los enlace
                e2t = cuenta_viogeomBz_A(red, L, i1, j1, k1, 1.0) + cuenta_viogeomBx_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&red[i1][j1][k1].zb.radio, &red[i2][j2][k2].xb.radio); //Re-Cambio los Bond_x
                break;
            }
            case 10: // site - site
            {
                e1t = error_geom_sitios_A(red, L, i1, j1, k1, 1.0) + error_geom_sitios_A(red, L, i2, j2, k2, 1.0);
                swap(&(red[i1][j1][k1].sitio.radio), &(red[i2][j2][k2].sitio.radio)); //Cambio los sitios
                e2t = error_geom_sitios_A(red, L, i1, j1, k1, 1.0) + error_geom_sitios_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&(red[i1][j1][k1].sitio.radio), &(red[i2][j2][k2].sitio.radio)); //Re-Cambio los sitios
                break;
            }
            case 11: // site - site
            {
                e1t = error_geom_sitios_A(red, L, i1, j1, k1, 1.0) + error_geom_sitios_A(red, L, i2, j2, k2, 1.0);
                swap(&(red[i1][j1][k1].sitio.radio), &(red[i2][j2][k2].sitio.radio)); //Cambio los sitios
                e2t = error_geom_sitios_A(red, L, i1, j1, k1, 1.0) + error_geom_sitios_A(red, L, i2, j2, k2, 1.0);
                if (e1t < e2t)
                    swap(&(red[i1][j1][k1].sitio.radio), &(red[i2][j2][k2].sitio.radio)); //Re-Cambio los sitios
                break;
            }
            default:
                break;
        }
    } // fin de un batido = 3*L*L*L swaps
}

/** Fucion principal Batidos **/
double NBatidos_T2_RG_MC_Seq(RED red, const int L, double nBATIDOS, unsigned long int *nBatidos_Real, int randi[55], int *j_ran) {
    int gsit_viola, sit_viola;
    int continuar, nviol_en_sit;
    double nbat;

    nviol_en_sit = Num_violacGeom_A(red, L, 1.0, &sit_viola); //Num_viola_tipo2(red,L);

    gsit_viola += sit_viola;

    nbat = 0.0;
    continuar = nBATIDOS == 0.0 ? gsit_viola > 0 : nbat <= nBATIDOS;
    continuar = 1;
    while (continuar) {
        batidotEq_MC_Seq(red, L, randi, j_ran);
        nbat = nbat + 1.0;
        if (fmod(nbat, 20.0f) == 0.0) {
            nviol_en_sit = Num_violacGeom_A(red, L, 1.0, &sit_viola); //Num_viola_tipo2(red,L);
            gsit_viola = sit_viola;
            *nBatidos_Real = nbat;
            printf("Batido No. %lu, con %d sitios con violaciones\n", *nBatidos_Real, gsit_viola);
        }

        continuar = nBATIDOS == 0.0 ? gsit_viola > 0 : nbat <= nBATIDOS;
        if (gsit_viola == 0)
            continuar = 0;
    }
    *nBatidos_Real = nbat;

    return (0);
}