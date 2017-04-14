#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "pred_MC.h"
#include "Imgs_Gplot_Rmol.h"

//#define NBAT 20.0 //10000.0
//#define BATIDOS(x) x*(4.0*L*L*L) //No. total de sitios + enlaces

int main(int argc, char *argv[]) {
    //Variables
    RED red;
    int L;
    double xmb = 26.0; //400.0;//26.0;
    double xms = 44.0; //42.0; la más difícil
    double sigma = 6.0; //50.0; //6.0; //20;
    double Translape, Ri;
    int k_aps;
    int NCLUSTERS, CLUSTER_SIZE;
    int NBatidos;
    unsigned long int nBatidos_Real;
    int nthreads;
    cadena name_red, name_3DGnuplot_conec, name_3DGnuplot_size, main_name, name_aux;
    int flag, i, j;

    //Variables Tiempo
    double tini1, tfin1, tfin2, tfin3;
    double tbat1, tbat2;

    int j_ran;
    int randi[55];

    int nviol_en_sit;
    int enl_vacios, sit_viola;

    double f0 = 0.5;

    enl_vacios = 0;
    NBatidos = 0;

    if (argc == 6) {//Batidos
        L = atoi(argv[1]);
        xmb = atof(argv[2]);
        xms = atof(argv[3]);
        sigma = atof(argv[4]);
        f0 = atof(argv[5]);
        NBatidos = 0;
    } else {
        fprintf(stderr, "\nParametros invalidos\n");
        fprintf(stderr, "\nCreacion\nUso: %s  L xmb xms sigma f0\ni.e. %s 50 26 46 6 0.4\n\n", argv[0], argv[0]);
        exit(-1);
    }

    //Creacion
    tini1 = time(NULL);

    /** Inicializacion **/
    srand((unsigned int) time(NULL));
    Gen_Rand(&global_rand);
    randomizar(randi, &j_ran);
    Red_Translape(xms, xmb, sigma, &Translape, &Ri);

    /** Creacion **/
    red = Inicializa_Red_RG_CV(L, xmb, xms, sigma, 1, 1, f0, randi, &j_ran); //sin violaciones tipo 2

    tfin1 = time(NULL);
    printf("\nTiempo de creacion con RG:%.2lf segs\n\n", tfin1 - tini1);

    nviol_en_sit = Num_violacSitios(red, L, &sit_viola); //Violaciones tipo 1
    printf("No.Violaciones en Sitios (Tipo 1): %d, Num Sit c/viol: %d\n", nviol_en_sit, sit_viola);
    nviol_en_sit = Num_violacGeom_A(red, L, 1.0, &sit_viola); //Num_viola_tipo2(red,L);
    printf("No.Violaciones geometricasA (Tipo 2): %d, Num Sit c/viol: %d\n", nviol_en_sit, sit_viola);
    nviol_en_sit = Num_viola_tipo2(red, L, &enl_vacios);
    printf("No.Violaciones geometricasB (Tipo 2):%d, enlaces vacios: %d\n", nviol_en_sit, enl_vacios);

    /*Guardar red sin batidos*/
    sprintf(main_name, "L%d_xmb%1.0f_xms%1.0f_s%1.0f_f0%1.0f_s_%dbat", L, xmb, xms, sigma, f0 * 10.0, NBatidos);
    sprintf(name_red, "archivos/%s", main_name);
    tfin2 = time(NULL);
    printf("\nTiempo total de creación mas Validaciones:%f segs\n\n", tfin2 - tini1);
    sprintf(name_3DGnuplot_conec, "Gplot3D_C_%s", main_name);
    sprintf(name_3DGnuplot_size, "Gplot3D_S_%s", main_name);
    tbat1 = time(NULL);

    ImgGnuplot3D_x_conectividad(red, name_3DGnuplot_conec, L, xms, sigma);
    genera_script_gplot3D_x_conec(name_3DGnuplot_conec);

    ImgGnuplot3D_x_size(red, name_3DGnuplot_size, L, xms, sigma);
    genera_script_gplot3D_x_size(name_3DGnuplot_size);



    NBatidos_T2_RG_MC_Seq(red, L, (double) NBatidos, &nBatidos_Real, randi, &j_ran);

    nviol_en_sit = Num_violacSitios(red, L, &sit_viola); //Violaciones tipo 1
    printf("\nDESPUES DE BATIDOS\nNo.Violaciones en Sitios antes Certificacion (Tipo 1): %d, Num Sit c/viol: %d\n", nviol_en_sit, sit_viola);
    nviol_en_sit = Num_violacGeom_A(red, L, 1.0, &sit_viola); //Num_viola_tipo2(red,L);
    printf("No.Violaciones geometricasA (Tipo 2): %d, Num Sit c/viol: %d\n", nviol_en_sit, sit_viola);
    nviol_en_sit = Num_viola_tipo2(red, L, &enl_vacios);
    printf("No.Violaciones geometricasB (Tipo 2):%d, enlaces vacios: %d\n", nviol_en_sit, enl_vacios);

    tbat2 = time(NULL);
    printf("\nTiempo de %lu Batidos:%6.0f segs\n\n", nBatidos_Real, tbat2 - tbat1);

    sprintf(main_name, "L%d_xmb%1.0f_xms%1.0f_s%1.0f_c_%lubat", L, xmb, xms, sigma, nBatidos_Real);
    sprintf(name_red, "archivos/%s", main_name);
    Red_guarda(red, L, name_red);
    sprintf(name_3DGnuplot_conec, "Gplot3D_C_%s", main_name);
    sprintf(name_3DGnuplot_size, "Gplot3D_S_%s", main_name);

    ImgGnuplot3D_x_conectividad(red, name_3DGnuplot_conec, L, xms, sigma);
    genera_script_gplot3D_x_conec(name_3DGnuplot_conec);

    ImgGnuplot3D_x_size(red, name_3DGnuplot_size, L, xms, sigma);
    genera_script_gplot3D_x_size(name_3DGnuplot_size);
    system("gnuplot script3Dsize");
    return 0;
}