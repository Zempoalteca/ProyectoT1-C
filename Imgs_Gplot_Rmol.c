#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "pred_MC.h"
#include "Imgs_Gplot_Rmol.h"

void genera_script_gplot3D_x_conec(cadena name_corte3DGnuplot) {
    FILE *ap;

    ap = fopen("script3Dconec", "w");
    if (ap == NULL)
        printf("ERROR.... no se pudo crear archivo script3D\n");
    else {
        fprintf(ap, "set hidden3d\n");
        fprintf(ap, "set output \"%s.eps\"\n", name_corte3DGnuplot);
        fprintf(ap, "set terminal postscript eps enhanced color size 4,4.5 font 16\n");
        fprintf(ap, "splot \"Big_%s.xyz\" u 2:3:4 with points title \"5-6 bonds\",", name_corte3DGnuplot);
        fprintf(ap, " \"Medium_%s.xyz\" u 2:3:4 with points title \"3-4 bonds\",", name_corte3DGnuplot);
        fprintf(ap, " \"Small_%s.xyz\" u 2:3:4 with points title \"1-2 bonds\"\n\n", name_corte3DGnuplot);
        fclose(ap);
        printf("\nScript3D generado\n");
        system("gnuplot script3Dconec");
    }
}

void genera_script_gplot3D_x_size(cadena name_corte3DGnuplot) {
    FILE *ap;

    ap = fopen("script3Dsize", "w");
    if (ap == NULL)
        printf("ERROR.... no se pudo crear archivo script3D\n");
    else {
        fprintf(ap, "set hidden3d\n");
        fprintf(ap, "set output \"%s.eps\"\n", name_corte3DGnuplot);
        fprintf(ap, "set terminal postscript eps enhanced color size 4,4.5 font 16\n");
        fprintf(ap, "splot \"Big_%s.xyz\" u 2:3:4 with points title \"big\",", name_corte3DGnuplot);
        fprintf(ap, " \"Medium_%s.xyz\" u 2:3:4 with points title \"medium\",", name_corte3DGnuplot);
        fprintf(ap, " \"Small_%s.xyz\" u 2:3:4 with points title \"small\"\n\n", name_corte3DGnuplot);
        fclose(ap);
        printf("\nScript3D generado\n");
        system("gnuplot script3Dsize");
    }
}

int cuentaEnlaces(int i, int j, int k, int L, RED red) {
    int enlaces = 0;
    if (val_enlace(i, j, k, DERECHA, red, L) > 0.0) enlaces++;
    if (val_enlace(i, j, k, IZQUIERDA, red, L) > 0.0) enlaces++;
    if (val_enlace(i, j, k, ABAJO, red, L) > 0.0) enlaces++;
    if (val_enlace(i, j, k, ARRIBA, red, L) > 0.0) enlaces++;
    if (val_enlace(i, j, k, FONDO, red, L) > 0.0) enlaces++;
    if (val_enlace(i, j, k, FRENTE, red, L) > 0.0) enlaces++;

    return enlaces;
}

void ImgGnuplot3D_x_conectividad(RED red, cadena nom_arch, int L, double xms, double sigma) {
    int i, j, k, f, neG = 0, neM = 0, neS = 0, enl_vacios = 0, nviol_en_sit;
    FILE *ap_arch_small;
    FILE *ap_arch_medium;
    FILE *ap_arch_big;
    cadena nom_arch_small;
    cadena nom_arch_medium;
    cadena nom_arch_big;

    sprintf(nom_arch_small, "Small_%s.xyz", nom_arch);
    sprintf(nom_arch_medium, "Medium_%s.xyz", nom_arch);
    sprintf(nom_arch_big, "Big_%s.xyz", nom_arch);

    ap_arch_small = fopen(nom_arch_small, "w");
    ap_arch_medium = fopen(nom_arch_medium, "w");
    ap_arch_big = fopen(nom_arch_big, "w");

    fprintf(ap_arch_small, "tipo x y z\n");
    fprintf(ap_arch_medium, "tipo x y z\n");
    fprintf(ap_arch_big, "tipo x y z\n");

    printf("L= %d\n", L);
    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++)
            for (k = 0; k < L; k++) {
                f = cuentaEnlaces(i, j, k, L, red);
                if (f > 0) {
                    //printf("ne:%d ",f);
                    if (f <= 2) {
                        neS++;
                        fprintf(ap_arch_small, "S %d %d %d\n", i, j, k);
                    } else
                        if (f <= 4) {
                        neM++;
                        fprintf(ap_arch_medium, "H %d %d %d\n", i, j, k);
                    } else {
                        neG++;
                        fprintf(ap_arch_big, "N %d %d %d\n", i, j, k);
                    }
                }
            }
    fclose(ap_arch_small);
    fclose(ap_arch_medium);
    fclose(ap_arch_big);
    printf("\nArchivos 3D Creados para GnuPlot:\nArchivo1:%s\nArchivo2:%s\nArchivo3:%s\n", nom_arch_small, nom_arch_medium, nom_arch_big);
}

void ImgGnuplot3D_x_size(RED red, cadena nom_arch, int L, double xms, double sigma) {
    int i, j, k, f;
    FILE *ap_arch_small;
    FILE *ap_arch_medium;
    FILE *ap_arch_big;
    cadena nom_arch_small;
    cadena nom_arch_medium;
    cadena nom_arch_big;

    sprintf(nom_arch_small, "Small_%s.xyz", nom_arch);
    sprintf(nom_arch_medium, "Medium_%s.xyz", nom_arch);
    sprintf(nom_arch_big, "Big_%s.xyz", nom_arch);

    ap_arch_small = fopen(nom_arch_small, "w");
    ap_arch_medium = fopen(nom_arch_medium, "w");
    ap_arch_big = fopen(nom_arch_big, "w");

    fprintf(ap_arch_small, "tipo x y z\n");
    fprintf(ap_arch_medium, "tipo x y z\n");
    fprintf(ap_arch_big, "tipo x y z\n");

    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++)
            for (k = 0; k < L; k++) {
                f = cuentaEnlaces(i, j, k, L, red);
                if (f > 0) {
                    //printf("ne:%d ",f);
                    if (red[i][j][k].sitio.radio <= (xms - 0.43 * sigma))
                        fprintf(ap_arch_small, "S %d %d %d\n", i, j, k);
                    else
                        if (((xms - 0.43 * sigma) < red[i][j][k].sitio.radio) && (red[i][j][k].sitio.radio <= (xms + 0.43 * sigma)))
                        fprintf(ap_arch_medium, "H %d %d %d\n", i, j, k);
                    else
                        fprintf(ap_arch_big, "N %d %d %d\n", i, j, k);
                }
            }
    printf("\n");
    fclose(ap_arch_small);
    fclose(ap_arch_medium);
    fclose(ap_arch_big);
    printf("\nArchivos 3D Creados para GnuPlot:\nArchivo1:%s\nArchivo2:%s\nArchivo3:%s\n", nom_arch_small, nom_arch_medium, nom_arch_big);
}