
#ifndef _imgs_gplot_rmol_h
#define _imgs_gplot_rmol_h

#define EJEX 1
#define EJEY 2
#define EJEZ 3

void genera_script_gplot3D_x_conec(cadena name_3DGnuplot_conec);
void genera_script_gplot3D_x_size(cadena name_3DGnuplot_size);
void ImgGnuplot3D_x_conectividad(RED red, cadena nom_arch, int L, double xms, double sigma);
void ImgGnuplot3D_x_size(RED red, cadena nom_arch, int L, double xms, double sigma);

#endif