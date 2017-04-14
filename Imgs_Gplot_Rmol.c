#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "pred_MC.h"
#include "Imgs_Gplot_Rmol.h"


/*void genera_script_gplot2D(cadena name_corte2DGnuplot){
	FILE *ap;

	ap=fopen("script2D","w");
	if(ap == NULL)
		printf("ERROR.... no se pudo crear archivo script2D\n");
	else{
		fprintf(ap,"set hidden3d\n");
		fprintf(ap,"set output \"%s.eps\"\n",name_corte2DGnuplot);
		fprintf(ap,"set terminal postscript eps enhanced color\n");
		fprintf(ap,"plot \"Big_%s.xy\" u 2:3 with points title \"big\",",name_corte2DGnuplot);
		fprintf(ap," \"Medium_%s.xy\" u 2:3 with points title \"medium\",",name_corte2DGnuplot);
		fprintf(ap," \"Small_%s.xy\" u 2:3 with points title \"small\"\n\n",name_corte2DGnuplot);
		fclose(ap);
		printf("\nScript2D generado\n");
	}
}*/

void genera_script_gplot3D_x_conec(cadena name_corte3DGnuplot)
{FILE *ap;
    
    ap=fopen("script3Dconec","w");
    if(ap == NULL)
        printf("ERROR.... no se pudo crear archivo script3D\n");
    else {
        fprintf(ap,"set hidden3d\n");
        fprintf(ap,"set output \"%s.eps\"\n",name_corte3DGnuplot);
        fprintf(ap,"set terminal postscript eps enhanced color size 4,4.5 font 16\n");
        fprintf(ap,"splot \"Big_%s.xyz\" u 2:3:4 with points title \"5-6 bonds\",",name_corte3DGnuplot);
        fprintf(ap," \"Medium_%s.xyz\" u 2:3:4 with points title \"3-4 bonds\",",name_corte3DGnuplot);
        fprintf(ap," \"Small_%s.xyz\" u 2:3:4 with points title \"1-2 bonds\"\n\n",name_corte3DGnuplot);
        fclose(ap);
        printf("\nScript3D generado\n");
        system("gnuplot script3Dconec");
    }
}
void genera_script_gplot3D_x_size(cadena name_corte3DGnuplot)
{FILE *ap;
    
    ap=fopen("script3Dsize","w");
    if(ap == NULL)
        printf("ERROR.... no se pudo crear archivo script3D\n");
    else {
        fprintf(ap,"set hidden3d\n");
        fprintf(ap,"set output \"%s.eps\"\n",name_corte3DGnuplot);
        fprintf(ap,"set terminal postscript eps enhanced color size 4,4.5 font 16\n");
        fprintf(ap,"splot \"Big_%s.xyz\" u 2:3:4 with points title \"big\",",name_corte3DGnuplot);
        fprintf(ap," \"Medium_%s.xyz\" u 2:3:4 with points title \"medium\",",name_corte3DGnuplot);
        fprintf(ap," \"Small_%s.xyz\" u 2:3:4 with points title \"small\"\n\n",name_corte3DGnuplot);
        fclose(ap);
        printf("\nScript3D generado\n");
        system("gnuplot script3Dsize");
    }
}

/*void genera_script_gplot3D(cadena name_corte3DGnuplot)
{FILE *ap;

  ap=fopen("script3D","w");
  if(ap == NULL)
    printf("ERROR.... no se pudo crear archivo script3D\n");
  else {
    fprintf(ap,"set hidden3d\n");
    fprintf(ap,"set output \"%s.eps\"\n",name_corte3DGnuplot);
    fprintf(ap,"set terminal postscript eps enhanced color size 4,4.5 font 16\n");
    fprintf(ap,"splot \"Big_%s.xyz\" u 2:3:4 with points title \"big\",",name_corte3DGnuplot);
    fprintf(ap," \"Medium_%s.xyz\" u 2:3:4 with points title \"medium\",",name_corte3DGnuplot);
    fprintf(ap," \"Small_%s.xyz\" u 2:3:4 with points title \"small\"\n\n",name_corte3DGnuplot);
    fclose(ap);
    printf("\nScript3D generado\n");
  }
}*/

/*void min_max_corte(RED red,double *min,double *max,int Z,int L,int eje)
{int i,j;
 double mi,ma;
   mi=red[0][0][0].sitio.radio;
   ma=mi;
   switch(eje)
   {  case EJEZ:{
   for(i=0;i<L;i++)
      for(j=0;j<L;j++)
      {
        if(red[i][j][Z].sitio.radio < mi)
          mi=red[i][j][Z].sitio.radio;
        else
        if(red[i][j][Z].sitio.radio > ma)
          ma=red[i][j][Z].sitio.radio;
      }
                 break;}
      case EJEX:{
   for(i=0;i<L;i++)
      for(j=0;j<L;j++)
      {
        if(red[i][Z][j].sitio.radio < mi)
          mi=red[i][Z][j].sitio.radio;
        else
        if(red[i][Z][j].sitio.radio > ma)
          ma=red[i][Z][j].sitio.radio;
      }
               break;}
       case EJEY:{
   for(i=0;i<L;i++)
      for(j=0;j<L;j++)
      {
        if(red[Z][i][j].sitio.radio < mi)
          mi=red[Z][i][j].sitio.radio;
        else
        if(red[Z][i][j].sitio.radio > ma)
          ma=red[Z][i][j].sitio.radio;
      }
                break; }
   }
   *min=mi; *max=ma;

}*/
/*
El formato para imágenes sería con la siguiente clasificación de sitios:
tipo    intervalo
S       sitio<=(xms-0.43*sigma)
H      (xms-0.43*sigma)<sitio<=(xms+0.43*sigma)
N      (xms+0.43*sigma)<sitio
las letras obedecen a que en la lógica del Rasmol S,H y N son tipos de
átomos.
Ahora bién, el formato del archivo texto sería:

tipo x y z

Por último, el nombre del archivo tiene que ir conla extensión "xyz"
*/

/*int hayEnlace(double renlace)
{
    if(renlace > 0.0)
        return 1;
    else
        return 0;
}*/


int cuentaEnlaces(int i,int j,int k,int L,RED red)
{
    int enlaces=0;
 /*   
    enlaces = enlaces + hayEnlace(val_enlace(i,j,k,DERECHA,red,L));
    enlaces = enlaces + hayEnlace(val_enlace(i,j,k,IZQUIERDA,red,L));
    enlaces = enlaces + hayEnlace(val_enlace(i,j,k,ABAJO,red,L));
    enlaces = enlaces + hayEnlace(val_enlace(i,j,k,ARRIBA,red,L));
    enlaces = enlaces + hayEnlace(val_enlace(i,j,k,FONDO,red,L));
    enlaces = enlaces + hayEnlace(val_enlace(i,j,k,FRENTE,red,L));
 */   
    if(val_enlace(i,j,k,DERECHA,red,L)>0.0) enlaces++;
    if(val_enlace(i,j,k,IZQUIERDA,red,L)>0.0) enlaces++;
    if(val_enlace(i,j,k,ABAJO,red,L)>0.0) enlaces++;
    if(val_enlace(i,j,k,ARRIBA,red,L)>0.0) enlaces++;
    if(val_enlace(i,j,k,FONDO,red,L)>0.0) enlaces++;
    if(val_enlace(i,j,k,FRENTE,red,L)>0.0) enlaces++;

    return enlaces;
}

/*void ImgGnuplot3D_x_2conectividad(RED red,cadena nom_arch,int L,double xms, double sigma)
{
    int i,j,k,f;
    char c;
    FILE *ap_arch_small;
    FILE *ap_arch_medium;
    FILE *ap_arch_big;
    cadena nom_arch_small;
    cadena nom_arch_medium;
    cadena nom_arch_big;
    

  /*  sprintf(nom_arch_small,"%s_Small.xyz",nom_arch);
    sprintf(nom_arch_medium,"%s_Medium.xyz",nom_arch);
    sprintf(nom_arch_big,"%s_Big.xyz",nom_arch);
 /
    sprintf(nom_arch_small,"Small_%s.xyz",nom_arch);
    sprintf(nom_arch_medium,"Medium_%s.xyz",nom_arch);
    sprintf(nom_arch_big,"Big_%s.xyz",nom_arch);

    ap_arch_small=fopen(nom_arch_small,"w");
    ap_arch_medium=fopen(nom_arch_medium,"w");
    ap_arch_big=fopen(nom_arch_big,"w");
    
    fprintf(ap_arch_small,"tipo x y z\n");
    fprintf(ap_arch_medium,"tipo x y z\n");
    fprintf(ap_arch_big,"tipo x y z\n");
    
    for(i=0;i<L;i++)
        for(j=0;j<L;j++)
            for(k=0;k<L;k++){
                f=cuentaEnlaces(i,j,k,L,red);
                if( f > 0 ){
                    if( f<=2 ){
                       fprintf(ap_arch_small,"S %d %d %d\n",i,j,k);// 1 o 2 enlaces
                    }
                    else{
                        if( f<=4 ){
                           fprintf(ap_arch_medium,"H %d %d %d\n",i,j,k);// 3 o 4 enlaces
                        }
                       else {
                           if( f<=6 ){
                        	fprintf(ap_arch_big,"N %d %d %d\n",i,j,k);// 5 o 6 enlaces
                           }
                       }
                    }
                }
            }
    fclose(ap_arch_small);
    fclose(ap_arch_medium);
    fclose(ap_arch_big);
    //printf("\nArchivos 3D Creados para GnuPlot:\nArchivo1:%s\nArchivo2:%s\nArchivo3:%s\n",nom_arch_small,nom_arch_medium,nom_arch_big);
}*/

void ImgGnuplot3D_x_conectividad(RED red,cadena nom_arch,int L,double xms, double sigma)
{  int i,j,k,f,neG=0,neM=0,neS=0,enl_vacios=0,nviol_en_sit;
    FILE *ap_arch_small;
    FILE *ap_arch_medium;
    FILE *ap_arch_big;
    cadena nom_arch_small;
    cadena nom_arch_medium;
    cadena nom_arch_big;
    
    sprintf(nom_arch_small,"Small_%s.xyz",nom_arch);
    sprintf(nom_arch_medium,"Medium_%s.xyz",nom_arch);
    sprintf(nom_arch_big,"Big_%s.xyz",nom_arch);
    
    ap_arch_small=fopen(nom_arch_small,"w");
    ap_arch_medium=fopen(nom_arch_medium,"w");
    ap_arch_big=fopen(nom_arch_big,"w");
    
    fprintf(ap_arch_small,"tipo x y z\n");
    fprintf(ap_arch_medium,"tipo x y z\n");
    fprintf(ap_arch_big,"tipo x y z\n");

printf("L= %d\n",L);
    for(i=0;i<L;i++)
        for(j=0;j<L;j++)
            for(k=0;k<L;k++)
            {
                f=cuentaEnlaces(i,j,k,L,red);
                if(f>0){
                  //printf("ne:%d ",f);
                   if(f <= 2 ){
                       neS++;
                       fprintf(ap_arch_small,"S %d %d %d\n",i,j,k);
						}
                   else
                       if(f <= 4){
									neM++;
                           fprintf(ap_arch_medium,"H %d %d %d\n",i,j,k);
								}
                       else{
								neG++;
                        fprintf(ap_arch_big,"N %d %d %d\n",i,j,k);
							 }
                }
                //else printf("cero enlaces en %d,%d,%d ",i,j,k);
            }
  //  printf("neG: %d ,neM:%d, neS:%d\n",neG,neM,neS);
    fclose(ap_arch_small);
    fclose(ap_arch_medium);
    fclose(ap_arch_big);
    printf("\nArchivos 3D Creados para GnuPlot:\nArchivo1:%s\nArchivo2:%s\nArchivo3:%s\n",nom_arch_small,nom_arch_medium,nom_arch_big);
}

void ImgGnuplot3D_x_size(RED red,cadena nom_arch,int L,double xms, double sigma)
{  int i,j,k,f;
    FILE *ap_arch_small;
    FILE *ap_arch_medium;
    FILE *ap_arch_big;
    cadena nom_arch_small;
    cadena nom_arch_medium;
    cadena nom_arch_big;
    
    sprintf(nom_arch_small,"Small_%s.xyz",nom_arch);
    sprintf(nom_arch_medium,"Medium_%s.xyz",nom_arch);
    sprintf(nom_arch_big,"Big_%s.xyz",nom_arch);
    
    ap_arch_small=fopen(nom_arch_small,"w");
    ap_arch_medium=fopen(nom_arch_medium,"w");
    ap_arch_big=fopen(nom_arch_big,"w");
    
    fprintf(ap_arch_small,"tipo x y z\n");
    fprintf(ap_arch_medium,"tipo x y z\n");
    fprintf(ap_arch_big,"tipo x y z\n");
    
    for(i=0;i<L;i++)
        for(j=0;j<L;j++)
            for(k=0;k<L;k++)
            {
                f=cuentaEnlaces(i,j,k,L,red);
                if(f>0){
                  //printf("ne:%d ",f);
                if(red[i][j][k].sitio.radio <=(xms-0.43*sigma) )
                    fprintf(ap_arch_small,"S %d %d %d\n",i,j,k);
                else
                    if(((xms-0.43*sigma) < red[i][j][k].sitio.radio ) && (red[i][j][k].sitio.radio <= (xms+0.43*sigma) ))
                        fprintf(ap_arch_medium,"H %d %d %d\n",i,j,k);
                    else
                        fprintf(ap_arch_big,"N %d %d %d\n",i,j,k);
                }
                //else printf("cero enlaces en %d,%d,%d ",i,j,k);
            }
    printf("\n");
    fclose(ap_arch_small);
    fclose(ap_arch_medium);
    fclose(ap_arch_big);
    printf("\nArchivos 3D Creados para GnuPlot:\nArchivo1:%s\nArchivo2:%s\nArchivo3:%s\n",nom_arch_small,nom_arch_medium,nom_arch_big);
}


/*void ImgGnuplot3D(RED red,cadena nom_arch,int L,double xms, double sigma)
{  int i,j,k;
    FILE *ap_arch_small;
    FILE *ap_arch_medium;
    FILE *ap_arch_big;
    cadena nom_arch_small;
    cadena nom_arch_medium;
    cadena nom_arch_big;
    
    sprintf(nom_arch_small,"Small_%s.xyz",nom_arch);
    sprintf(nom_arch_medium,"Medium_%s.xyz",nom_arch);
    sprintf(nom_arch_big,"Big_%s.xyz",nom_arch);
    
    ap_arch_small=fopen(nom_arch_small,"w");
    ap_arch_medium=fopen(nom_arch_medium,"w");
    ap_arch_big=fopen(nom_arch_big,"w");
    
    fprintf(ap_arch_small,"tipo x y z\n");
    fprintf(ap_arch_medium,"tipo x y z\n");
    fprintf(ap_arch_big,"tipo x y z\n");
    
    for(i=0;i<L;i++)
        for(j=0;j<L;j++)
            for(k=0;k<L;k++)
            {
                if(red[i][j][k].sitio.radio <=(xms-0.43*sigma) )
                    fprintf(ap_arch_small,"S %d %d %d\n",i,j,k);
                else
                    if(((xms-0.43*sigma) < red[i][j][k].sitio.radio ) && (red[i][j][k].sitio.radio <= (xms+0.43*sigma) ))
                        fprintf(ap_arch_medium,"H %d %d %d\n",i,j,k);
                    else
                        fprintf(ap_arch_big,"N %d %d %d\n",i,j,k);
            }
    fclose(ap_arch_small);
    fclose(ap_arch_medium);
    fclose(ap_arch_big);
    printf("\nArchivos 3D Creados para GnuPlot:\nArchivo1:%s\nArchivo2:%s\nArchivo3:%s\n",nom_arch_small,nom_arch_medium,nom_arch_big);
}*/

/*void ImgRasmol3D(RED red,cadena nom_arch,int L,double xms, double sigma)
{  int i,j,k;
   char c;
   FILE *ap_arch;

   ap_arch=fopen(nom_arch,"w");
  
   fprintf(ap_arch,"tipo x y z\n");
   for(i=0;i<L;i++)
    for(j=0;j<L;j++)
      for(k=0;k<L;k++)
      {
        if(red[i][j][k].sitio.radio <=(xms-0.43*sigma) )
             c='S';
        else
           if(((xms-0.43*sigma) < red[i][j][k].sitio.radio ) && (red[i][j][k].sitio.radio <= (xms+0.43*sigma) ))
             c='H';
           else
               c= 'N';
        fprintf(ap_arch,"%c %d %d %d\n",c,i,j,k);
     }
  fclose(ap_arch);
  printf("\nArchivo 3D Creado para Rasmol:\nArchivo:%s\n",nom_arch);
}*/

/*
El formato para imágenes sería con la siguiente clasificación de sitios:
tipo    intervalo
S       sitio<=(xms-0.43*sigma)
H      (xms-0.43*sigma)<sitio<=(xms+0.43*sigma)
N      (xms+0.43*sigma)<sitio
las letras obedecen a que en la lógica del Rasmol S,H y N son tipos de
átomos.
Ahora bién, el formato del archivo texto sería:

tipo x y z

Por ultimo, el nombre del archivo tiene que ir con la extensión "xyz"
*/
/*void ImgRasmol2D(RED red,cadena name_corte,int Z,int L,double xms, double sigma,int eje)
{  int i,j;
   char c;
   FILE *ap_arch;


   ap_arch=fopen(name_corte,"w");
   fprintf(ap_arch,"tipo x y\n");
   switch(eje)
   {  case EJEY: {
   for(i=0;i<L;i++)
      for(j=0;j<L;j++)
      {
        if(red[Z][i][j].sitio.radio >  (xms+0.43*sigma) )
             c= 'N';
        else 
           if(((xms-0.43*sigma) < red[Z][i][j].sitio.radio ) && (red[Z][i][j].sitio.radio <= (xms+0.43*sigma) ))
             c='H';
           else
               c= 'S';
        fprintf(ap_arch,"%c %d %d\n",c,i,j);
     }
                      break;}
       case EJEX: {
   for(i=0;i<L;i++)
      for(j=0;j<L;j++)
      {
        if(red[i][Z][j].sitio.radio >  (xms+0.43*sigma) )
             c= 'N';
        else
           if(((xms-0.43*sigma) < red[i][Z][j].sitio.radio ) && (red[i][Z][j].sitio.radio <= (xms+0.43*sigma) ))
             c='H';
           else
               c= 'S';
          fprintf(ap_arch,"%c %d %d\n",c,i,j);
     }
                      break;}
        case EJEZ: {
   for(i=0;i<L;i++)
      for(j=0;j<L;j++)
      {
        if(red[i][j][Z].sitio.radio >  (xms+0.43*sigma) )
             c= 'N';
        else
           if(((xms-0.43*sigma) < red[i][j][Z].sitio.radio ) && (red[i][j][Z].sitio.radio <= (xms+0.43*sigma) ))
             c='H';
           else
               c= 'S';
          fprintf(ap_arch,"%c %d %d\n",c,i,j);
     }
                   break;}
   }
  fprintf(ap_arch,"\n");
  fclose(ap_arch);
  printf("\nArchivo 2D Creado para Rasmol:\nArchivo:%s, corte en:%d\n",name_corte,Z);
}*/

/*void ImgGnuplot2D(RED red,cadena nom_arch,int Z,int L,double xms, double sigma,int eje)
{  int i,j;

   FILE *ap_arch_small;
   FILE *ap_arch_medium;
   FILE *ap_arch_big;
   cadena nom_arch_small;
   cadena nom_arch_medium;
   cadena nom_arch_big;


   sprintf(nom_arch_small,"Small_%s.xy",nom_arch);
   sprintf(nom_arch_medium,"Medium_%s.xy",nom_arch);
   sprintf(nom_arch_big,"Big_%s.xy",nom_arch);

   ap_arch_small=fopen(nom_arch_small,"w");
   if(ap_arch_small == NULL)
      printf("Error, apertura archivo %s\n",nom_arch_small);
   ap_arch_medium=fopen(nom_arch_medium,"w");
   if(ap_arch_medium == NULL)
      printf("Error, apertura archivo %s\n",nom_arch_medium);
   ap_arch_big=fopen(nom_arch_big,"w");
   if(ap_arch_big == NULL)
      printf("Error, apertura archivo %s\n",nom_arch_big);
   fprintf(ap_arch_small,"tipo x y\n");

   fprintf(ap_arch_medium,"tipo x y\n");

   fprintf(ap_arch_big,"tipo x y\n");

   switch(eje)
   {  case EJEY: {
   for(i=0;i<L;i++)
      for(j=0;j<L;j++)
      {
        if(red[Z][i][j].sitio.radio >  (xms+0.43*sigma) )
            fprintf(ap_arch_big,"N %d %d\n",i,j);
        else 
           if(((xms-0.43*sigma) < red[Z][i][j].sitio.radio ) && (red[Z][i][j].sitio.radio <= (xms+0.43*sigma) ))
                fprintf(ap_arch_medium,"H %d %d\n",i,j);
           else
        	fprintf(ap_arch_small,"S %d %d\n",i,j);
     }
                      break;}
       case EJEX: {
   for(i=0;i<L;i++)
      for(j=0;j<L;j++)
      {
        if(red[i][Z][j].sitio.radio >  (xms+0.43*sigma) )
            fprintf(ap_arch_big,"N %d %d\n",i,j);
        else
           if(((xms-0.43*sigma) < red[i][Z][j].sitio.radio ) && (red[i][Z][j].sitio.radio <= (xms+0.43*sigma) ))
                fprintf(ap_arch_medium,"H %d %d\n",i,j);
           else
        	fprintf(ap_arch_small,"S %d %d\n",i,j);
     }
                      break;}
        case EJEZ: {
   for(i=0;i<L;i++)
      for(j=0;j<L;j++)
      {
        if(red[i][j][Z].sitio.radio >  (xms+0.43*sigma) )
            fprintf(ap_arch_big,"N %d %d\n",i,j);
        else
           if(((xms-0.43*sigma) < red[i][j][Z].sitio.radio ) && (red[i][j][Z].sitio.radio <= (xms+0.43*sigma) ))
                fprintf(ap_arch_medium,"H %d %d\n",i,j);
           else
        	fprintf(ap_arch_small,"S %d %d\n",i,j);
     }
                   break;}
   }

  fprintf(ap_arch_small,"\n");
  fprintf(ap_arch_medium,"\n");
  fprintf(ap_arch_big,"\n");
  fclose(ap_arch_small);
  fclose(ap_arch_medium);
  fclose(ap_arch_big);
  printf("\nArchivos 2D Creados para GnuPlot:\nArchivo1:%s\nArchivo2:%s\nArchivo3:%s\n",nom_arch_small,nom_arch_medium,nom_arch_big);
}*/



