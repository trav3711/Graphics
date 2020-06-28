#include <FPT.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "M2d_matrix_tools.c"

int numpoints;
double swidth, sheight;
double x[1500], y[1500];
int numpolys;
int psize[1000];
int con[1000][20];
double red[1000], grn[1000], blu[1000];
double a[3][3];

double find_min(double a[], int n){
	int min = a[0];
	for(int i = 1; i < n; i++){
		if(a[i] < min){min = a[i];}
	}
	return min;
}

double find_max(double a[], int n){
	int max = a[0];
	for(int i = 0; i < n; i++){
		if(max < a[i]){max = a[i];}
	}
	return max;
}


void translate(double x[], double y[], int newx, int newy){
	double xmin = find_min(x, numpoints);
	double ymin = find_min(y, numpoints);
	double xmax = find_max(x, numpoints);
	double ymax = find_max(y, numpoints);
	double dx = newx - (xmin + (xmax-xmin)/2.0);
	double dy = newy - (ymin + (ymax-ymin)/2.0);
	
	M2d_make_translation(a, dx, dy);
	M2d_mat_mult_points(x, y, a, x, y, numpoints);
}

void scale(double x[], double y[]){
	double xmin = find_min(x, numpoints);
	double xmax = find_max(x, numpoints);
	double ymin = find_min(y, numpoints);
	double ymax = find_max(y, numpoints);
	double sx, sy;

	if((xmax-xmin) >= (ymax-ymin)){
		sx = (swidth-100)/(xmax-xmin);
		sy = sx;
	}else{
		sy = (sheight-100)/(ymax-ymin);
		sx = sy;
	}
	M2d_make_scaling(a, sx, sy);
	M2d_mat_mult_points(x, y, a, x, y, numpoints);
}

void rotate(double x[], double y[], double angle){
	translate(x, y, 0, 0);
	M2d_make_rotation_radians(a, angle);
	M2d_mat_mult_points(x, y, a, x, y, numpoints);
	translate(x, y, sheight/2, swidth/2);
}

void draw(){
	for(int i = 0; i < numpolys; i++){
		double tempx[psize[i]];
		double tempy[psize[i]];
		for(int j = 0; j < psize[i]; j++){
			int index = con[i][j];
			tempx[j] = x[index];
			tempy[j] = y[index];
			
		}
		G_rgb(red[i], grn[i], blu[i]);
		G_fill_polygon(tempx, tempy, psize[i]);
	}
}

int main(int argc, char **argv){

	swidth = 750; sheight = 750;
	G_init_graphics(swidth, sheight);

	FILE *fp;
	int numobjects = argc-1;

	int key1, key2 = 0;

	while(key1 != 'q'){
		key1 = G_wait_key();//here

		fp = fopen(argv[key1-48], "r"); 

		if(!fp && key1 != 'q'){
			printf("please choose another input \n");
			exit(0);
		}
		fscanf(fp, "%d", &numpoints);
	
		for(int i = 0; i <numpoints; i++){
			fscanf(fp, "%lf %lf", &x[i], &y[i]);
		}
		
		fscanf(fp, "%d", &numpolys);

		for(int i = 0; i < numpolys; i++){
			fscanf(fp, "%d", &psize[i]);
			for(int j = 0; j < psize[i]; j++){
				fscanf(fp, "%d", &con[i][j]);
			}
		}
		for(int i = 0; i < numpoints; i++){
			fscanf(fp, "%lf %lf %lf", &red[i], &grn[i], &blu[i]);
		}

		translate(x, y, 0, 0);
		scale(x, y);
		translate(x, y, swidth/2, sheight/2);
		draw();

		do {
			if(key2 == 'q'){exit(0);}

			G_rgb(1, 1, 1);
			G_clear();
			
			rotate(x, y, 0.03);
			
			draw();
			key2 = G_wait_key();
		}while(key1 == key2);
	}

	G_close();
}
