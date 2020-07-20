#include <stdio.h>
#include <FPT.h>
#include <stdlib.h>
#include "M2d_matrix_toolsJeff.c"

// draw image(s)
//click n save window
//find center of mass of the window
//clip object(s) on a pershape basis

//algorithm
/*
good to good - nothing
good to bad - record intersection
bad to good - record intersection and keep destination
*/

int numpoints[10];
double swidth, sheight;
double x[10][1500], y[10][1500];
int numpolys[10];
int psize[10][1000];
int con[10][1000][20];
double red[10][1000], grn[10][1000], blu[10][1000];

double scalars[10][2];
double dxy[10][4];
double rot[10][3][3];
//double master[10][3][3];

double window[2][100];
int wn;
//int xcenter, ycenter;
//double x_inter, y_inter;

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

void create_bounding_box(double *x, double *y, int c){
	double xmin = find_min(x, numpoints[c]);
	double xmax = find_max(x, numpoints[c]);
	double ymin = find_min(y, numpoints[c]);
	double ymax = find_max(y, numpoints[c]);

	double xmid = (xmin + (xmax - xmin)/2.0);
	double ymid = (ymin + (ymax - ymin)/2.0);

	dxy[c][0] = 0 - xmid;
	dxy[c][1] = 0 - ymid;
	dxy[c][2] = (swidth/2);
	dxy[c][3] = (sheight/2);

	if((xmax-xmin) >= (ymax-ymin)){
		scalars[c][0] = (swidth-100)/(xmax-xmin);
		scalars[c][1] = scalars[c][0];
	}else{
		scalars[c][1] = (sheight-100)/(ymax-ymin);
		scalars[c][0] = scalars[c][1];
	}
}

void translate(double x[], double y[], int dx, int dy, int c){
	double a[3][3];
	M2d_make_translation(a, dx, dy);
	M2d_mat_mult_points(x, y, a, x, y, numpoints[c]);

}

void scale(double x[], double y[], int c){
	double a[3][3];
	M2d_make_scaling(a, scalars[c][0], scalars[c][1]);
	M2d_mat_mult_points(x, y, a, x, y, numpoints[c]);
}

void make_rot_mat(double x[], double y[], double angle, int c){
	double tran[3][3];
	double r[3][3];

	M2d_make_translation(tran, (-swidth/2), (-sheight/2));
	M2d_make_rotation_radians(r, angle);
	M2d_mat_mult(rot[c], r, tran);

	M2d_make_translation(tran, sheight/2, sheight/2);
	M2d_mat_mult(rot[c], tran, rot[c]);
}

int click_n_save(){
    double P[2];

    G_rgb(0, 1, 0.5);
    G_fill_rectangle(0, 0, swidth, 20);

    G_rgb(0, 0, 0);
    G_wait_click(P);
    int i , j;
		i = 0;

    while(P[1] > 20){
			//j = (i+1) % size;

        window[0][i] = P[0];//the x
        window[1][i] = P[1];//the y
        G_circle(window[0][i], window[1][i], 2);
        if(i > 0){
            //G_line(window[0][j], window[1][j], window[0][i], window[1][i]);
        }
        i++;
        G_wait_click(P);
    }

    if (i >= 1) {
        //G_line(window[0][0], window[1][0], window[0][j], window[1][j]);
				//window[0][i-1] = windor(int m = 0; m < numpolys[index]){

	}else{return -1;}
		for(int k = 0; k < 10; k++){
			//printf("%lf %lf \n", window[0][k], window[1][k]);
		}
    return i;//returns number of vertices
}

double triangle_signed_area (double A[2], double B[2], double C[2])
{
  double a ;

  a = A[0]*B[1] + B[0]*C[1] + C[0]*A[1]
    - A[1]*B[0] - B[1]*C[0] - C[1]*A[0] ;

  a = a/2 ;

  return a ;
}

double intersect_2_lines (double A[2], double B[2],
                          double C[2], double D[2],
                          double intersection[2])
{
  double areaACD, areaCBD, aquad, u ;

  // order matters
  areaACD = triangle_signed_area (A,C,D) ;
  areaCBD = triangle_signed_area (C,B,D) ;
  aquad   = areaACD + areaCBD ;

  if (aquad == 0) return 0 ;

  u = areaACD/aquad ;

  intersection[0] = A[0] + u*(B[0] - A[0]) ;
  intersection[1] = A[1] + u*(B[1] - A[1]) ;

  return 1 ;
}

int Clip_Polygon_Against_Line(
															double a, double b, double c,
															double *polyx, double *polyy, int size,
															double wx1, double wx2,
															double wy1, double wy2,
															double *resx, double *resy){

	int num, j, i;
	num = 0;
	double x1, x2, y1, y2, x21, y21, den, t, x_inter, y_inter;
	double s1, s2;
	double p1[2], p2[2], w1[2], w2[2], inter[2];


	for(i = 0; i < size; i++){
		j = (i+1) % size;

		x1 = polyx[i]; y1 = polyy[i];
		x2 = polyx[j]; y2 = polyy[j];

		s1 = (a*x1 + b*y1 + c);
		s2 = (a*x2 + b*y2 + c);

		if((s1 >= 0) && (s2 >= 0)){
			//bad to bad so do nothing
		} else if((s1 < 0) && (s2 < 0)){
			//good to good
			resx[num] = x2 ; resy[num] = y2 ; num++ ;
		} else {

			p1[0] = x1; p1[1] = y1;
			p2[0] = x2; p2[1] = y2;
			w1[0] = wx1; w1[1] = wy1;
			w2[0] = wx2; w2[1] = wy2;
			intersect_2_lines(p1,p2,w1,w2,inter);




				/*x21 = x2 - x1 ; y21 = y2 - y1 ;
        den = a*x21 + b*y21 ;
        if (den == 0) continue ; // do nothing-should never happen
        t = -(a*x1 + b*y1 + c)/den ;*/
				//printf("t = %lf, den = %lf\n", t, den);
        //x_inter = x1 + t*x21 ;
        //y_inter = y1 + t*y21 ;
				x_inter = inter[0];
				y_inter = inter[1];

			if(s1 < 0){
				resx[num] = x_inter; resy[num] = y_inter; num++;
			} else {// end for k
				resx[num] = x_inter; resy[num] = y_inter; num++;
				resx[num] = x2; resy[num] = y2; num++;
			}
		}
	}
	//returns size of new poly
	return num;
}

int  Clip_Polygon_Against_Convex_Window (
      double *px,  double *py, int psize,
      double *wx,  double *wy, int wsize)

{
   double nx[100],ny[100],  a,b,c,  cwx,cwy ;
   int i,k,m ;
   // find center of mass of window
   cwx = 0.0 ; cwy = 0.0 ;
   for (k = 0 ; k < wsize ; k++) {
     cwx += wx[k] ; cwy += wy[k] ;
   }
   cwx /= wsize ; cwy /= wsize ;

   // clip the polygon against each edge of the window
   for (k = 0 ; k < wsize ; k++) {
      m = k+1 ; if (m == wsize) { m = 0 ; }

      // ax + by + c = 0 is eqn of this window edge
			a = wy[m] - wy[k] ;
      b = wx[k] - wx[m] ;
      c = -(a*wx[k] + b*wy[k]);

      // but we need for ax + by + c < 0 to reflect "inside"
      if (a*cwx + b*cwy + c > 0) {
				a = -a ; b = -b ; c = -c ;
      }

      psize = Clip_Polygon_Against_Line (a,b,c,px,py,psize,wx[k],wx[m],wy[k],wy[m],nx,ny) ;


     // copy back in preparation for next pass
     for (i = 0 ; i < psize ; i++) {
			 px[i] = nx[i]; py[i] = ny[i];
     }
     //G_fill_polygon(px,py,psize) ;

   }
   return psize ;
}


/*int make_matrix(int i){
	double tran1[3][3]; double tran2[3][3];
	double scale[3][3]; double rot[3][3];
	double angle = 0.03;

	M2d_make_translation(tran1, dxy[i][0], dxy[i][1]);
	M2d_make_translation(tran2, dxy[i][2], dxy[i][3]);
	M2d_make_scaling(scale, scalars[i][0], scalars[i][1]);
	M2d_make_rotation_radians(rot, angle);

	M2d_mat_mult(master[i], tran1, scale);
	M2d_mat_mult(master[i], master[i], rot);
	M2d_mat_mult(master[i], master[i], tran2);

	create_bounding_box(x[i], y[i], c);

	translate(x[i], y[i], dxy[i][0], dxy[i][1], i);
	scale(x[i], y[i], i);
	translate(x[i], y[i], dxy[i][2], dxy[i][3], i);

	make_rot_mat(x[i], y[i], 0.03, i);
	return 1;
}*/

void draw_first(int c){
	for(int i = 0; i < numpolys[c]; i++){
		double tempx[psize[c][i]];
		double tempy[psize[c][i]];
		for(int j = 0; j < psize[c][i]; j++){
			int index = con[c][i][j];
			tempx[j] = x[c][index];
			tempy[j] = y[c][index];
		}
		G_rgb(red[c][i], grn[c][i], blu[c][i]);
		G_fill_polygon(tempx, tempy, psize[c][i]);
	}
}

//*x is the x's for a given poly
//*y is the y's for a given poly
//n is the size of those two arrays
void draw(double *x, double *y, int n){
		int num;
		double tempx[n]; double tempy[n];
		for(int k = 0; k < n; k++){
			tempx[k] = x[k];
			tempy[k] = y[k];
		}
		num = Clip_Polygon_Against_Convex_Window(tempx, tempy, n, window[0], window[1], wn) ;
		G_fill_polygon(tempx, tempy, num);
}

int main(int argc, char **argv){

	swidth = 750; sheight = 750;
	G_init_graphics(swidth, sheight);

	FILE *fp;
	int numobjects = argc-1;
	for(int c = 0; c < numobjects; c++){
		//gets all of the stuff
		//-----------------------------------------------------
		fp = fopen(argv[c+1], "r");
		if(!fp){
			printf("please choose another input \n");
			exit(0);
		}
		fscanf(fp, "%d", &numpoints[c]);

		for(int i = 0; i < numpoints[c]; i++){
			fscanf(fp, "%lf %lf", &x[c][i], &y[c][i]);
		}

		fscanf(fp, "%d", &numpolys[c]);

		for(int i = 0; i < numpolys[c]; i++){
			fscanf(fp, "%d", &psize[c][i]);
			for(int j = 0; j < psize[c][i]; j++){
				fscanf(fp, "%d", &con[c][i][j]);
			}
		}
		for(int i = 0; i < numpoints[c]; i++){
			fscanf(fp, "%lf %lf %lf", &red[c][i], &grn[c][i], &blu[c][i]);
		}
		//--------------------------------------------------------
		create_bounding_box(x[c], y[c], c);

		translate(x[c], y[c], dxy[c][0], dxy[c][1], c);
		scale(x[c], y[c], c);
		translate(x[c], y[c], dxy[c][2], dxy[c][3], c);

		make_rot_mat(x[c], y[c], 0.03, c);
		//make_matrix(c);
	}

	int i, key, done, index;
	key = 0;
	done = 1;
	key = G_wait_key();
	index = key-48;
	//first draw
  draw_first(index);
  wn = click_n_save();

	//i is the object
	//j is the polygon
	//k is the connection
	while(done == 1){
		if(key == 'q'){exit(0);}
		i = key - 48;
		G_clear();
		for(int j = 0; j < numpolys[i]; j++){
			G_rgb(0,0,0);
			//G_clear();
			double thisx[psize[i][j]]; double thisy[psize[i][j]];
			for(int k = 0; k < psize[i][j]; k++){
				thisx[k] = x[i][con[i][j][k]]; thisy[k] = y[i][con[i][j][k]];
			}
			G_rgb(red[i][j], grn[i][j], blu[i][j]);
			draw(thisx, thisy, psize[i][j]);
		}
		M2d_mat_mult_points(x[i], y[i], rot[i], x[i], y[i], numpoints[i]);
		G_polygon(window[0], window[1], wn);

		key = G_wait_key();
		G_rgb(1, 1, 1);
		G_clear();
	}
	G_close();
}
