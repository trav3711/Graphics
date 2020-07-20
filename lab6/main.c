//0, 1, 2, ... select object
//t - translate
//r - x,y,z - rotate
//c - change sign
//h = tan(theta/2)

#include <FPT.h>
#include "M3d_matrix_tools.c"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double swidth, sheight;

int numobjects;
int numpolys[10], numpoints[10];
double x[10][15000], y[10][15000], z[10][15000];
//double drawx[10][15000], drawy[10][15000], drawz[10][15000];
double red[10], grn[10], blu[10];
double xbar[10][15000], ybar[10][15000];
int psize[10][10000];
int con[10][15000][10];
int viewangle = 40;

double xaverage[10], yaverage[10], zaverage[10];

int sign = 1;
int key;
int action = 't';
int onum = 0;
double rotation_angle = .05;

double xmin, ymin, zmin, xmax, ymax, zmax;
double xcenter, ycenter, zcenter;

int backface;
int s = 1;

typedef
struct {
  int objnum;
  int polynum;
  double dist;
}
POLY;

POLY zdist[150000];

int printa(double *a, int n){
  for(int i = 0; i < n; i++){
    //if(a[i] == inf)
    printf("%d - %lf\n",i,a[i]);
  }
  printf("\n");
}

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

int average(int n){
	double xsum = 0; double ysum = 0; double zsum = 0;
	for(int i = 0; i < numpoints[n]; i++){
		xsum += x[n][i];
		ysum += y[n][i];
		zsum += z[n][i];
	}
	xaverage[n] = xsum/numpoints[n];
	yaverage[n] = ysum/numpoints[n];
	zaverage[n] = zsum/numpoints[n];
}

int translate_scale(){
	double temp[4][4];
	double temp2[4][4];

	//centers and scales
	for(int i = 0; i < numobjects; i++){
		M3d_make_identity(temp);
		M3d_make_identity(temp2);
    double sf;
    int n = numpoints[i];

		average(i);
		M3d_make_translation(temp, -xaverage[i], -yaverage[i], -zaverage[i]);

    xmax = find_max(x[i], n);
    ymax = find_max(y[i], n);
    zmax = find_max(z[i], n);
    xmin = find_min(x[i], n);
    ymin = find_min(y[i], n);
    zmin = find_min(z[i], n);

		double xdifference = xmax - xmin;
    double ydifference = ymax - ymin;
    double zdifference = zmax - zmin;
	  if ((xdifference > ydifference) && (xdifference > zdifference)) {
			sf = xdifference;
	  }
	  else if ((ydifference > xdifference) && (ydifference > zdifference)) {
      sf = ydifference;
	  }
	  else {
      sf = zdifference;
	  }
    if(sf == 0){
      sf = 2;
    }
    M3d_make_scaling(temp2, 1/sf, 1/sf, 1/sf);
		M3d_mat_mult(temp2, temp, temp2);
		M3d_make_translation(temp, xaverage[i], yaverage[i], zaverage[i]+2);
		M3d_mat_mult(temp2, temp, temp2);
		M3d_mat_mult_points(x[i], y[i], z[i], temp2, x[i], y[i], z[i], n);

    //M3d_make_translation(temp2, 0, )
    average(i);
  }
}

int create_2d_lists(){
  double x1, y1;

	for(int i = 0; i < numobjects; i++){
    //int nsize = backface(i);

		for(int j = 0; j < numpoints[i]; j++){


      if(z[i][j] == 0){
        x1 = z[i][j]; y1 = y[i][j];
      } else {
        x1 = x[i][j]/z[i][j];
        y1 = y[i][j]/z[i][j];
      }
      xbar[i][j] = (400/tan(viewangle*(M_PI/180)))*(x1)+400;
  		ybar[i][j] = (400/tan(viewangle*(M_PI/180)))*(y1)+400;
		}
	}
}

int backface_elimination(double *px, double *py, double *pz, int size){
  double dot;
  double xcenter, ycenter, zcenter;
  double Ax, Ay, Az, Bx, By, Bz;
  double Nx, Ny, Nz, Ex, Ey, Ez;

  xcenter = ycenter = zcenter = 0 ;
  for(int i = 0; i < size; i++){
    xcenter += px[i];
    ycenter += py[i];
    zcenter += pz[i];
  }
  xcenter/=size; ycenter/=size; zcenter/=size;

  Ax = px[0] - xcenter; Bx = px[1] - xcenter;
  Ay = py[0] - ycenter; By = py[1] - ycenter;
  Az = pz[0] - zcenter; Bz = pz[1] - zcenter;

  Nx = (Ay*Bz) - (By*Az); Ny = (Bx*Az) - (Ax*Bz); Nz = (Ax*By) - (Bx*Ay);
  Ex = -xcenter; Ey = -ycenter; Ez = -zcenter;

  dot = (Nx*Ex) + (Ny*Ey) + (Nz*Ez);
  //printf("dot = %lf\n", dot);

  if(s == 1){
    if(dot >= 0){
      return 1;
    } else {return 0;}
  } else {
    if(dot < 0){
      return 1;
    } else {return 0;}
  }
}

int read_objects(int argc, char **argv){
  numobjects = argc-1;
  FILE *fp;
  for(int i = 0; i < numobjects; i++){
    fp = fopen(argv[i+1], "r");
    if(!fp) {exit(0);}
    fscanf(fp, "%d", &numpoints[i]);

    for(int j = 0; j < numpoints[i]; j++){
      fscanf(fp, "%lf %lf %lf", &x[i][j], &y[i][j], &z[i][j]);
      //printf("z at %d = %lf\n",i, z[i][j]);
    }
    fscanf(fp, "%d", &numpolys[i]);
		//printf("numpolys at %d = %d\n", i, numpolys[i]);
    for(int j = 0; j < numpolys[i]; j++){
      fscanf(fp, "%d", &psize[i][j]);
      for(int k = 0; k < psize[i][j]; k++){
        fscanf(fp, "%d", &con[i][j][k]);
      }
    }
  }
	return 1;
}

int draw_wire_frame(int onum){
  backface = 1;
  double px[20],py[20],pz[20] ;

	//printf("numpolys = %d\n", numpolys[onum]);
	for(int i = 0; i < numpolys[onum]; i++){
		double thisx[psize[onum][i]],thisy[psize[onum][i]];

		for(int j = 0; j < psize[onum][i]; j++){
			thisx[j] = xbar[onum][con[onum][i][j]];
			thisy[j] = ybar[onum][con[onum][i][j]];

			px[j] = x[onum][con[onum][i][j]];
			py[j] = y[onum][con[onum][i][j]];
			pz[j] = z[onum][con[onum][i][j]];

		}
    //printa(px, 4);

    //check_vectors(thisx, thisy, psize[onum][i], backface);
		//backface = backface_elimination(px,py,pz,4) ;
    if(backface == 1){
      //G_fill_polygon(thisx, thisy, psize[onum][i]);
      //G_rgb(0,0,0);
      G_polygon(thisx, thisy, psize[onum][i]);
    } else {
      //printf("hidden");
    }
	}
}

int compare (const void *p, const void *q){
  POLY *a, *b ;

  a = (POLY*)p ;
  b = (POLY*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0 ;
}

int paint(){
  double xcenter, ycenter, zcenter;
  double thisx[10], thisy[10];
  int c = 0;

  xcenter = ycenter = zcenter = 0 ;
  //printf("z = %lf\n", zcenter);
  create_2d_lists();
  G_rgb(0,0,0) ;
  G_clear() ;
  for(int i = 0; i < numobjects; i++){
    for(int j = 0; j < numpolys[i]; j++){
      for(int k = 0; k < psize[i][j]; k++){  //printa(px, 4);

    //check_vectors(thisx, thisy, psize[onum][i], backface);
		//backface = backface_elimination(px,py,pz,4) ;
        xcenter += x[i][con[i][j][k]];
        ycenter += y[i][con[i][j][k]];
        zcenter += z[i][con[i][j][k]];
        //printf("z at %d = %lf\n", i, z[i][con[i][j][k]]);

      }
      xcenter/=psize[i][j]; ycenter/=psize[i][j]; zcenter/=psize[i][j];
      zdist[c].dist = zcenter;
      zdist[c].objnum = i;
      zdist[c].polynum = j;
      c++;
    }
    //G_rgb((i*i*i),(100-i),i) ;
    //draw_wire_frame(i);
  }
  //printf("c = %d\n", c);
  qsort(zdist, c, sizeof(POLY), compare);
  //printf("c = %d", c);
  for(int i = 0; i < c; i++){
    //printf("dist at %d = %lf\n", i, zdist[i].dist);
  }


  for(int i = c-1; i >= 0; i--){
    //printf("i = %d\n", i);

    //printf("%d - %lf\n", i, zdist[i].dist);
    int objnum = zdist[i].objnum;
    int polynum = zdist[i].polynum;
    for(int j = 0; j < psize[objnum][polynum]; j++){
        thisx[j] = xbar[objnum][con[objnum][polynum][j]];
        thisy[j] = ybar[objnum][con[objnum][polynum][j]];
    }
    G_rgb(1, 1, 1);
    //printa(thisx, psize[objnum][polynum]);
    //printf()
    G_fill_polygon(thisx, thisy, psize[objnum][polynum]);
    //G_circle(400, 400, 10);
    G_rgb(0,0,0);
    G_polygon(thisx, thisy, psize[objnum][polynum]);
  }


}





int main(int argc, char **argv ){
	read_objects(argc, argv);
  translate_scale();
  swidth = 800; sheight = 800;
  G_init_graphics(swidth, sheight);

  //print_all();

  double t[4][4], rx[4][4], ry[4][4], rz[4][4], V[4][4];

  while (1) {
    M3d_make_identity(V);
    //rotation_angle = sign*rotation_angle;
    M3d_make_rotation_aboutx(rx, cos(rotation_angle), sin(rotation_angle));
    M3d_make_rotation_abouty(ry, cos(rotation_angle), sin(rotation_angle));
    M3d_make_rotation_aboutz(rz, cos(rotation_angle), sin(rotation_angle));

    int q = G_wait_key() ;

    if (q == 'q') {
      exit(0) ;

    } else if (q == 'c') {
      sign = -sign ;
      rotation_angle = -rotation_angle;

    } else if (q == 't') {
      action = q ;

    } else if (q == 'r') {
      action = q ;

    } else if (('0' <= q) && (q <= '9')) {
      key = q - '0' ;
      if (key < numobjects) { onum = key ; }
      //average(onum);

    } else if(q == 's'){
      //printf("backface = %d\n", backface);
      if(s == 1){s = 0;} else {s = 1;}

      //printf("backface = %d\n", backface);

    } else if ((q == 'x') && (action == 't')) {
      M3d_make_translation(V, sign*0.1, 0, 0);
      //average(onum);

    } else if ((q == 'y') && (action == 't')) {
      M3d_make_translation(V, 0, sign*0.1, 0);
      //average(onum);

    } else if ((q == 'z') && (action == 't')) {
      //printf("enter\n");
      M3d_make_translation(V, 0, 0, sign*0.1);
      //average(onum);

    } else if ((q == 'x') && (action == 'r')) {

      M3d_make_translation(t, -xaverage[onum], -yaverage[onum], -zaverage[onum]);
      M3d_mat_mult(V, rx, t);
      M3d_make_translation(t, xaverage[onum], yaverage[onum], zaverage[onum]);
      M3d_mat_mult(V, t, V);

    } else if ((q == 'y') && (action == 'r')) {
      M3d_make_translation(t, -xaverage[onum], -yaverage[onum], -zaverage[onum]);
      M3d_mat_mult(V, ry, t);
      M3d_make_translation(t, xaverage[onum], yaverage[onum], zaverage[onum]);
      M3d_mat_mult(V, t, V);

    } else if ((q == 'z') && (action == 'r')) {
      M3d_make_translation(t, -xaverage[onum], -yaverage[onum], -zaverage[onum]);
      M3d_mat_mult(V, rz, t);
      M3d_make_translation(t, xaverage[onum], yaverage[onum], zaverage[onum]);
      M3d_mat_mult(V, t, V);

    } else {
      printf("no action\n") ;
    }



    M3d_mat_mult_points (x[onum],y[onum],z[onum],  V,
			 x[onum],y[onum],z[onum],numpoints[onum]+1) ;
    average(onum);

    //backface();

    paint();

    /*create_2d_lists(onum);

    G_rgb(0,0,0) ;
    G_clear() ;
    G_rgb(0,0,1) ;
    //    draw_all_objects() ;
    draw_wire_frame(onum);*/
  }
	G_wait_key();
}
