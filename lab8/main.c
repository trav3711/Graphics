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
int frame = 1;
double hither = 5;
double yon = 100;

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

int inita(double *a, int n){
  for(int i = 0; i < n; i++){
    a[i] = 0;
  }
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

double vector_length(double vx, double vy, double vz){
  double length;
  double sum = pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0);
  length = sqrt(sum);
  return length;
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

    average(i);
  }
}

int create_2d_lists(){
  double x1, y1;

	for(int i = 0; i < numobjects; i++){

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
int vector_length(double vx, double vy, double vz){
  double length = 0;
  length = sqrt(pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0));
  return length;
}
  Ax = px[0] - xcenter; Bx = px[1] - xcenter;
  Ay = py[0] - ycenter; By = py[1] - ycenter;
  Az = pz[0] - zcenter; Bz = pz[1] - zcenter;

  Nx = (Ay*Bz) - (By*Az); Ny = (Bx*Az) - (Ax*Bz); Nz = (Ax*By) - (Bx*Ay);
  Ex = -xcenter; Ey = -ycenter; Ez = -zcenter;

  dot = (Nx*Ex) + (Ny*Ey) + (Nz*Ez);

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
    }
    fscanf(fp, "%d", &numpolys[i]);
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

	for(int i = 0; i < numpolys[onum]; i++){
		double thisx[psize[onum][i]],thisy[psize[onum][i]];

		for(int j = 0; j < psize[onum][i]; j++){
			thisx[j] = xbar[onum][con[onum][i][j]];
			thisy[j] = ybar[onum][con[onum][i][j]];

			px[j] = x[onum][con[onum][i][j]];
			py[j] = y[onum][con[onum][i][j]];
			pz[j] = z[onum][con[onum][i][j]];

		}
    if(backface == 1){
      G_polygon(thisx, thisy, psize[onum][i]);
    } else {
      //printf("hidden");
    }
	}
}

double set_color(double *px, double *py, double *pz, int size){
  double sourcex, sourcey, sourcez;
  double Intensity, ambient, diffusemax, specular, specpow;
  double NuLu, RuEu, NuEu;

  double A[3], B[3], center[3], N[3], E[3], R[3], L[3];
  double Nlength, Elength, Llength;

  inita(A, 3); inita(B, 3); inita(center, 3);
  inita(N, 3); inita(E, 3); inita(R, 3); inita(L, 3);

  sourcex = 100; sourcey = 200; sourcez = 0;
  Nlength = Elength = Llength = 0;
  NuLu = RuEu = NuEu = 0;
  ambient = 0.2; diffusemax = 0.5; specpow = 50;
  specular = 1.0 - ambient - diffusemax;

  for(int i = 0; i < size; i++){
    center[0] += px[i];
    center[1] += py[i];
    center[2] += pz[i];
  }
  center[0]/=size; center[1]/=size; center[2]/=size;

  A[0] = px[0] - center[0]; B[0] = px[1] - center[0];
  A[1] = py[0] - center[1]; B[1] = py[1] - center[1];
  A[2] = pz[0] - center[2]; B[2] = pz[1] - center[2];

  N[0] = (A[1]*B[2]) - (B[1]*A[2]); N[1] = (B[0]*A[2]) - (A[0]*B[2]); N[2] = (A[0]*B[1]) - (B[0]*A[1]);
  E[0] = -center[0]; E[1] = -center[1]; E[2] = -center[2];
  L[0] = sourcex-center[0]; L[1] = sourcey-center[1]; L[2] = sourcez-center[2];

  Nlength = vector_length(N[0], N[1], N[2]);
  Elength = vector_length(E[0], E[1], E[2]);
  Llength = vector_length(L[0], L[1], L[2]);

  N[0] /= Nlength; N[1] /= Nlength; N[2] /= Nlength;
  E[0] /= Elength; E[1] /= Elength; E[2] /= Elength;
  L[0] /= Llength; L[1] /= Llength; L[2] /= Llength;

  NuLu = (N[0]*L[0]) + (N[1]*L[1]) + (N[2]*L[2]);
  NuEu = (N[0]*E[0]) + (N[1]*E[1]) + (N[2]*E[2]);

  if((NuLu < 0) && (NuEu < 0)){
    N[0]*=-1; N[1]*=-1; N[2]*=-1;
    NuEu = -NuEu;
    NuLu = -NuLu;
  } else if(((NuLu < 0) && (NuEu > 0)) || ((NuLu > 0) && (NuEu < 0))){
    return ambient;
  }
  R[0] = (2*(NuLu)*N[0]) - L[0];
  R[1] = (2*(NuLu)*N[1]) - L[1];
  R[2] = (2*(NuLu)*N[2]) - L[2];

  RuEu = (R[0]*E[0]) + (R[1]*E[1]) + (R[2]*E[2]);
  Intensity = ambient + (diffusemax*NuLu) + (specular*pow(RuEu, specpow));
  return Intensity;

}

int clip_poly_on_plane(double a, double b, double c, double d,
                          double *px, double *py, double *pz, int size,
                        double *resx, double *resy, double *resz)
{
  int num, i, j;
  double x1, y1, z1, x2, y2, z2, xint, yint, zint;
  double s1, s2, x21, y21, z21;
  double constants, terms, t;

  x1 = y1 = z1 = 0;
  x2 = y2 = z2 = 0;
  x21 = y21 = z21 = 0;
  xint = yint = zint = 0;
  s1 = s2 = 0;

  for(int i = 0; i < size; i++){
    j = (i + 1) % 3;

    // load up segment to be clipped
    x1 = px[i]; x2 = px[j];
    y1 = py[i]; y2 = py[j];
    z1 = pz[i]; z2 = pz[j];

    // clip line segment (x1,y1)-(x2,y2) against line
    s1 = (a*x1 + b*y1 + c*z1 + d) ;
    s2 = (a*x2 + b*y2 + c*z2 + d) ;

    if ((s1 >= 0) && (s2 >= 0)) {
       // out to out, do nothing
    } else if ((s1 < 0) && (s2 < 0)) {
       // in to in
       resx[num] = x2 ; resy[num] = y2 ; num++ ;
    } else {
        // one is in, the other out, so find the intersection
        //xint = x2 + t*x21
        x21 = x2-x1; y21 =y2-y1; z21 = z2-z1;
        constants = (a*x2) + (b*y2) + (c*z2) + d;
        terms = (a*x21) + (b*y21) + (c*z21);
        t = -constants/terms;

        xint = x2 + (t*x21);
        yint = y2 + (t*y21);
        zint = z2 + (t*z21);

        if (s1 < 0) {
          // in to out
          resx[num] = xint; resy[num] = yint; resz[num] = zint; num++ ;
        } else  {
          // out to in
          resx[num] = xint; resy[num] = yint; resz[num] = zint; num++;
          resx[num] = x2; resy[num] = y2; resz[num] = z2; num++ ;
        }
    }
  }
  return num;
}

int clip_poly_on_every_plane(double *polyx, double *polyy, double *polyz, int size){
  //double yonwin[3][4]; //inita(yonwin);
  //double hitherwin[3][4]; //inita(hitherwin);
  double nx[100], ny[100], nz[100];
  //double thisx[3], thisy[3], thisz[3];
  int wn = 6;
  double a[wn], b[wn], c[wn], d[wn];

  inita(a, wn); inita(b, wn); inita(c, wn); inita(d, wn);

  /*double yonheight = tan(viewangle*(M_PI/180))*yon;
  double hitherheight = tan(viewangle*(M_PI/180))*hither;

  yonwin[0][0] = yonheight;
  yonwin[1][0] = yonheight;
  yonwin[2][0] = yon;

  yonwin[0][1] = -yonheight;
  yonwin[1][1] = yonheight;
  yonwin[2][1] = yon;

  yonwin[0][2] = -yonheight;
  yonwin[1][2] = -yonheight;
  yonwin[2][2] = yon;

  yonwin[0][3] = yonheight;
  yonwin[1][3] = -yonheight;
  yonwin[2][3] = yon;

  hitherwin[0][0] = hitherheight;
  hitherwin[1][0] = hitherheight;
  hitherwin[2][0] = hither;

  hitherwin[0][1] = -hitherheight;
  hitherwin[1][1] = hitherheight;
  hitherwin[2][1] = hither;
// one is in, the other out, so find the intersection
  hitherwin[0][2] = -hitherheight;
  hitherwin[1][2] = -hitherheight;
  hitherwin[2][2] = hither;

  hitherwin[0][3] = hitherheight;
  hitherwin[1][3] = -hitherheight;
  hitherwin[2][3] = hither;*/

  double v = tan(viewangle)*yon;

  a[0] = 0; a[1] =  1; a[2] =  0; a[3] = 1; a[4] = 0;   a[5] = 0;
  b[0] = 1; b[1] =  0; b[2] =  1; b[3] = 0; b[4] = 0;   b[5] = 0;
  c[0] = v; c[1] = -v; c[2] = -v; c[3] = v; c[4] = 1;   c[5] = 1;
  d[0] = 0; d[1] =  0; d[2] =  0; d[3] = 0; d[4] = yon; d[5] = hither;

  for(int i = 0; i < wn; i++){

    size = clip_poly_on_plane(a[i], b[i], c[i], d[i], polyx, polyy, polyz, size, nx, ny, nz);

    for(int j = 0; j < size; j ++){
      polyx[i] = nx[i]; polyy[i] = ny[i]; polyz[i] = nz[i];
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
  double thisx[10], thisy[10], rgb[3];
  double px[20],py[20],pz[20] ;
  int c = 0;

  xcenter = ycenter = zcenter = 0 ;
  inita(rgb, 3);
  create_2d_lists();
  G_rgb(0,0,0) ;
  G_clear() ;
  for(int i = 0; i < numobjects; i++){

    for(int j = 0; j < numpolys[i]; j++){

      xcenter = ycenter = zcenter = 0 ;
      for(int k = 0; k < psize[i][j]; k++){
        xcenter += x[i][con[i][j][k]];
        ycenter += y[i][con[i][j][k]];
        zcenter += z[i][con[i][j][k]];
      }
      xcenter/=psize[i][j]; ycenter/=psize[i][j]; zcenter/=psize[i][j];
      zdist[c].dist = zcenter;

      zdist[c].objnum = i;
      zdist[c].polynum = j;
      c++;
    }
  }
  qsort(zdist, c, sizeof(POLY), compare);
  for(int i = 0; i < c; i++){
  }

  for(int i = c-1; i >= 0; i--){
    int objnum = zdist[i].objnum;
    int polynum = zdist[i].polynum;
    for(int j = 0; j < psize[objnum][polynum]; j++){
        thisx[j] = xbar[objnum][con[objnum][polynum][j]];
        thisy[j] = ybar[objnum][con[objnum][polynum][j]];

        px[j] = x[objnum][con[objnum][polynum][j]];
        py[j] = y[objnum][con[objnum][polynum][j]];
        pz[j] = z[objnum][con[objnum][polynum][j]];
    }
    int size = psize[objnum][polynum];

    int nsize = clip_poly_on_every_plane(px, py, pz, size);

    if(objnum == 0){
      rgb[0] = 1; rgb[1] = 0; rgb[2] = 0;
    } else if(objnum == 1){
      rgb[0] = 0; rgb[1] = 1; rgb[2] = 0;
    } else if( objnum == 2){
      rgb[0] = 0; rgb[1] = 0; rgb[2] = 1;
    } else{G_rgb(1, 1, 1);}
    double intensity, colorratio, colorratiobase;

    colorratiobase = 0.7;

    intensity = set_color(px, py, pz, nsize);

    colorratio = intensity/colorratiobase;
    G_rgb(rgb[0]*colorratio, rgb[1]*colorratio, rgb[2]*colorratio);

    G_fill_polygon(thisx, thisy, nsize);
    if(frame > 0){
      G_rgb(0,0,0);
      G_polygon(thisx, thisy, nsize);
    }
  }
}

int main(int argc, char **argv ){
	read_objects(argc, argv);
  translate_scale();
  swidth = 800; sheight = 800;
  G_init_graphics(swidth, sheight);

  double t[4][4], rx[4][4], ry[4][4], rz[4][4], V[4][4];

  while (1) {
    M3d_make_identity(V);
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

    } else if(q == 's'){
      if(s == 1){s = 0;} else {s = 1;}

    } else if(q == 'w'){
      frame*=-1;
    } else if ((q == 'x') && (action == 't')) {
      M3d_make_translation(V, sign*0.1, 0, 0);

    } else if ((q == 'y') && (action == 't')) {
      M3d_make_translation(V, 0, sign*0.1, 0);

    } else if ((q == 'z') && (action == 't')) {
      M3d_make_translation(V, 0, 0, sign*0.1);

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
