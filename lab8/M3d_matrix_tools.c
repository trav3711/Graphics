#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846

/*

 ( x')          (x)
 ( y')  =   M * (y)
 ( z')          (z)
 ( 1 )          (1)

instead of (x',y',z',1) = (x,y,z,1) * M

*/



int M3d_print_mat (double a[4][4])
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           printf(" %12.4lf ",a[r][c]) ;
      }
      printf("\n") ;
  }

  return 1 ;
}





int M3d_copy_mat (double a[4][4], double b[4][4])
// a = b
{
  for(int r=0;r<4;r++){
    for(int c=0;c<4;c++){
      a[r][c]=b[r][c];

    }
  }
}





int M3d_make_identity (double a[4][4])
// a = I
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           if (r == c) a[r][c] = 1.0 ;
               else    a[r][c] = 0.0 ;
      }
  }

  return 1 ;
}





int M3d_make_translation (double a[4][4], double dx, double dy, double dz)
{
  M3d_make_identity(a) ;
  a[0][3] =  dx ;  a[1][3] = dy ; a[2][3]= dz ;
  return 1 ;
}





int M3d_make_scaling (double a[4][4], double sx, double sy, double sz)
{
  M3d_make_identity(a);
  a[0][0]=sx; a[1][1]=sy; a[2][2]=sz;
}


int M3d_make_rotation_aboutx (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known

{
  M3d_make_identity(a);
  a[1][1]=cs; a[1][2]=-sn;
  a[2][1]=sn; a[2][2]= cs;
}

int M3d_make_rotation_abouty (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known

{
  M3d_make_identity(a);
  a[0][0]=cs; a[0][2]=-sn;
  a[2][0]=sn; a[2][2]= cs;
}

int M3d_make_rotation_aboutz (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known

{
  M3d_make_identity(a);
  a[0][0]=cs; a[0][1]=-sn;
  a[1][0]=sn; a[1][1]= cs;
}

//already know cs and sn, so do not need angle



int M3d_mat_mult (double result[4][4], double a[4][4], double b[4][4])
// res = a * b
// this is SAFE, i.e. the user can make a call such as
// M2d_mat_mult(p,  p,q) or M2d_mat_mult(p,  q,p) or  M2d_mat_mult(p, p,p)
{
  double sum=0;
  double res[4][4];

  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int c=0; c<4;c++){
	sum+=(a[i][c]*b[c][j]);
      }
      res[i][j]=sum;
      sum=0;
    }
  }
  for(int i=0;i<4; i++){
    for(int j=0;j<4;j++){
      result[i][j]=res[i][j];
    }
  }
}

//multiplies two 3*3 and puts into a new 3*3


//safe means it should still work even if some elements are identical
//so as in example above, p*q can be put into p, use a temporary matrix then use copy

//basically 3x1s

int M3d_mat_mult_pt (double P[3],   double m[4][4], double Q[3])
// P = m*Q
// SAFE, user may make a call like M2d_mat_mult_pt (W, m,W) ;
{
  double u,v,t;

  u=m[0][0]*Q[0] +m[0][1]*Q[1]+m[0][2]*Q[2]+m[0][3];
  v=m[1][0]*Q[0] +m[1][1]*Q[1]+m[1][2]*Q[2]+m[1][3];
  t=m[2][0]*Q[0] +m[2][1]*Q[1]+m[2][2]*Q[2]+m[2][3];

  P[0] = u;
  P[1] = v;
  P[2] = t;


  return 1;

}

//takes whole matrix and multiplies by one point
//



int M3d_mat_mult_points (double *X, double *Y, double *Z,
                         double m[4][4],
                         double *x, double *y,double *z, int numpoints)
// |X0 X1 X2 ...|       |x0 x1 x2 ...|
// |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
// |Z0 Z1 Z2 ...|       |z0 z1 z2 ...\
// | 1  1  1 ...|       | 1  1  1 ...|

// SAFE, user may make a call like M2d_mat_mult_points (x,y, m, x,y, n) ;
{

  double u,v,t;

  int i;

  for (i=0;i<numpoints;i++){
    u=m[0][0]*x[i] +m[0][1]*y[i]+m[0][2]*z[i]+m[0][3];
    v=m[1][0]*x[i] +m[1][1]*y[i]+m[1][2]*z[i]+m[1][3];
    t=m[2][0]*x[i] +m[2][1]*y[i]+m[2][2]*z[i]+m[2][3];

    X[i]=u;
    Y[i]=v;
    Z[i]=t;



  }

  return 1;




}


//multiplies matrix by all the data and gives new data
//numpoints is important here, depends on how many numbers there are

/*
int main (){
  double x[5] = {1,5,2,4,6};
  double y[5] = {2,7,9,3,2};
  int numpoints=5;


  double  a[3][3] = {2, 1, 3,
		     3, 5, 6,
                     1, 1, 1};


  M2d_mat_mult_points (x, y,a,x,y,numpoints);

  printf("%lf, %lf, %lf, %lf, %lf\n",x[0],x[1],x[2],x[3],x[4]);
  printf("%lf, %lf, %lf, %lf, %lf\n",y[0],y[1],y[2],y[3],y[4]);




}
*/
