#include <stdio.h>
#include <math.h>


/*

 ( x')          (x)
 ( y')  =   M * (y)  
 ( 1 )          (1)

instead of (x',y',1) = (x,y,1) * M  

*/



int M2d_print_mat (double a[3][3])
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           printf(" %12.4lf ",a[r][c]) ;
      }
      printf("\n") ;
  }

  return 1 ;
} 


int M2d_copy_mat (double a[3][3], double b[3][3])
// a = b
{
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			b[i][j] = a[i][j];
		}
	}
	return 1;
} 


int M2d_make_identity (double a[3][3])
// a = I
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           if (r == c) a[r][c] = 1.0 ;
               else    a[r][c] = 0.0 ;
      }
  }

  return 1 ;
} 


int M2d_make_translation (double a[3][3], double dx, double dy)
{
  M2d_make_identity(a) ;
  a[0][2] =  dx ;  
  a[1][2] = dy ;

  return 1 ;
}

int M2d_make_scaling (double a[3][3], double sx, double sy)
{
	M2d_make_identity(a);
	a[0][0] = sx;
	a[1][1] = sy;

	return 1;
}


int M2d_make_rotation_radians (double a[3][3],  double r)
//just making the rotation matrix in radians
{
	M2d_make_identity(a);

	a[0][0] = cos(r);
	a[0][1] = (-1)*sin(r);
	a[1][0] = sin(r);
	a[1][1] = cos(r);

	return 1;
}


int M2d_make_rotation_degrees (double a[3][3],  double d)
{
	double r = (d*M_PI)/180.0;
	M2d_make_rotation_radians(a, r);

	return 1;
}


int M2d_make_rotation_cs (double a[3][3], double cs, double sn)
// this one assumes cosine and sine are already known
{
	M2d_make_identity(a);

	a[0][0] = cs;
	a[0][1] = sn;
	a[0][1] = (-1)*sn;
	a[1][1] = cs;

	return 1;
}

int M2d_mat_mult (double res[3][3], double a[3][3], double b[3][3])
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// M2d_mat_mult(p,  p,q) or M2d_mat_mult(p,  q,p) or  M2d_mat_mult(p, p,p)
{
	double buff[3][3];

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			for(int c = 0; c < 3; c++){
				buff[i][j] += (a[i][c]*b[c][j]); 
			}
		}
	}

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			res[i][j] = buff[i][j];
		}
	}

	return 1;
}


int M2d_mat_mult_pt (double P[2],   double m[3][3], double Q[2])
// P = m*Q
// SAFE, user may make a call like M2d_mat_mult_pt (W, m,W) ;
{
	double result[3][3] = { {0, 0, 0},
							{0, 0, 0},
							{0, 0, 0} };	
	
	double temp[3][3];
	M2d_make_identity(temp);
	temp[0][2] = Q[0];
	temp[1][2] = Q[1];

	M2d_mat_mult(result, m, temp);
	
	P[0] = result[0][2];
	P[1] = result[1][2];

	return 1;
}

int M2d_mat_mult_points (double *X, double *Y,
                         double m[3][3],
                         double *x, double *y, int numpoints)
// |X0 X1 X2 ...|       |x0 x1 x2 ...|
// |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
// | 1  1  1 ...|       | 1  1  1 ...|

// SAFE, user may make a call like M2d_mat_mult_points (x,y, m, x,y, n) ;
{
	for(int i = 0; i < numpoints; i++){
		double temp[2] = {x[i], y[i]};
		double result[2];

		M2d_mat_mult_pt(result, m , temp);

		X[i] = result[0]; Y[i] = result[1];
	}

	return 1;
}
