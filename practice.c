#include "M2d_matrix_tools.c"

//double x[5];
//double y[5];
double a[3][3];

int main(){
	double x[5] = {2, 2, 6, 8, 6};
	double y[5] = {0, -2, 0, 0, -2};

	for(int i = 0; i < 5; i++){
		printf("(%lf, %lf)", x[i], y[i]);
	}
	
	M2d_make_translation(a, 0, 2);
	M2d_mat_mult_points(x, y, a, x, y, 5);

	M2d_make_rotation_degrees(a, 45);
	M2d_mat_mult_points(x, y, a, x, y, 5);

	M2d_make_scaling(a, 1/2, 1/2);
	M2d_mat_mult_points(x, y, a, x, y, 5);

	M2d_make_rotation_degrees(a, (-5));
	M2d_mat_mult_points(x, y, a, x, y, 5);

	M2d_make_translation(a, -2, -1);
	M2d_mat_mult_points(x, y, a, x, y, 5);

	for(int i = 0; i < 5; i++){
		printf("(%lf, %lf)", x[i], y[i]);
	}
}

