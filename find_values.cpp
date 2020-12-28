#include "find_values.h"
#include "io.h"
#include <iostream>

using namespace std;
#define EPS 1e-50

/*
void three_diag(double* mat, int n) {
    int i; int j; int k;
	double norm_a;
	double norm;

	for (i = 0; i < n - 2; ++i)
	{
		norm_a = 0;
		for (j = i + 2; j < n; ++j)
			norm_a += mat[ind(i, j)] * mat[ind(i, j)];

		norm = sqrt(mat[ind(i +1, i)] * mat[ind(i + 1, i)] + norm_a);
		mat[ind(i + 1, i)] -= norm;

		norm_a = 1.0 / sqrt(mat[ind(i + 1, i)] * mat[ind(i + 1, i)] + norm_a);
		for (j = i + 1; j < n; ++j)
			mat[ind(j, i)] *= norm_a;

		for (j = i + 1; j < n; ++j)
		{
			norm_a = 0;
			for (k = i + 1; k < n; ++k)
				norm_a += mat[ind(k, i)] * mat[ind(k, j)];

			norm_a *= 2.0;
			for (k = i + 1; k < n; ++k)
				mat[ind(k, j)] -= norm_a * mat[ind(k, i)];
		}

		for (j = i + 1; j < n; ++j)
		{
			norm_a = 0;
			for (k = i + 1; k < n; ++k)
				norm_a += mat[ind(k, i)] * mat[ind(j, k)];

			norm_a *= 2.0;
			for (k = i + 1; k < n; ++k)
				mat[ind(j, k)] -= norm_a * mat[ind(k, i)];
		}

		mat[ind(i + 1, i)] = norm;
		for (j = i + 2; j < n; ++j)
			mat[ind(j, i)] = 0;
	}
}*/



void three_diag(double* mat, int n) {
    double x, y, r,  cos, sin;
	for (int i = 1; i < n - 1; ++i) {
		for (int j = i + 1; j < n; ++j) {
			x = mat[ind(i - 1, i)];
			y = mat[ind(i - 1, j)];
			r = sqrt(x * x + y * y);
			
			if (fabs(y) < EPS) {
				continue;
			}
			if (fabs(r) < EPS ){
				continue;
			}

			cos = x / r;
			sin = -1 * y / r;

		 	mat[ind(i - 1, i)] = r;
			mat[ind(i - 1, j)] = 0;

			for (int k = i + 1; k < n; ++k) {
				if (k == j) continue;
				x = mat[ind(i, k)];
				y = mat[ind(j, k)];
				mat[ind(i, k)] = x * cos - y * sin;
				mat[ind(j, k)] = x * sin + y * cos;
			}

			x = mat[ind(i, i)];
			y = mat[ind(j, j)];
			r = mat[ind(i, j)];
			mat[ind(i, i)] = (x * cos - r * sin) * cos - (r * cos - y * sin) * sin;
			mat[ind(i, j)] = (x * cos - r * sin) * sin + (r * cos - y * sin) * cos;
			mat[ind(j, j)] = (x * sin + r * cos) * sin + (r * sin + y * cos) * cos;
		}
	}
}

int n_(double* mat, int n, double l) {
	int res = 0;
	double a;
	a = mat[0] - l;
	if (a < 0) res = 1;

	for (int i = 1; i < n; ++i) {
		if (fabs(a) < EPS){
			//cout << "FASB " << fabs(a) << endl;
			return -1;
		}
		a = mat[ind(i, i)] - l - mat[ind(i - 1, i)] * mat[ind(i - 1, i)] / a;
		if (a < 0) ++res;
	}
	return res;
}

void find_values(double* mat, int n, double* res, double eps) {
    int count = 0;
    long double from, to, no;
    three_diag(mat, n);
	//cout << "DIAG: " << endl;
	//print(mat, n, n, n, cout);
	no = norm(mat, n);
    from = -no - 2 * eps; to = no + 2 * eps;

    int count_left = n_(mat, n, from);
    int count_vals = n_(mat, n, to) - count_left;

    int i = 0;
    long double cur_from = from, cur_to = to, cur_count = 0, cur_m = 0;
	double cur_sgn = 0;
	int tmp = 0, done = 0;
    while (i < count_vals) {   
		done = 0;
        while (!done && (cur_to - cur_from) > eps) {
            cur_m = (cur_from + cur_to) / 2.0;
            //cout << cur_from << " " << cur_to << " " << (cur_to - cur_from)  << endl;
            //cout << "CUR M " << cur_m << endl;
			tmp = n_(mat, n, cur_m);
			//cout << "CUR M N_ " << cur_m << " " << tmp << endl;
			if (tmp == -1) {
				//cout << "SDF " << eps << " " << cur_from << " " << cur_m << " " << cur_to << " " << fabs(cur_to - cur_from) << endl;
				cur_sgn = 0;
				if (cur_m + 10 * eps < to) {
					cur_sgn = 1;
					cout << "R" << endl;
				} else {
					cur_sgn = -1;
					cout << "L" << endl;
				}
				while (tmp == -1) {
					cur_m += cur_sgn * eps;
					tmp = n_(mat, n, cur_m);
				}
				done = 1;
			} else if (tmp < i + 1){
                cur_from = cur_m;
            } else {
                cur_to = cur_m;
            }
            count++;
        }

        cur_m = 0.5 * (cur_from + cur_to);
		cur_count = n_(mat, n, cur_to) - n_(mat, n, cur_from);
		//cout << "CUR COUNT " << cur_count << " " << cur_m << endl;

		for (int j = 0; i + j < n && j < cur_count; ++j) {
			res[i + j] = cur_m;	
		}

		i += cur_count;

		cur_from = cur_m;
		cur_to = to;
    }
    cout << "count: " << count << endl;
}