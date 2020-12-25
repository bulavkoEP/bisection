#include "find_values.h"
#define EPS 1e-20

void three_diag(double* mat, int n) {
    double x, y, r, s, cos, sin;
	for (int i = 1; i < n - 1; ++i) {
		for (int j = i + 1; j < n; ++j) {
			x = mat[i * n + i - 1];
			y = mat[j * n + i - 1];
			r = sqrt(x * x + y * y);
			if (fabs(y) < EPS) continue;

			cos = x / r;
			sin = -y / r;

		 	mat[i * n + i - 1] = r;
			mat[(i - 1) * n + i] = r;
			mat[j * n + i - 1] = 0;
			mat[(i - 1) * n + j] = 0;

			for (int k = i + 1; k < n; ++k) {
				if (k == j) continue;
				x = mat[i * n + k];
				y = mat[j * n + k];
				mat[k * n + i] = mat[i * n + k] = x * cos - y * sin;
				mat[k * n + j] = mat[j * n + k] = x * cos + y * sin;
			}

			x = mat[i * n + i];
			y = mat[j * n + j];
			r = mat[i * n + j];
			s = mat[j * n + i];

			mat[i * n + i] = (x * cos - s * sin) * cos - (r * cos - y * sin) * sin;
			mat[j * n + i] = (x * cos - s * sin) * sin + (r * cos - y * sin) * cos;
			mat[i * n + j] = mat[j * n + i];
			mat[j * n + j] = (x * sin + s * cos) * sin + (r * sin + y * cos) * cos;
		}
	}
}

int n_(double* mat, int n, double l) {
	int res = 0;
	double a;
	a = mat[0] - l;
	if (a < 0) res = 1;

	for (int i = 1; i < n; ++i) {
		a = mat[i * n + i] - l - mat[i * n + i - 1] * mat[(i - 1) * n + i] / a;
		if (a < 0) ++res;
	}
	return res;
}

void find_values(double* mat, int n, double* res, double eps) {
    int count = 0;
    double from, to, no;
    three_diag(mat, n);
	no = norm(mat, n);
    from = -no - 0.1; to = no + 0.1;
	//cout << "SdfSD" << endl;

    int count_left = n_(mat, n, from);
    int count_vals = n_(mat, n, to) - count_left;

    int i = 0;
    double cur_from = from, cur_to = to, cur_count = 0, cur_m = 0;
    while (i < count_vals) {   
        while (cur_to - cur_from > eps) {
            cur_m = 0.5 * (cur_from + cur_to);
            //cout << cur_from << " " << cur_to << " " << cur_to - cur_from << endl;
            //cout << "CUR M " << cur_m << endl;
            if (n_(mat, n, cur_m) < i + 1){
                cur_from = cur_m;
                //cout << "H" << endl;
            } else {
                cur_to = cur_m;
                //cout << "R" << endl;
            }
            count++;
        }

        cur_m = 0.5 * (cur_from + cur_to);
		cur_count = n_(mat, n, cur_to) - n_(mat, n, cur_from);

		for (int j = 0; j < cur_count; ++j)
			res[i + j] = cur_m;

		i += cur_count;

		cur_from = cur_m;
		cur_to = to;
    }
    cout << "count: " << count << endl;
}