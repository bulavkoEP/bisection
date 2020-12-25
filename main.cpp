#include "io.h"
#include "find_values.h"
#ifndef MATR
#include "functions.h"
#define MATR
#endif

#include <time.h>
#include <cmath>

int main(int argc, char * argv[]) 
{
    string name = "";
    int n, m, k;
    double eps;
    if (argc < 4) 
        return -1;
    
    try {
        float f;
        n = stoi(argv[1]);
        f = stod(argv[1]);
        if (abs(f - (float) n) > 1e-15) {
            cout << "incorrect argument" << endl;
            return 0;
        }
        m = stoi(argv[2]);
        f = stod(argv[2]);
        if (abs(f - (float) m) > 1e-15) {
            cout << "incorrect argument" << endl;
            return 0;
        }
        eps = stod(argv[3]);

        k = stoi(argv[4]);
        f = stof(argv[4]);
        if (abs(f - (float) k) > 1e-15) {
            cout << "incorrect argument" << endl;
            return 0;
        }
    } catch (const invalid_argument& ia) {
        cout << "incorrect argument" << endl;
        return 0;
    }
    if (k < 0 || k > 5 || n <= 0 || m <= 0 || eps <= 0)
    {
        cout << "incorrect argument" << endl;
        return 0;
    }
    double* mat = new double[n * n];
    if (argc == 6 && k == 0) {
        name = argv[5];
        read_matrix_from_file(mat, n, name);
    } 
    else if (argc == 5 || k != 0) {
        generate_matrix_from_formula(mat, n, k);
    }
    cout << "Matrix:" << endl;
    print(mat, n, n, m, cout);
    cout << endl;

    double* res = new double[n];

    double trace_init = 0;
    double norm_init = 0;
    for (int i = 0; i < n; ++i) trace_init += mat[i * n + i];
    for (int i = 0; i < n * n; ++i) norm_init += mat[i] * mat[i];
    norm_init = sqrt(norm_init);

    clock_t t_start = clock();
    cout << "JSDF" << endl;
    find_values(mat, n, res, eps);

    cout << "time: " << (double) (clock() - t_start) / CLOCKS_PER_SEC << endl;

    cout << "RES: " << endl;
    print(res, n, 1, m, cout);
    cout << endl;

    double trace_after = 0, norm_after = 0; 
    for (int i = 0; i < n; ++i) trace_after += res[i];
    for (int i = 0; i < n; ++i) norm_after += res[i] * res[i];
    norm_after = sqrt(norm_after);
    cout << "norm 1: " << abs(trace_init - trace_after) << endl;
    cout << "norm 2: " << abs(norm_init - norm_after) << endl;

    delete[] mat; delete[] res;
}
