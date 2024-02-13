#include <string>
#include <random>
#include <chrono>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <iostream>
#include <stdlib.h>

#include "pm_wc.hpp"
#include "lcew.hpp"
#include "fast_mm.hpp"

using namespace std;

template <class RNG>
vector<vector<bool>> random_bool_mat(int n, double density, RNG &rng)
{
    bernoulli_distribution b(density);
    vector<vector<bool>> res(n, vector<bool>(n, false));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            res[i][j] = b(rng);

    return res;
}

vector<vector<bool>> mat_mul(vector<vector<bool>> &a, vector<vector<bool>> &b)
{
    size_t n = a.size();
    vector<vector<bool>> res(n, vector<bool>(n, false));
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            for (size_t k = 0; k < n; k++)
                res[i][j] = res[i][j] || (a[i][k] && b[k][j]);
    return res;
}

template <class RNG>
double test_mm(size_t n, double p, RNG &rng)
{
    vector<vector<bool>> a = random_bool_mat(n, p, rng);
    vector<vector<bool>> b = random_bool_mat(n, p, rng);
    SparseBoolMatrix a_s = SparseBoolMatrix::from_dense(a);
    SparseBoolMatrix b_s = SparseBoolMatrix::from_dense(b);

    chrono::steady_clock sc;
    auto check_dense = mat_mul(a, b);
    SparseBoolMatrix check = SparseBoolMatrix::from_dense(check_dense);

    auto start = sc.now();
    SparseBoolMatrix res = matrix_mult(a_s, b_s);
    auto end = sc.now();
    auto time_rep = static_cast<chrono::duration<double>>(end - start).count();
    assert(res == check);

    return time_rep;
}

const size_t nb_rep = 10;

/**
 * Generate timing data
 */
void run_tests(size_t max_size, size_t step)
{
    random_device rd;
    mt19937 rng(rd());
    for (float p = 1; p <= 10; p += 1)
    {
        string name_file = "results_" + to_string(p) + ".csv";
        ofstream outFile;
        outFile.open(name_file);

        if (outFile.is_open())
        {
            for (size_t n = step; n <= max_size; n += step)
            {
                for (size_t i = 0; i < nb_rep; i++)
                {
                    cout << n << "," << p << "," << i << endl;

                    double time_total = test_mm(n, p / n, rng);
                    outFile << n << "," << p << "," << time_total << endl;
                }
            }
            outFile.close();
        }
        else
        {
            cout << "Unable to open files: " << name_file << endl;
        }
    }
}

int main()
{
    size_t max_size = 2000;
    size_t step = 500;

    run_tests(max_size, step);
    return 0;
}