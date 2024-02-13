#include "fast_mm.hpp"
#include "lcew.hpp"
#include <cmath>

/**
 * Write a SparseBoolMatrix to a string
 * in which LCEW is the same as matrix multiplication.
 *
 * Appends a (dense) string representation of the sparse boolean matrix `a`
 * at the end of the vector `v`.
 * If `lhs` is true, `a` is written in row-major order,
 * and other in column-major order.
 */
void convert_to_string(vector<int> &v, const SparseBoolMatrix &a, bool lhs = true)
{
    int n = a.n;
    vector<int> res(n * n, DEFAULT_WILDCARD);
    for (auto &[i, j] : a.entries)
    {
        if (lhs)
        {
            res[n * i + j] = 1;
        }
        else
        {
            res[n * j + i] = 2;
        }
    }

    v.insert(v.end(), res.begin(), res.end());
}

SparseBoolMatrix matrix_mult(const SparseBoolMatrix &a, const SparseBoolMatrix &b)
{
    vector<int> txt;
    convert_to_string(txt, a, true);
    convert_to_string(txt, b, false);
    int n = a.n;
    // The way that we choose `t` here differs from the paper.
    // First, we use `m_in` as an upper bounds on G,
    // we replace `n+m_out` with `n` because we do not have an estimate of `m_out`,
    // and we add the multiplicative constant `10`
    // to balance the trade-off between preprocessing and query time beyond the big-O.
    int t = 50 * n * sqrt((a.entries.size() + b.entries.size()) / n) + 1;
    Lcew ds(txt, t);

    SparseBoolMatrix res;
    res.n = n;

    auto compute_diag = [&](int i, int j)
    {
        int offset = n * n;
        int l = 0;
        while (i + l < n && j + l < n)
        {
            int r = ds.lcew(n * (i + l), offset + n * (j + l));
            l += r / n;
            if (i + l < n && j + l < n)
                res.entries.emplace_back(i + l, j + l);

            l += 1;
        }
    };
    for (int i = 0; i < n; i++)
        compute_diag(i, 0);

    for (int j = 1; j < n; j++)
        compute_diag(0, j);

    std::sort(res.entries.begin(), res.entries.end());

    return res;
}

SparseBoolMatrix SparseBoolMatrix::from_dense(vector<vector<bool>> &v)
{
    SparseBoolMatrix res;
    size_t n = v.size();
    res.n = n;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            if (v[i][j])
                res.entries.emplace_back(i, j);

    return res;
}
