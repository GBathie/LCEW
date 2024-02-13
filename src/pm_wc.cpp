#include "pm_wc.hpp"
#include "ntt.hpp"
#include <iostream>
#include <cassert>

template <class T, class F>
vector<unsigned> vec_map(const vector<T> &v, F &&f)
{
    vector<unsigned> res(v.size());
    for (size_t i = 0; i < v.size(); i++)
    {
        res[i] = f(v[i]);
    }
    return res;
}

vector<bool> pm_wc(const vector<int> &pat, const vector<int> &text, const unordered_set<int> &wc)
{
    int n = text.size();
    int m = pat.size();
    // Flip t before FFT
    vector<int> t_rev(text.rbegin(), text.rend());
    auto wc_zero = [&](int c)
    { return wc.contains(c) ? 0 : c; };
    vector<unsigned> p = vec_map(pat, wc_zero);
    vector<unsigned> t = vec_map(t_rev, wc_zero);

    auto p3 = vec_map(p, [](int i)
                      { return i * i * i; });
    auto p2 = vec_map(p, [](int i)
                      { return i * i; });
    auto t3 = vec_map(t, [](int i)
                      { return i * i * i; });
    auto t2 = vec_map(t, [](int i)
                      { return i * i; });
    auto conv1 = conv(p3, t);
    auto conv2 = conv(p2, t2);
    auto conv3 = conv(p, t3);

    for (int i = 0; i < n; i++)
    {
        conv1[i] += -2 * conv2[i] + conv3[i];
    }

    vector<bool> res(n, false);
    for (int j = m - 1; j < n; j++)
    {
        res[n - j - 1] = conv1[j] == 0;
    }

    return res;
}
