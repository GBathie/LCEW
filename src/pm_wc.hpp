/**
 * Perform pattern matching in strings with wildcards
 * using FFT, follwing the work of [Clifford and Clifford].
 */

 #pragma once

#include <vector>
#include <unordered_set>
#include <random>

using std::unordered_set;
using std::vector;

/**
 * Find occurences of `p` in `t` using symbols in `wc` as wildcards.
 *
 * Returns a vector `A` of size `t.size()` s.t. `A[i]` is true if and only
 * if there is an occurrence of `p` in `t` starting at position `i`.
 */
vector<bool> pm_wc(const vector<int> &p, const vector<int> &t, const unordered_set<int> &wc);

void test_pm_wc(int it, std::mt19937 &rng);