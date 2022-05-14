#include <omp.h>

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <cmath>

#define GENCUT_R_PACKAGE

// [[Rcpp::plugins(openmp)]]

#ifndef SKIP_R_API
#include <Rcpp.h>
#define printf Rprintf
#endif