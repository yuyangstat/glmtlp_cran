/**********
    C++ Routines for initialization.

    Copyright (C) 2021  Chunlin Li, Yu Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.
**********/

#include <algorithm>
#include <cmath>
#include "utils.h"

inline double link(double mu, int family)
{
    double result = 0.0;
    switch (family)
    {
    case 1: /* gaussian */
        result = mu;
        break;

    case 2: /* binomial */
        result = std::log(mu / (1.0 - mu));
        break;
    }
    return result;
}

inline double dev_unit(double y, double mu, double w, int family)
{
    double result = 0.0;
    switch (family)
    {
    case 1: /* gaussian */
        result = (y - mu) * (y - mu) * w;
        break;

    case 2: /* binomial */
        result = -w * (y == 1.0 ? std::log(mu)
                                : std::log(1.0 - mu));
        break;
    }
    return result;
}

void initialize(double *b0,
                double *b,
                double *dev,
                double *w_sum,
                double *rw,
                double *x_sd,
                double *xwx,
                double *rho,
                const double *y,
                const double *X,
                double *w,
                const double *rho0,
                const int n,
                const int p,
                const int ntune,
                const int standard,
                const int family)
{
    /*
    family: 
    1: gaussian
    2: binomial
    */

    // initialize w and w_sum
    *w_sum = (double)n;
    double weight_sum = accumulate(w, n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        w[i] *= n / weight_sum;
    }

    // initialize b0, b
    // MAY USE BLAS?
    double y_mean = inner_prod(y, w, n, 0.0) / n;
    std::fill(b0, b0 + ntune, link(y_mean, family));
    std::fill(b, b + ntune * p, 0.0);

    // initialize dev
    double null_dev = 0.0;
    for (int i = 0; i < n; ++i)
    {
        rw[i] = (y[i] - y_mean) * w[i];
        null_dev += dev_unit(y[i], y_mean, w[i], family);
    }
    std::fill(dev, dev + ntune, null_dev);

    std::copy(rho0, rho0 + p, rho);

    if (standard)
    {
        for (int j = 0; j < p; ++j)
        {
            double x_mu = inner_prod(X + j * n, w, n, 0.0) / n;
            double wxss = inner_prod(X + j * n, X + j * n, w, n, 0.0) / n;
            double sd = sqrt(wxss - x_mu * x_mu);

            xwx[j] = wxss;
            x_sd[j] = sd;
            rho[j] *= sd;
        }
    }
    else
    {
        for (int j = 0; j < p; ++j)
        {
            xwx[j] = inner_prod(X + j * n, X + j * n, w, n, 0.0) / n;
            x_sd[j] = 1.0;
        }
    }
}
