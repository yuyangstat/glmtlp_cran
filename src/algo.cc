/**********
    C++ Routines for Generalized Linear Model Optimization.

    Copyright (C) 2020-2021  Chunlin Li, Yu Yang

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
#include <chrono>
#include <iostream>
#include <random>
#include "glmtlp.h"
#include "utils.h"

/**********
    Coordinate-wise Descent

    See <README.md> for details.

    References: 
    <https://link.springer.com/article/10.1023/A:1017501703105>
    <https://doi.org/10.1214/07-AOAS131>
**********/

void coordinate_descent(double *__restrict b0,
                        double *__restrict b,
                        double *__restrict rw,
                        double *__restrict eta,
                        const double *__restrict X,
                        const double w_sum,
                        const double *__restrict xwx,
                        const double *__restrict w,
                        const double *__restrict rho,
                        const double lambda,
                        const int n,
                        const int p,
                        const double delta,
                        const double tol,
                        const int cd_maxit,
                        int *__restrict it,
                        int *__restrict working_set,
                        const int working_len)
{
    double b_j;
    double b_j_working;
    double b_change;

    int iter = *it;
    // coordinate descent iteration
    for (; iter != cd_maxit; ++iter)
    {
        // absolute difference in solution b
        double difference = 0.0;

        // shuffle?

        for (int idx = 0; idx != working_len; ++idx)
        {
            int j = working_set[idx];
            b_j = b[j];

            // update b_j
            b_j_working = soft_thresh(
                inner_prod(rw, X + j * n, n, 0.0) / (n * xwx[j] * delta) + b_j,
                lambda * rho[j] / (xwx[j] * delta));
            b_change = b_j_working - b_j;

            if (b_change != 0.0)
            {
                difference = std::abs(b_change) > difference
                                 ? std::abs(b_change)
                                 : difference;

                // update r
                if (eta)
                {
                    for (int i = 0; i < n; ++i)
                    {
                        rw[i] -= w[i] * X[j * n + i] * b_change;
                        eta[i] += X[j * n + i] * b_change;
                    }
                }
                else
                {
                    for (int i = 0; i < n; ++i)
                    {
                        rw[i] -= w[i] * X[j * n + i] * b_change;
                    }
                }

                b[j] = b_j_working;
            }
        }

        // update intercept
        b_change = accumulate(rw, n, 0.0) / (w_sum * delta);
        *b0 += b_change;

        for (int i = 0; i < n; ++i)
        {
            rw[i] -= w[i] * b_change;
        }

        if (difference <= tol)
        {
            break;
        }
    }

    *it = iter;
}

void newton_raphson(double *__restrict b0,
                    double *__restrict b,
                    double *__restrict rw,
                    double *__restrict eta,
                    double w_sum,
                    double *__restrict xwx,
                    const double *__restrict y,
                    const double *__restrict X,
                    const double *__restrict w,
                    const double *__restrict rho,
                    const double lambda,
                    const int n,
                    const int p,
                    const double delta,
                    const double tol,
                    int *__restrict it,
                    const int nr_maxit,
                    const int cd_maxit,
                    const int *__restrict is_working,
                    int *__restrict working_set,
                    const int working_len)
{
    // newton-raphson iteration
    int iter = *it;
    int cd_it = 0;

    double *w_working = new double[n]();
    double *b_working = new double[p];
    std::copy(b, b + p, b_working);

    for (; iter != nr_maxit; ++iter)
    {
        // update w, rw, ...
        // this is only for logistic, replace
        for (int i = 0; i < n; ++i)
        {
            if (w[i] != 0.0)
            {
                rw[i] = 1.0 / (1.0 + std::exp(-*b0 - eta[i]));
                if (std::abs(rw[i] - 1.0) < tol)
                {
                    rw[i] = 1.0;
                    w_working[i] = tol;
                }
                else if (std::abs(rw[i]) < tol)
                {
                    rw[i] = 0.0;
                    w_working[i] = tol;
                }
                else
                {
                    w_working[i] = rw[i] * (1.0 - rw[i]);
                }
                rw[i] = (y[i] - rw[i]) * w[i];
            }
        }

        // update w_sum, xwx
        w_sum = accumulate(w_working, n, 0.0);
        for (int j = 0; j < p; ++j)
        {
            if (is_working[j])
            {
                xwx[j] = inner_prod(X + j * n, X + j * n, w_working, n, 0.0)/n;
            }
        }

        // weighted lasso over a working set
        coordinate_descent(b0, b_working, rw, eta, X, w_sum, xwx, w_working, rho,
                           lambda, n, p, delta, tol, cd_maxit,
                           &cd_it, working_set, working_len);

        // check convergence
        double difference = 0.0;
        for (int idx = 0; idx < working_len; ++idx)
        {
            int j = working_set[idx];
            double b_change = b_working[j] - b[j];
            difference = std::abs(b_change) > difference
                             ? std::abs(b_change)
                             : difference;
        }

        std::copy(b_working, b_working + p, b);
        if (difference <= tol)
            break;
    }

    delete[] w_working;
    delete[] b_working;
}
