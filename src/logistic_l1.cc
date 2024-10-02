/**********
    C++ Routine for Computing L1 Logistic Regression (dense design version).

    Copyright (C) 2021  Chunlin Li

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

/**********
    Minimize: 
    -sum_i [y_i log(mu_i) + (1 - y_i)log(1 - mu_i)]/n + lambda sum_j [rho_j |b_j|]

    Algorithm:
    ?????
    

    See <README.md> for details.
**********/

#include <cmath>    // exp, abs
#include <iostream> // printf
#include "algo.h"
#include "glmtlp.h"
#include "utils.h"

// logistic regression with lasso penalty: no observation weight

void logistic_l1_ssr(double *__restrict b0,
                     double *__restrict b,
                     double *__restrict rw,
                     double *__restrict eta,
                     const double *__restrict y,
                     const double *__restrict X,
                     double w_sum,
                     double *__restrict xwx,
                     double *__restrict w,
                     const double *__restrict rho,
                     const double *__restrict lambda,
                     const int nlambda,
                     const int n,
                     const int p,
                     const double delta,
                     const double tol,
                     const int nr_maxit,
                     const int cd_maxit,
                     double *__restrict dev)
{
    // active set, strong set, and working set
    int *ever_active = new int[p]();
    int *is_strong = new int[p]();
    int *is_working = new int[p]();

    // working coordinates
    int *working_set = new int[p];
    int working_len = 0;

    for (int k = 1; k != nlambda; ++k)
    {
        /* 
            strong rule: 
            |sum_i [x_{ij} r_i(lambda_k)]| < 2lambda_{k+1} - lambda_k
            where r(lambda_k) satisfies optimality condition 
        */
        // NEED TO HANDLE EXCEPTIONS (DIVERGENCE CASE)
        for (int j = 0; j < p; ++j)
        {
            // IF RHO > 0.
            if (std::abs(inner_prod(rw, X + j * n, n, 0.0)) / n >=
                (2.0 * lambda[k] - lambda[k - 1]) * rho[j])
            {
                is_strong[j] = 1;
            }
            else
            {
                is_strong[j] = 0;
            }
        }

        // initialize working set
        working_len = 0;
        if (k > 1)
        {
            for (int j = 0; j < p; ++j)
            {
                if (ever_active[j] || b[k * p + j] != 0.0)
                {
                    ever_active[j] = 1;
                    is_working[j] = 1;
                    working_set[working_len++] = j;
                }
                else
                {
                    is_working[j] = 0;
                }
            }
        }
        else
        {
            for (int j = 0; j < p; ++j)
            {
                if (is_strong[j])
                {
                    working_set[working_len++] = j;
                    is_working[j] = 1;
                }
            }
        }

        int it = 0;  // newton-raphson iteration number
        int kkt = 0; // kkt condition
        for (;;)
        {

            for (;;)
            {

                // newton-raphson over strong set
                newton_raphson(b0 + k, b + k * p, rw, eta, w_sum, xwx, y, X, w, rho,
                               lambda[k], n, p, delta, tol, &it, nr_maxit,
                               cd_maxit, is_working, working_set, working_len);

                // check KKT, update working set
                kkt = 1;
                for (int j = 0; j < p; ++j)
                {
                    if (is_strong[j] && !is_working[j] &&
                        std::abs(inner_prod(rw, X + j * n, n, 0.0)) / n >
                            lambda[k] * rho[j])
                    {
                        // KKT condition is violated
                        // add variable to working set
                        kkt = 0;
                        working_set[working_len++] = j;
                        is_working[j] = 1;
                    }
                }

                // termination
                if (kkt || it >= nr_maxit)
                    break;
            }

            if (kkt)
            {
                for (int j = 0; j < p; ++j)
                {
                    if (!is_working[j] && !is_strong[j] &&
                        std::abs(inner_prod(rw, X + j * n, n, 0.0)) / n >
                            lambda[k] * rho[j])
                    {
                        // KKT condition is violated
                        // add variable to working set
                        kkt = 0;
                        working_set[working_len++] = j;
                        is_working[j] = 1;
                    }
                }
            }

            // termination
            if (kkt || it >= nr_maxit)
                break;
        }

        if (kkt)
        {

            // deviance
            double deviance = 0.0;
            for (int i = 0; i < n; ++i)
            {
                if (w[i] != 0.0)
                {
                    deviance -= y[i] == 1.0 ? w[i] * std::log(1.0 - rw[i])
                                            : w[i] * std::log(1.0 + rw[i]);
                }
            }

            dev[k] = deviance;

            if (k != nlambda - 1)
            {
                if (deviance < 0.01 * dev[0])
                {
                    std::fill(dev + k + 1, dev + nlambda, NAN);
                    std::fill(b0 + k + 1, b0 + nlambda, NAN);
                    std::fill(b + k * p + p, b + nlambda * p, NAN);
                    break;
                }

                // estimate converges at lambda[k]
                // warm start
                std::copy(b + k * p, b + k * p + p, b + k * p + p);
                b0[k + 1] = b0[k];
            }
        }
        else
        {
            // NEED TO HANDLE DIVERGENT CASE (START FROM LAST KKT SOLUTION)

            //printf(
            //    "Warning: the coordinate descent algorithm does not converge (lambda = %f).\n",
            //    lambda[k]);
        }
    }

    // free spaces
    delete[] ever_active;
    delete[] is_strong;
    delete[] is_working;
    delete[] working_set;
}
