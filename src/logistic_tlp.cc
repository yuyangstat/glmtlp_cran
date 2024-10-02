/**********
    C++ Routine for Computing TLP Logistic Regression (dense design version).

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

#include <cmath>    // exp, abs
#include <iostream> // printf
#include "algo.h"
#include "glmtlp.h"
#include "utils.h"

void logistic_tlp_ssr(double *__restrict b0,
                      double *__restrict b,
                      double *__restrict rw,
                      double *__restrict eta,
                      const double *__restrict y,
                      const double *__restrict X,
                      const double *__restrict x_sd,
                      double w_sum,
                      double *__restrict xwx,
                      double *__restrict w,
                      const double *__restrict rho,
                      const double *__restrict lambda,
                      const int nlambda,
                      const double tau,
                      const int n,
                      const int p,
                      const double delta,
                      const double tol,
                      const int nr_maxit,
                      const int cd_maxit,
                      const int dc_maxit,
                      double *__restrict dev)
{
    double *rw_working = new double[n]; // working residual of d.c.
    double *eta_working = new double[n];

    // working coordinates
    int *is_working = new int[p]();
    int *working_set = new int[p];
    int working_len = 0;

    double *rho_working = new double[p];

    for (int k = 1; k < nlambda; ++k)
    {

        // call l1 logistic solver
        logistic_l1_ssr(&b0[k - 1], &b[(k - 1) * p], rw, eta, y, X, w_sum, xwx, w, rho,
                        &lambda[k - 1], 2, n, p, delta, tol, nr_maxit, cd_maxit, &dev[k - 1]);

        if (k != nlambda - 1)
        {
            b0[k + 1] = b0[k];
            std::copy(b + k * p, b + k * p + p, b + k * p + p);
        }

        // difference-of-convex program
        std::copy(rw, rw + n, rw_working);
        std::copy(eta, eta + n, eta_working);
        std::copy(rho, rho + p, rho_working);

        int it = 0;

        for (; it < dc_maxit; ++it)
        {
            int converge = 1; // indicate whether rho changes

            // update rho, tlp module
            for (int j = 0; j < p; ++j)
            {
                if (std::abs(b[k * p + j]) * x_sd[j] >= tau)
                {
                    converge = rho_working[j] == 0.0 ? converge : 0;
                    rho_working[j] = 0.0;
                }
                else
                {
                    converge = rho_working[j] == rho[j] ? converge : 0;
                    rho_working[j] = rho[j];
                }
            }

            // check convergence
            if (converge)
                break;

            // working set
            std::fill(is_working, is_working + p, 0);
            working_len = 0;
            for (int j = 0; j < p; ++j)
            {
                if (b[k * p + j] != 0.0)
                {
                    working_set[working_len++] = j;
                    is_working[j] = 1;
                }
            }

            // iterate on working set
            int it_nr = 0;
            for (;;)
            {
                // newton
                newton_raphson(b0 + k, b + k * p, rw_working, eta_working, w_sum, xwx, y, X, w, rho_working,
                               lambda[k], n, p, delta, tol, &it_nr, nr_maxit, cd_maxit, is_working,
                               working_set, working_len);

                // check kkt for all
                int kkt = 1;
                for (int j = 0; j < p; ++j)
                {
                    if (!is_working[j] &&
                        std::abs(inner_prod(rw_working, X + j * n, n, 0.0)) / n >
                            lambda[k] * rho_working[j])
                    {
                        // KKT condition is violated
                        // add variable to working set
                        kkt = 0;
                        working_set[working_len++] = j;
                        is_working[j] = 1;
                    }
                }

                if (kkt || it_nr >= nr_maxit)
                    break;
            }
        }

        // update deviance
        double deviance = 0.0;
        for (int i = 0; i < n; ++i)
        {
            if (w[i] != 0.0)
            {
                deviance -= w[i] * (y[i] == 1.0 ? std::log(1.0 - rw_working[i])
                                                : std::log(1.0 + rw_working[i]));
            }
        }
        
        dev[k] = deviance;
        
        if (deviance < 0.01 * dev[0] || std::isnan(deviance))
        {
            if (k != nlambda - 1)
            {
                std::fill(dev + k + 1, dev + nlambda, NAN);
                std::fill(b0 + k + 1, b0 + nlambda, NAN);
                std::fill(b + k * p + p, b + nlambda * p, NAN);
            }

            // for (; k < nlambda; ++k)
            // {
            //     dev[k] = deviance;
            // }
            break;
        }

        // IF DIVERGES, WARNING?

        //if (it >= dc_maxit)
        //{
        //printf("it = %d, lambda = %f.\n", it, lambda[k]);
        // printf(
        //     "Warning: the coordinate descent algorithm does not converge (lambda = %f).\n",
        //     lambda[k]);
        //}
    }

    // free spaces
    delete[] rw_working;
    delete[] eta_working;
    delete[] is_working;
    delete[] working_set;
    delete[] rho_working;
}
