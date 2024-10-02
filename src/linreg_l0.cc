/**********
    C++ Routine for Computing L0 Regression (dense design version).

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

#include <cmath>    
#include <queue>
#include "algo.h"
#include "glmtlp.h"
#include "utils.h"

void l0_projection(double *__restrict b0,
                   double *__restrict b,
                   double *__restrict b_working,
                   double *__restrict rw_working,
                   double b0_init,
                   const double *__restrict rw_init,
                   const double *__restrict X,
                   const double w_sum,
                   const double *__restrict xwx,
                   const double *__restrict w,
                   const double *__restrict rho,
                   const int *__restrict s,
                   const int ns,
                   const int n,
                   const int p,
                   const double delta,
                   const double tol,
                   const int cd_maxit,
                   double *__restrict dev)
{
    double *eta = NULL;

    int s_max = s[ns - 1]; // assume s is inceasing
    int *supp = new int[p];

    int supp_len = 0;
    for (int j = 0; j < p; ++j)
    {
        if (rho[j] == 0.0) 
        {
            supp[supp_len++] = j;
        }
    }
    
    std::priority_queue<std::pair<double, int>> queue;
    for (int j = 0; j < p; ++j)
    {
        if (std::abs(b_working[j]) > tol && rho[j] != 0.0)
        {
            queue.push(std::pair<double, int>(std::abs(b_working[j]), j));
        }
    }

    // upper bound of size of candidate model
    int df_max = s_max < (int) queue.size() ? s_max : queue.size();

    // sort the coefficients
    for (int df = 0; df < df_max; ++df)
    {
        supp[df + supp_len] = queue.top().second;
        queue.pop();
    }

    std::copy(rw_init, rw_init + n, rw_working);
    std::fill(b_working, b_working + p, 0.0);
    double b0_working = b0_init;

    for (int k = 0; k < ns; ++k)
    {
        int df = s[k];
        if (df > df_max) 
        {
            break;
        }

        int it = 0;
        
        coordinate_descent(&b0_working, b_working, rw_working, eta, X, w_sum, xwx, w, rho, 0.0,
                           n, p, delta, tol, cd_maxit, &it, supp, df + supp_len);

        double dev_working = 0.0;
        for (int i = 0; i < n; ++i)
        {
            if (w[i] != 0.0)
            {
                dev_working += rw_working[i] * rw_working[i] / w[i];
            }
        }

        if (dev_working < dev[k])
        {
            // save the result 
            b0[k] = b0_working;
            std::copy(b_working, b_working + p, b + k * p);
            dev[k] = dev_working;
        }
    }

    delete[] supp;
}


void linreg_l0_ssr(double *__restrict b0,
                   double *__restrict b,
                   double *__restrict rw,
                   double *__restrict eta,
                   const double *__restrict X,
                   const double *__restrict x_sd,
                   const double w_sum,
                   const double *__restrict xwx,
                   const double *__restrict w,
                   const double *__restrict rho,
                   const int *__restrict s,
                   const int ns,
                   const double *__restrict lambda,
                   const int nlambda,
                   const double tau,
                   const int n,
                   const int p,
                   const double delta,
                   const double tol,
                   const int cd_maxit,
                   const int dc_maxit,
                   double *__restrict dev)
{
    /*
        b0:         nlambda-dim intercept array;
        b:          (p * nlambda)-dim coefficient array;
        r:          n-dim working residual array;
        X:          (n * p)-dim design matrix array;
        w:          n-dim observation weight array;
        rho:        p-dim penalty factor array;
        lambda:     penalization parameter array;
        nlambda:    number of candidate lambda;
        n:          number of observations;
        p:          number of features;
        delta:      factor for majorized coordinate descent;
        tol:        tolerance error;
        cd_maxit:   maximum iteration of coordinate descent; 
    */
   
    double *rw_working = new double[n]; // working residual of d.c.
    double *rw_init = new double[n];
    std::copy(rw, rw + n, rw_init);

    // working coordinates
    int *is_working = new int[p];
    int *working_set = new int[p];
    int working_len = 0;

    double *rho_working = new double[p];

    // b_working[0:p-1] saved for warm start
    // b_working[p:2p-1] for working estimating coefficients
    double *b_working = new double[p * 2]();
    double b0_working[2] {b0[0], b0[0]};
    double dev_working[2] {dev[0], dev[0]};

    // save the initial guess of b0
    double b0_init = b0[0];

    for (int k = 1; k < nlambda; ++k)
    {
        // warm start
        std::copy(b_working, b_working+p, b_working + p);
        b0_working[1] = b0_working[0];

        // call lasso solver
        linreg_l1_ssr(b0_working, b_working, rw, eta, X, w_sum, xwx, w, rho,
                      &lambda[k - 1], 2, n, p, delta, tol, cd_maxit, dev_working);

        // save solutions for warm start
        if (k < nlambda - 1)
        {
            b0_working[0] = b0_working[1];
            std::copy(b_working + p, b_working + 2 * p, b_working);
            dev_working[0]=dev_working[1];
        }

        // difference-of-convex program
        std::copy(rw, rw + n, rw_working);
        std::copy(rho, rho + p, rho_working);

        int it = 0;

        for (; it < dc_maxit; ++it)
        {
            int converge = 1; // indicate whether rho changes

            // update rho
            for (int j = 0; j < p; ++j)
            {
                if (std::abs(b_working[p + j]) * x_sd[j] >= tau)
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
                if (b_working[p + j] != 0.0)
                {
                    working_set[working_len++] = j;
                    is_working[j] = 1;
                }
            }

            // iterate on working set
            int it_cd = 0;
            for (;;)
            {
                coordinate_descent(&b0_working[1], b_working + p, rw_working, eta, X, w_sum, xwx, w, rho_working,
                                   lambda[k], n, p, delta, tol, cd_maxit,
                                   &it_cd, working_set, working_len);

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

                if (kkt || it_cd >= cd_maxit)
                    break;
            }
        }

        // project to { b : ||b||_0 <= s }
        // iterate from r_init and b = 0

        // NEED TO OPTIMIZE THIS!
        l0_projection(b0, b, b_working + p, rw_working, b0_init, rw_init,
                      X, w_sum, xwx, w, rho, s, ns, n, p, delta, tol, cd_maxit, dev);
    }

    delete[] rw_working;
    delete[] rw_init;
    delete[] is_working;
    delete[] working_set;
    delete[] rho_working;
    delete[] b_working;
}

