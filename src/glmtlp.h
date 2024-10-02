/**********
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

#ifndef _GLMTLP_H_
#define _GLMTLP_H_

void linreg_l1_ssr(double *__restrict b0,
                   double *__restrict b,
                   double *__restrict rw,
                   double *__restrict eta,
                   const double *__restrict X,
                   const double w_sum,
                   const double *__restrict xwx,
                   const double *__restrict w,
                   const double *__restrict rho,
                   const double *__restrict lambda,
                   const int nlambda,
                   const int n,
                   const int p,
                   const double delta,
                   const double tol,
                   const int cd_maxit,
                   double *__restrict dev);

void linreg_tlp_ssr(double *__restrict b0,
                    double *__restrict b,
                    double *__restrict rw,
                    double *__restrict eta,
                    const double *__restrict X,
                    const double *__restrict x_sd,
                    const double w_sum,
                    const double *__restrict xwx,
                    const double *__restrict w,
                    const double *__restrict rho,
                    const double *__restrict lambda,
                    const int nlambda,
                    const double tau,
                    const int n,
                    const int p,
                    const double delta,
                    const double tol,
                    const int cd_maxit,
                    const int dc_maxit,
                    double *__restrict dev); 


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
                   double *__restrict dev);

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
                     double *__restrict dev);

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
                      double *__restrict dev);

void logistic_l0_ssr(double *__restrict b0,
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
                     const int *__restrict s,
                     const int ns,
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
                     double *__restrict dev);

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
                const int family);

#endif
