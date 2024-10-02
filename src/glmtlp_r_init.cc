/**********
    R Interface.

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

#define R_NO_REMAP

#include "glmtlp.h"
#include "R.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"

int Rf_nrows(SEXP s);
int Rf_ncols(SEXP s);
int LENGTH(SEXP x);
#define nrows Rf_nrows
#define ncols Rf_ncols
R_len_t Rf_length(SEXP x);
#define length Rf_length

extern "C"
{
    SEXP gaussian_l1(SEXP y_vec,
                     SEXP X_mat,
                     SEXP weights,
                     SEXP penalty_factor,
                     SEXP lambda_vec,
                     SEXP delta_val,
                     SEXP standardize,
                     SEXP tol,
                     SEXP cd_maxit)
    {
        // input
        double *y = REAL(y_vec);
        double *X = REAL(X_mat);
        double *w = REAL(weights);
        double *rho0 = REAL(penalty_factor);
        double *lambda = REAL(lambda_vec);
        double delta = *REAL(delta_val);
        int standard = *INTEGER(standardize);

        int n = nrows(X_mat);
        int p = ncols(X_mat);
        int nlambda = length(lambda_vec);

        // output
        SEXP intercept = PROTECT(Rf_allocVector(REALSXP, nlambda));
        SEXP beta = PROTECT(Rf_allocVector(REALSXP, p * nlambda));
        SEXP deviance = PROTECT(Rf_allocVector(REALSXP, nlambda));

        double *b0 = REAL(intercept);
        double *b = REAL(beta);
        double *dev = REAL(deviance);

        // used for computing, need initialization
        double w_sum;
        double *rw = new double[n];
        double *x_sd = new double[p];
        double *xwx = new double[p];
        double *rho = new double[p];
        double *eta = NULL;

        initialize(b0, b, dev, &w_sum, rw, x_sd, xwx, rho, y, X, w, rho0, n, p, nlambda, standard, 1);

        // call lasso solver
        linreg_l1_ssr(b0,
                      b,
                      rw,
                      eta,
                      X,
                      w_sum,
                      xwx,
                      w,
                      rho,
                      lambda,
                      nlambda,
                      n,
                      p,
                      delta,
                      *REAL(tol),
                      *INTEGER(cd_maxit),
                      dev);

        delete[] x_sd;
        delete[] xwx;
        delete[] rho;
        delete[] rw;

        SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));

        SET_VECTOR_ELT(result, 0, intercept);
        SET_VECTOR_ELT(result, 1, beta);
        SET_VECTOR_ELT(result, 2, deviance);

        UNPROTECT(4);
        return result;
    }

    SEXP gaussian_tlp(SEXP y_vec,
                      SEXP X_mat,
                      SEXP weights,
                      SEXP penalty_factor,
                      SEXP lambda_vec,
                      SEXP tau,
                      SEXP delta_val,
                      SEXP standardize,
                      SEXP tol,
                      SEXP cd_maxit,
                      SEXP dc_maxit)
    {
        // input
        double *y = REAL(y_vec);
        double *X = REAL(X_mat);
        double *w = REAL(weights);
        double *rho0 = REAL(penalty_factor);
        double *lambda = REAL(lambda_vec);
        double delta = *REAL(delta_val);
        int standard = *INTEGER(standardize);

        int n = nrows(X_mat);
        int p = ncols(X_mat);
        int nlambda = length(lambda_vec);

        // output
        SEXP intercept = PROTECT(Rf_allocVector(REALSXP, nlambda));
        SEXP beta = PROTECT(Rf_allocVector(REALSXP, p * nlambda));
        SEXP deviance = PROTECT(Rf_allocVector(REALSXP, nlambda));

        double *b0 = REAL(intercept);
        double *b = REAL(beta);
        double *dev = REAL(deviance);

        // used for computing, need initialization
        double w_sum;
        double *rw = new double[n];
        double *x_sd = new double[p];
        double *xwx = new double[p];
        double *rho = new double[p];
        double *eta = NULL;

        initialize(b0, b, dev, &w_sum, rw, x_sd, xwx, rho, y, X, w, rho0, n, p, nlambda, standard, 1);

        // call tlp solver
        linreg_tlp_ssr(b0,
                       b,
                       rw,
                       eta,
                       X,
                       x_sd,
                       w_sum,
                       xwx,
                       w,
                       rho,
                       lambda,
                       nlambda,
                       *REAL(tau),
                       n,
                       p,
                       delta,
                       *REAL(tol),
                       *INTEGER(cd_maxit),
                       *INTEGER(dc_maxit),
                       dev);

        delete[] x_sd;
        delete[] xwx;
        delete[] rho;
        delete[] rw;

        SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
        SET_VECTOR_ELT(result, 0, intercept);
        SET_VECTOR_ELT(result, 1, beta);
        SET_VECTOR_ELT(result, 2, deviance);

        UNPROTECT(4);
        return result;
    }

    SEXP gaussian_l0(SEXP y_vec,
                     SEXP X_mat,
                     SEXP weights,
                     SEXP penalty_factor,
                     SEXP kappa_vec,
                     SEXP lambda_vec,
                     SEXP tau,
                     SEXP delta_val,
                     SEXP standardize,
                     SEXP tol,
                     SEXP cd_maxit,
                     SEXP dc_maxit)
    {
        // input
        double *y = REAL(y_vec);
        double *X = REAL(X_mat);
        double *w = REAL(weights);
        double *rho0 = REAL(penalty_factor);
        double *lambda = REAL(lambda_vec);
        double delta = *REAL(delta_val);
        int *kappa = INTEGER(kappa_vec);
        int standard = *INTEGER(standardize);

        int n = nrows(X_mat);
        int p = ncols(X_mat);
        int nkappa = length(kappa_vec);
        int nlambda = length(lambda_vec);

        // output
        SEXP intercept = PROTECT(Rf_allocVector(REALSXP, nkappa));
        SEXP beta = PROTECT(Rf_allocVector(REALSXP, p * nkappa));
        SEXP deviance = PROTECT(Rf_allocVector(REALSXP, nkappa));

        double *b0 = REAL(intercept);
        double *b = REAL(beta);
        double *dev = REAL(deviance);

        // used for computing, need initialization
        double w_sum;
        double *rw = new double[n];
        double *x_sd = new double[p];
        double *xwx = new double[p];
        double *rho = new double[p];
        double *eta = NULL;

        initialize(b0, b, dev, &w_sum, rw, x_sd, xwx, rho, y, X, w, rho0, n, p, nkappa, standard, 1);

        linreg_l0_ssr(b0,
                      b,
                      rw,
                      eta,
                      X,
                      x_sd,
                      w_sum,
                      xwx,
                      w,
                      rho,
                      kappa,
                      nkappa,
                      lambda,
                      nlambda,
                      *REAL(tau),
                      n,
                      p,
                      delta,
                      *REAL(tol),
                      *INTEGER(cd_maxit),
                      *INTEGER(dc_maxit),
                      dev);

        delete[] x_sd;
        delete[] xwx;
        delete[] rho;
        delete[] rw;

        SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));

        SET_VECTOR_ELT(result, 0, intercept);
        SET_VECTOR_ELT(result, 1, beta);
        SET_VECTOR_ELT(result, 2, deviance);

        UNPROTECT(4);
        return result;
    }

    SEXP logistic_l1(SEXP y_vec,
                     SEXP X_mat,
                     SEXP weights,
                     SEXP penalty_factor,
                     SEXP lambda_vec,
                     SEXP delta_val,
                     SEXP standardize,
                     SEXP tol,
                     SEXP nr_maxit,
                     SEXP cd_maxit)
    {
        // input
        double *y = REAL(y_vec);
        double *X = REAL(X_mat);
        double *w = REAL(weights);
        double *rho0 = REAL(penalty_factor);
        double *lambda = REAL(lambda_vec);
        double delta = *REAL(delta_val);
        int standard = *INTEGER(standardize);

        int n = nrows(X_mat);
        int p = ncols(X_mat);
        int nlambda = length(lambda_vec);

        // output
        SEXP intercept = PROTECT(Rf_allocVector(REALSXP, nlambda));
        SEXP beta = PROTECT(Rf_allocVector(REALSXP, p * nlambda));
        SEXP deviance = PROTECT(Rf_allocVector(REALSXP, nlambda));

        double *b0 = REAL(intercept);
        double *b = REAL(beta);
        double *dev = REAL(deviance);

        // used for computing, need initialization
        double w_sum;
        double *rw = new double[n];
        double *x_sd = new double[p];
        double *xwx = new double[p];
        double *rho = new double[p];
        double *eta = new double[n]();

        initialize(b0, b, dev, &w_sum, rw, x_sd, xwx, rho, y, X, w, rho0, n, p, nlambda, standard, 2);

        logistic_l1_ssr(b0,
                        b,
                        rw,
                        eta,
                        y,
                        X,
                        w_sum,
                        xwx,
                        w,
                        rho,
                        lambda,
                        nlambda,
                        n,
                        p,
                        delta,
                        *REAL(tol),
                        *INTEGER(nr_maxit),
                        *INTEGER(cd_maxit),
                        dev);

        delete[] x_sd;
        delete[] xwx;
        delete[] eta;
        delete[] rho;
        delete[] rw;

        SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));

        SET_VECTOR_ELT(result, 0, intercept);
        SET_VECTOR_ELT(result, 1, beta);
        SET_VECTOR_ELT(result, 2, deviance);

        UNPROTECT(4);
        return result;
    }

    SEXP logistic_tlp(SEXP y_vec,
                      SEXP X_mat,
                      SEXP weights,
                      SEXP penalty_factor,
                      SEXP lambda_vec,
                      SEXP tau,
                      SEXP delta_val,
                      SEXP standardize,
                      SEXP tol,
                      SEXP nr_maxit,
                      SEXP cd_maxit,
                      SEXP dc_maxit)
    {
        // input
        double *y = REAL(y_vec);
        double *X = REAL(X_mat);
        double *w = REAL(weights);
        double *rho0 = REAL(penalty_factor);
        double *lambda = REAL(lambda_vec);
        double delta = *REAL(delta_val);
        int standard = *INTEGER(standardize);

        int n = nrows(X_mat);
        int p = ncols(X_mat);
        int nlambda = length(lambda_vec);

        // output
        SEXP intercept = PROTECT(Rf_allocVector(REALSXP, nlambda));
        SEXP beta = PROTECT(Rf_allocVector(REALSXP, p * nlambda));
        SEXP deviance = PROTECT(Rf_allocVector(REALSXP, nlambda));

        double *b0 = REAL(intercept);
        double *b = REAL(beta);
        double *dev = REAL(deviance);

        // used for computing, need initialization
        double w_sum;
        double *rw = new double[n];
        double *x_sd = new double[p];
        double *xwx = new double[p];
        double *rho = new double[p];
        double *eta = new double[n]();

        initialize(b0, b, dev, &w_sum, rw, x_sd, xwx, rho, y, X, w, rho0, n, p, nlambda, standard, 2);

        logistic_tlp_ssr(b0,
                         b,
                         rw,
                         eta,
                         y,
                         X,
                         x_sd,
                         w_sum,
                         xwx,
                         w,
                         rho,
                         lambda,
                         nlambda,
                         *REAL(tau),
                         n,
                         p,
                         delta,
                         *REAL(tol),
                         *INTEGER(nr_maxit),
                         *INTEGER(cd_maxit),
                         *INTEGER(dc_maxit),
                         dev);

        delete[] x_sd;
        delete[] xwx;
        delete[] eta;
        delete[] rho;
        delete[] rw;

        SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));

        SET_VECTOR_ELT(result, 0, intercept);
        SET_VECTOR_ELT(result, 1, beta);
        SET_VECTOR_ELT(result, 2, deviance);

        UNPROTECT(4);
        return result;
    }

    SEXP logistic_l0(SEXP y_vec,
                     SEXP X_mat,
                     SEXP weights,
                     SEXP penalty_factor,
                     SEXP kappa_vec,
                     SEXP lambda_vec,
                     SEXP tau,
                     SEXP delta_val,
                     SEXP standardize,
                     SEXP tol,
                     SEXP nr_maxit,
                     SEXP cd_maxit,
                     SEXP dc_maxit)
    {
        // input
        double *y = REAL(y_vec);
        double *X = REAL(X_mat);
        double *w = REAL(weights);
        double *rho0 = REAL(penalty_factor);
        double *lambda = REAL(lambda_vec);
        double delta = *REAL(delta_val);
        int *kappa = INTEGER(kappa_vec);
        int standard = *INTEGER(standardize);

        int n = nrows(X_mat);
        int p = ncols(X_mat);
        int nkappa = length(kappa_vec);
        int nlambda = length(lambda_vec);

        // output
        SEXP intercept = PROTECT(Rf_allocVector(REALSXP, nkappa));
        SEXP beta = PROTECT(Rf_allocVector(REALSXP, p * nkappa));
        SEXP deviance = PROTECT(Rf_allocVector(REALSXP, nkappa));

        double *b0 = REAL(intercept);
        double *b = REAL(beta);
        double *dev = REAL(deviance);

        // used for computing, need initialization
        double w_sum;
        double *rw = new double[n];
        double *x_sd = new double[p];
        double *xwx = new double[p];
        double *rho = new double[p];
        double *eta = new double[n]();

        initialize(b0, b, dev, &w_sum, rw, x_sd, xwx, rho, y, X, w, rho0, n, p, nkappa, standard, 2);

        logistic_l0_ssr(b0,
                        b,
                        rw,
                        eta,
                        y,
                        X,
                        x_sd,
                        w_sum,
                        xwx,
                        w,
                        rho,
                        kappa,
                        nkappa,
                        lambda,
                        nlambda,
                        *REAL(tau),
                        n,
                        p,
                        delta,
                        *REAL(tol),
                        *INTEGER(nr_maxit),
                        *INTEGER(cd_maxit),
                        *INTEGER(dc_maxit),
                        dev);

        delete[] x_sd;
        delete[] xwx;
        delete[] eta;
        delete[] rho;
        delete[] rw;

        SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));

        SET_VECTOR_ELT(result, 0, intercept);
        SET_VECTOR_ELT(result, 1, beta);
        SET_VECTOR_ELT(result, 2, deviance);

        UNPROTECT(4);
        return result;
    }

} // extern "C"

static const R_CallMethodDef CallEntries[] = {
    {"gaussian_l0", (DL_FUNC)&gaussian_l0, 12},
    {"gaussian_l1", (DL_FUNC)&gaussian_l1, 9},
    {"gaussian_tlp", (DL_FUNC)&gaussian_tlp, 11},
    {"logistic_l0", (DL_FUNC)&logistic_l0, 13},
    {"logistic_l1", (DL_FUNC)&logistic_l1, 10},
    {"logistic_tlp", (DL_FUNC)&logistic_tlp, 12},
    {NULL, NULL, 0}};

void R_init_glmtlp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
