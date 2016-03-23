#pragma once

#include "utils.hpp"

#include <algorithm>

namespace bobyqa_detail {

template <class Function>
void prelim(
    const Function &function,
    const long n,
    const long npt,
    double *const x,
    const double *const xl,
    const double *const xu,
    const double rhobeg,
    const long maxfun,
    double *const xbase,
    double *xpt,
    double *const fval,
    double *const gopt,
    double *const hq,
    double *const pq,
    double *bmat,
    double *zmat,
    const long ndim,
    const double *const sl,
    const double *const su,
    long& nf,
    long& kopt
)
{
    /* Local variables */
    long ih, nfm;
    long nfx = 0, ipt = 0, jpt = 0;
    double fbeg = 0, diff = 0, temp = 0, stepa = 0, stepb = 0;
    long itemp;

    /*     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the */
    /*       same as the corresponding arguments in SUBROUTINE BOBYQA. */
    /*     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU */
    /*       are the same as the corresponding arguments in BOBYQB, the elements */
    /*       of SL and SU being set in BOBYQA. */
    /*     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but */
    /*       it is set by PRELIM to the gradient of the quadratic model at XBASE. */
    /*       If XOPT is nonzero, BOBYQB will change it to its usual value later. */
    /*     NF is maintaned as the number of calls of CALFUN so far. */
    /*     KOPT will be such that the least calculated value of F so far is at */
    /*       the point XPT(KOPT,.)+XBASE in the space of the variables. */

    /*     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, */
    /*     BMAT and ZMAT for the first iteration, and it maintains the values of */
    /*     NF and KOPT. The vector X is also changed by PRELIM. */

    /*     Set some constants. */

    /* Parameter adjustments */
    const long zmat_dim1 = npt;
    const long zmat_offset = 1 + zmat_dim1;
    zmat -= zmat_offset;
    const long xpt_dim1 = npt;
    const long xpt_offset = 1 + xpt_dim1;
    xpt -= xpt_offset;
    const long bmat_dim1 = ndim;
    const long bmat_offset = 1 + bmat_dim1;
    bmat -= bmat_offset;

    /* Function Body */
    const double rhosq = rhobeg * rhobeg;
    const double recip = 1.0 / rhosq;
    const long np = n + 1;

    /*     Set XBASE to the initial vector of variables, and set the initial */
    /*     elements of XPT, BMAT, HQ, PQ and ZMAT to zero. */

    for (long j = 1; j <= n; ++j) {
        xbase[j] = x[j];
        for (long k = 1; k <= npt; ++k) {
            xpt[k + j * xpt_dim1] = 0.0;
        }
        for (long i = 1; i <= ndim; ++i) {
            bmat[i + j * bmat_dim1] = 0.0;
        }
    }
    const long ih_n = n * np / 2;
    for (long ih = 1; ih <= ih_n; ++ih) {
        hq[ih] = 0.0;
    }
    for (long k = 1; k <= npt; ++k) {
        pq[k] = 0.0;
        const long j_n = npt - np;
        for (long j = 1; j <= j_n; ++j) {
            zmat[k + j * zmat_dim1] = 0.0;
        }
    }

    /*     Begin the initialization procedure. NF becomes one more than the number */
    /*     of function values so far. The coordinates of the displacement of the */
    /*     next initial interpolation point from XBASE are set in XPT(NF+1,.). */

    nf = 0;
L50:
    nfm = nf;
    nfx = nf - n;
    ++(nf);
    if (nfm <= n << 1) {
        if (nfm >= 1 && nfm <= n) {
            stepa = rhobeg;
            if (su[nfm] == 0.0) {
                stepa = -stepa;
            }
            xpt[nf + nfm * xpt_dim1] = stepa;
        } else if (nfm > n) {
            stepa = xpt[nf - n + nfx * xpt_dim1];
            stepb = -(rhobeg);
            if (sl[nfx] == 0.0) {
                stepb = std::min(2.0 * rhobeg, su[nfx]);
            }
            if (su[nfx] == 0.0) {
                stepb = std::max(-2.0 * rhobeg, sl[nfx]);
            }
            xpt[nf + nfx * xpt_dim1] = stepb;
        }
    } else {
        itemp = (nfm - np) / n;
        jpt = nfm - itemp * n - n;
        ipt = jpt + itemp;
        if (ipt > n) {
            itemp = jpt;
            jpt = ipt - n;
            ipt = itemp;
        }
        xpt[nf + ipt * xpt_dim1] = xpt[ipt + 1 + ipt * xpt_dim1];
        xpt[nf + jpt * xpt_dim1] = xpt[jpt + 1 + jpt * xpt_dim1];
    }

    /*     Calculate the next value of F. The least function value so far and */
    /*     its index are required. */

    for (long j = 1; j <= n; ++j) {
        x[j] = std::min(std::max(xl[j], xbase[j] + xpt[nf + j * xpt_dim1]), xu[j]);
        if (xpt[nf + j * xpt_dim1] == sl[j]) {
            x[j] = xl[j];
        }
        if (xpt[nf + j * xpt_dim1] == su[j]) {
            x[j] = xu[j];
        }
    }
    const double f = function(n, x + 1);
    fval[nf] = f;
    if (nf == 1) {
        fbeg = f;
        kopt = 1;
    } else if (f < fval[kopt]) {
        kopt = nf;
    }

    /*     Set the nonzero initial elements of BMAT and the quadratic model in the */
    /*     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions */
    /*     of the NF-th and (NF-N)-th interpolation points may be switched, in */
    /*     order that the function value at the first of them contributes to the */
    /*     off-diagonal second derivative terms of the initial quadratic model. */

    if (nf <= (n << 1) + 1) {
        if (nf >= 2 && nf <= n + 1) {
            gopt[nfm] = (f - fbeg) / stepa;
            if (npt < nf + n) {
                bmat[nfm * bmat_dim1 + 1] = -1.0 / stepa;
                bmat[nf + nfm * bmat_dim1] = 1.0 / stepa;
                bmat[npt + nfm + nfm * bmat_dim1] = -0.5 * rhosq;
            }
        } else if (nf >= n + 2) {
            ih = nfx * (nfx + 1) / 2;
            temp = (f - fbeg) / stepb;
            diff = stepb - stepa;
            hq[ih] = 2.0 * (temp - gopt[nfx]) / diff;
            gopt[nfx] = (gopt[nfx] * stepb - temp * stepa) / diff;
            if (stepa * stepb < 0.0) {
                if (f < fval[nf - n]) {
                    fval[nf] = fval[nf - n];
                    fval[nf - n] = f;
                    if (kopt == nf) {
                        kopt = nf - n;
                    }
                    xpt[nf - n + nfx * xpt_dim1] = stepb;
                    xpt[nf + nfx * xpt_dim1] = stepa;
                }
            }
            bmat[nfx * bmat_dim1 + 1] = -(stepa + stepb) / (stepa * stepb);
            bmat[nf + nfx * bmat_dim1] = -0.5 / xpt[nf - n + nfx * xpt_dim1];
            bmat[nf - n + nfx * bmat_dim1] = -bmat[nfx * bmat_dim1 + 1] - bmat[nf + nfx * bmat_dim1];
            zmat[nfx * zmat_dim1 + 1] = sqrt_2 / (stepa * stepb);
            zmat[nf + nfx * zmat_dim1] = sqrt_0_5 / rhosq;
            zmat[nf - n + nfx * zmat_dim1] = -zmat[nfx * zmat_dim1 + 1] - zmat[nf + nfx * zmat_dim1];
        }

        /*     Set the off-diagonal second derivatives of the Lagrange functions and */
        /*     the initial quadratic model. */

    } else {
        ih = ipt * (ipt - 1) / 2 + jpt;
        zmat[nfx * zmat_dim1 + 1] = recip;
        zmat[nf + nfx * zmat_dim1] = recip;
        zmat[ipt + 1 + nfx * zmat_dim1] = -recip;
        zmat[jpt + 1 + nfx * zmat_dim1] = -recip;
        temp = xpt[nf + ipt * xpt_dim1] * xpt[nf + jpt * xpt_dim1];
        hq[ih] = (fbeg - fval[ipt + 1] - fval[jpt + 1] + f) / temp;
    }
    if (nf < npt && nf < maxfun) {
        goto L50;
    }

}

} // namespace bobyqa_detail
