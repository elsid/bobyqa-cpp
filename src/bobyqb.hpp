#pragma once

#include "altmov.hpp"
#include "prelim.hpp"
#include "rescue.hpp"
#include "trsbox.hpp"

namespace bobyqa_detail {

template <class Function>
double bobyqb(
    const Function &function,
    const long n,
    const long npt,
    double *const x,
    const double *const xl,
    const double *const xu,
    const double rhobeg,
    const double rhoend,
    const long maxfun,
    double *const xbase,
    double *xpt,
    double *const fval,
    double *const xopt,
    double *const gopt,
    double *const hq,
    double *const pq,
    double *bmat,
    double *zmat,
    const long ndim,
    double *const sl,
    double *const su,
    double *const xnew,
    double *const xalt,
    double *const d__,
    double *const vlag,
    double *const w
)
{
    /* Local variables */
    double f = 0;
    long ih, nf, jp;
    double dx;
    double den = 0, dsq = 0, rho = 0, sum = 0, diff = 0, beta = 0, gisq = 0;
    long knew = 0;
    double temp, suma, sumb, bsum, fopt;
    long kopt = 0;
    double curv;
    long ksav;
    double gqsq = 0, dist = 0, sumw = 0, sumz = 0, diffa = 0, diffb = 0, diffc = 0, hdiag = 0;
    long kbase;
    double alpha = 0, delta = 0, adelt = 0, denom = 0, fsave = 0, bdtol = 0, delsq = 0;
    long nresc, nfsav;
    double ratio = 0, dnorm = 0, vquad = 0, pqold = 0;
    long itest;
    double sumpq, scaden;
    double errbig, cauchy = 0, fracsq, biglsq, densav;
    double bdtest;
    double crvmin, frhosq;
    double distsq;
    long ntrits;
    double xoptsq;



    /*     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN */
    /*       are identical to the corresponding arguments in SUBROUTINE BOBYQA. */
    /*     XBASE holds a shift of origin that should reduce the contributions */
    /*       from rounding errors to values of the model and Lagrange functions. */
    /*     XPT is a two-dimensional array that holds the coordinates of the */
    /*       interpolation points relative to XBASE. */
    /*     FVAL holds the values of F at the interpolation points. */
    /*     XOPT is set to the displacement from XBASE of the trust region centre. */
    /*     GOPT holds the gradient of the quadratic model at XBASE+XOPT. */
    /*     HQ holds the explicit second derivatives of the quadratic model. */
    /*     PQ contains the parameters of the implicit second derivatives of the */
    /*       quadratic model. */
    /*     BMAT holds the last N columns of H. */
    /*     ZMAT holds the factorization of the leading NPT by NPT submatrix of H, */
    /*       this factorization being ZMAT times ZMAT^T, which provides both the */
    /*       correct rank and positive semi-definiteness. */
    /*     NDIM is the first dimension of BMAT and has the value NPT+N. */
    /*     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively. */
    /*       All the components of every XOPT are going to satisfy the bounds */
    /*       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when */
    /*       XOPT is on a constraint boundary. */
    /*     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the */
    /*       vector of variables for the next call of CALFUN. XNEW also satisfies */
    /*       the SL and SU constraints in the way that has just been mentioned. */
    /*     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW */
    /*       in order to increase the denominator in the updating of UPDATE. */
    /*     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT. */
    /*     VLAG contains the values of the Lagrange functions at a new point X. */
    /*       They are part of a product that requires VLAG to be of length NDIM. */
    /*     W is a one-dimensional array that is used for working space. Its length */
    /*       must be at least 3*NDIM = 3*(NPT+N). */

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
    const long np = n + 1;
    const long nptm = npt - np;
    const long nh = n * np / 2;

    /*     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, */
    /*     BMAT and ZMAT for the first iteration, with the corresponding values of */
    /*     of NF and KOPT, which are the number of calls of CALFUN so far and the */
    /*     index of the interpolation point at the trust region centre. Then the */
    /*     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is */
    /*     less than NPT. GOPT will be updated if KOPT is different from KBASE. */

    prelim(function, n, npt, x, xl, xu, rhobeg, maxfun, xbase,
            xpt + xpt_offset, fval, gopt, hq, pq, bmat + bmat_offset,
            zmat + zmat_offset, ndim, sl, su, nf, kopt);
    xoptsq = 0.0;
    for (long i__ = 1; i__ <= n; ++i__) {
        xopt[i__] = xpt[kopt + i__ * xpt_dim1];
        xoptsq += square(xopt[i__]);
    }
    fsave = fval[1];
    if (nf < npt) {
        DEBUG_LOG("Return from BOBYQA because the objective function has been called "
                  "max_f_evals times.\n");
        goto L720;
    }
    kbase = 1;

    /*     Complete the settings that are required for the iterative procedure. */

    rho = rhobeg;
    delta = rho;
    nresc = nf;
    ntrits = 0;
    diffa = 0.0;
    diffb = 0.0;
    itest = 0;
    nfsav = nf;

    /*     Update GOPT if necessary before the first iteration and after each */
    /*     call of RESCUE that makes a call of CALFUN. */

L20:
    if (kopt != kbase) {
        ih = 0;
        for (long j = 1; j <= n; ++j) {
            for (long i__ = 1; i__ <= j; ++i__) {
                ++ih;
                if (i__ < j) {
                    gopt[j] += hq[ih] * xopt[i__];
                }
                gopt[i__] += hq[ih] * xopt[j];
            }
        }
        if (nf > npt) {
            for (long k = 1; k <= npt; ++k) {
                temp = 0.0;
                for (long j = 1; j <= n; ++j) {
                    temp += xpt[k + j * xpt_dim1] * xopt[j];
                }
                temp = pq[k] * temp;
                for (long i__ = 1; i__ <= n; ++i__) {
                    gopt[i__] += temp * xpt[k + i__ * xpt_dim1];
                }
            }
        }
    }

    /*     Generate the next point in the trust region that provides a small value */
    /*     of the quadratic model subject to the constraints on the variables. */
    /*     The long NTRITS is set to the number "trust region" iterations that */
    /*     have occurred since the last "alternative" iteration. If the length */
    /*     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to */
    /*     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW. */

L60:
    trsbox(n, npt, xpt + xpt_offset, xopt, gopt, hq, pq, sl,
            su, delta, xnew, d__, w, w + np - 1, w + np + n - 1,
            w + np + (n << 1) - 1, w + np + n * 3 - 1, &dsq, &crvmin);
    dnorm = std::min(delta, std::sqrt(dsq));
    if (dnorm < 0.5 * rho) {
        ntrits = -1;
        distsq = square(10.0 * rho);
        if (nf <= nfsav + 2) {
            goto L650;
        }

        /*     The following choice between labels 650 and 680 depends on whether or */
        /*     not our work with the current RHO seems to be complete. Either RHO is */
        /*     decreased or termination occurs if the errors in the quadratic model at */
        /*     the last three interpolation points compare favourably with predictions */
        /*     of likely improvements to the model within distance HALF*RHO of XOPT. */

        errbig = std::max(std::max(diffa, diffb), diffc);
        frhosq = rho * .125 * rho;
        if (crvmin > 0.0 && errbig > frhosq * crvmin) {
            goto L650;
        }
        bdtol = errbig / rho;
        for (long j = 1; j <= n; ++j) {
            bdtest = bdtol;
            if (xnew[j] == sl[j]) {
                bdtest = w[j];
            }
            if (xnew[j] == su[j]) {
                bdtest = -w[j];
            }
            if (bdtest < bdtol) {
                curv = hq[(j + j * j) / 2];
                for (long k = 1; k <= npt; ++k) {
                    curv += pq[k] * square(xpt[k + j * xpt_dim1]);
                }
                bdtest += 0.5 * curv * rho;
                if (bdtest < bdtol) {
                    goto L650;
                }
            }
        }
        goto L680;
    }
    ++ntrits;

    /*     Severe cancellation is likely to occur if XOPT is too far from XBASE. */
    /*     If the following test holds, then XBASE is shifted so that XOPT becomes */
    /*     zero. The appropriate changes are made to BMAT and to the second */
    /*     derivatives of the current model, beginning with the changes to BMAT */
    /*     that do not depend on ZMAT. VLAG is used temporarily for working space. */

L90:
    if (dsq <= xoptsq * .001) {
        fracsq = xoptsq * .25;
        sumpq = 0.0;
        for (long k = 1; k <= npt; ++k) {
            sumpq += pq[k];
            sum = -0.5 * xoptsq;
            for (long i__ = 1; i__ <= n; ++i__) {
                sum += xpt[k + i__ * xpt_dim1] * xopt[i__];
            }
            w[npt + k] = sum;
            temp = fracsq - 0.5 * sum;
            for (long i__ = 1; i__ <= n; ++i__) {
                w[i__] = bmat[k + i__ * bmat_dim1];
                vlag[i__] = sum * xpt[k + i__ * xpt_dim1] + temp * xopt[i__];
                const long ip = npt + i__;
                for (long j = 1; j <= i__; ++j) {
                    bmat[ip + j * bmat_dim1] = bmat[ip + j * bmat_dim1] + w[
                        i__] * vlag[j] + vlag[i__] * w[j];
                }
            }
        }

        /*     Then the revisions of BMAT that depend on ZMAT are calculated. */

        for (long jj = 1; jj <= nptm; ++jj) {
            sumz = 0.0;
            sumw = 0.0;
            for (long k = 1; k <= npt; ++k) {
                sumz += zmat[k + jj * zmat_dim1];
                vlag[k] = w[npt + k] * zmat[k + jj * zmat_dim1];
                sumw += vlag[k];
            }
            for (long j = 1; j <= n; ++j) {
                sum = (fracsq * sumz - 0.5 * sumw) * xopt[j];
                for (long k = 1; k <= npt; ++k) {
                    sum += vlag[k] * xpt[k + j * xpt_dim1];
                }
                w[j] = sum;
                for (long k = 1; k <= npt; ++k) {
                    bmat[k + j * bmat_dim1] += sum * zmat[k + jj * zmat_dim1];
                }
            }
            for (long i__ = 1; i__ <= n; ++i__) {
                const long ip = i__ + npt;
                temp = w[i__];
                for (long j = 1; j <= i__; ++j) {
                    bmat[ip + j * bmat_dim1] += temp * w[j];
                }
            }
        }

        /*     The following instructions complete the shift, including the changes */
        /*     to the second derivative parameters of the quadratic model. */

        ih = 0;
        for (long j = 1; j <= n; ++j) {
            w[j] = -0.5 * sumpq * xopt[j];
            for (long k = 1; k <= npt; ++k) {
                w[j] += pq[k] * xpt[k + j * xpt_dim1];
                xpt[k + j * xpt_dim1] -= xopt[j];
            }
            for (long i__ = 1; i__ <= j; ++i__) {
                ++ih;
                hq[ih] = hq[ih] + w[i__] * xopt[j] + xopt[i__] * w[j];
                bmat[npt + i__ + j * bmat_dim1] = bmat[npt + j + i__ *
                    bmat_dim1];
            }
        }
        for (long i__ = 1; i__ <= n; ++i__) {
            xbase[i__] += xopt[i__];
            xnew[i__] -= xopt[i__];
            sl[i__] -= xopt[i__];
            su[i__] -= xopt[i__];
            xopt[i__] = 0.0;
        }
        xoptsq = 0.0;
    }
    if (ntrits == 0) {
        goto L210;
    }
    goto L230;

    /*     XBASE is also moved to XOPT by a call of RESCUE. This calculation is */
    /*     more expensive than the previous shift, because new matrices BMAT and */
    /*     ZMAT are generated from scratch, which may include the replacement of */
    /*     interpolation points whose positions seem to be causing near linear */
    /*     dependence in the interpolation conditions. Therefore RESCUE is called */
    /*     only if rounding errors have reduced by at least a factor of two the */
    /*     denominator of the formula for updating the H matrix. It provides a */
    /*     useful safeguard, but is not invoked in most applications of BOBYQA. */

L190:
    nfsav = nf;
    kbase = kopt;
    rescue(function, n, npt, xl, xu, maxfun, xbase, xpt +  xpt_offset,
            fval, xopt, gopt, hq, pq, bmat + bmat_offset,
            zmat + zmat_offset, ndim, sl, su, nf, delta,
            kopt, vlag, w - 2, w + n + np - 1, w + ndim + np - 1);

    /*     XOPT is updated now in case the branch below to label 720 is taken. */
    /*     Any updating of GOPT occurs after the branch below to label 20, which */
    /*     leads to a trust region iteration as does the branch to label 60. */

    xoptsq = 0.0;
    if (kopt != kbase) {
        for (long i__ = 1; i__ <= n; ++i__) {
            xopt[i__] = xpt[kopt + i__ * xpt_dim1];
            xoptsq += square(xopt[i__]);
        }
    }
    if (nf < 0) {
        nf = maxfun;
        DEBUG_LOG("Return from BOBYQA because the objective function has been called "
                  "max_f_evals times.\n");
        goto L720;
    }
    nresc = nf;
    if (nfsav < nf) {
        nfsav = nf;
        goto L20;
    }
    if (ntrits > 0) {
        goto L60;
    }

    /*     Pick two alternative vectors of variables, relative to XBASE, that */
    /*     are suitable as new positions of the KNEW-th interpolation point. */
    /*     Firstly, XNEW is set to the point on a line through XOPT and another */
    /*     interpolation point that minimizes the predicted value of the next */
    /*     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL */
    /*     and SU bounds. Secondly, XALT is set to the best feasible point on */
    /*     a constrained version of the Cauchy step of the KNEW-th Lagrange */
    /*     function, the corresponding value of the square of this function */
    /*     being returned in CAUCHY. The choice between these alternatives is */
    /*     going to be made when the denominator is calculated. */

L210:
    altmov(n, npt, xpt + xpt_offset, xopt, bmat + bmat_offset, zmat + zmat_offset,
            ndim, sl, su, kopt, knew, adelt, xnew,
            xalt, alpha, cauchy, w, w + np - 1, w + ndim);
    for (long i__ = 1; i__ <= n; ++i__) {
        d__[i__] = xnew[i__] - xopt[i__];
    }

    /*     Calculate VLAG and BETA for the current choice of D. The scalar */
    /*     product of D with XPT(K,.) is going to be held in W(NPT+K) for */
    /*     use when VQUAD is calculated. */

L230:
    for (long k = 1; k <= npt; ++k) {
        suma = 0.0;
        sumb = 0.0;
        sum = 0.0;
        for (long j = 1; j <= n; ++j) {
            suma += xpt[k + j * xpt_dim1] * d__[j];
            sumb += xpt[k + j * xpt_dim1] * xopt[j];
            sum += bmat[k + j * bmat_dim1] * d__[j];
        }
        w[k] = suma * (0.5 * suma + sumb);
        vlag[k] = sum;
        w[npt + k] = suma;
    }
    beta = 0.0;
    for (long jj = 1; jj <= nptm; ++jj) {
        sum = 0.0;
        for (long k = 1; k <= npt; ++k) {
            sum += zmat[k + jj * zmat_dim1] * w[k];
        }
        beta -= sum * sum;
        for (long k = 1; k <= npt; ++k) {
            vlag[k] += sum * zmat[k + jj * zmat_dim1];
        }
    }
    dsq = 0.0;
    bsum = 0.0;
    dx = 0.0;
    for (long j = 1; j <= n; ++j) {
        dsq += square(d__[j]);
        sum = 0.0;
        for (long k = 1; k <= npt; ++k) {
            sum += w[k] * bmat[k + j * bmat_dim1];
        }
        bsum += sum * d__[j];
        jp = npt + j;
        for (long i__ = 1; i__ <= n; ++i__) {
            sum += bmat[jp + i__ * bmat_dim1] * d__[i__];
        }
        vlag[jp] = sum;
        bsum += sum * d__[j];
        dx += d__[j] * xopt[j];
    }
    beta = dx * dx + dsq * (xoptsq + dx + dx + 0.5 * dsq) + beta - bsum;
    vlag[kopt] += 1.0;

    /*     If NTRITS is zero, the denominator may be increased by replacing */
    /*     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if */
    /*     rounding errors have damaged the chosen denominator. */

    if (ntrits == 0) {
        denom = square(vlag[knew]) + alpha * beta;
        if (denom < cauchy && cauchy > 0.0) {
            for (long i__ = 1; i__ <= n; ++i__) {
                xnew[i__] = xalt[i__];
                d__[i__] = xnew[i__] - xopt[i__];
            }
            cauchy = 0.0;
            goto L230;
        }
        if (denom <= 0.5 * square(vlag[knew])) {
            if (nf > nresc) {
                goto L190;
            }
            DEBUG_LOG("Return from BOBYQA because of much cancellation in a denominator.\n");
            goto L720;
        }

        /*     Alternatively, if NTRITS is positive, then set KNEW to the index of */
        /*     the next interpolation point to be deleted to make room for a trust */
        /*     region step. Again RESCUE may be called if rounding errors have damaged */
        /*     the chosen denominator, which is the reason for attempting to select */
        /*     KNEW before calculating the next value of the objective function. */

    } else {
        delsq = delta * delta;
        scaden = 0.0;
        biglsq = 0.0;
        knew = 0;
        for (long k = 1; k <= npt; ++k) {
            if (k == kopt) {
                goto L350;
            }
            hdiag = 0.0;
            for (long jj = 1; jj <= nptm; ++jj) {
                hdiag += square(zmat[k + jj * zmat_dim1]);
            }
            den = beta * hdiag + square(vlag[k]);
            distsq = 0.0;
            for (long j = 1; j <= n; ++j) {
                distsq += square(xpt[k + j * xpt_dim1] - xopt[j]);
            }
            temp = std::max(1.0, square(distsq / delsq));
            if (temp * den > scaden) {
                scaden = temp * den;
                knew = k;
                denom = den;
            }
            biglsq = std::max(biglsq, temp * square(vlag[k]));
L350:
            ;
        }
        if (scaden <= 0.5 * biglsq) {
            if (nf > nresc) {
                goto L190;
            }
            DEBUG_LOG("Return from BOBYQA because of much cancellation in a denominator.\n");
            goto L720;
        }
    }

    /*     Put the variables for the next calculation of the objective function */
    /*       in XNEW, with any adjustments for the bounds. */


    /*     Calculate the value of the objective function at XBASE+XNEW, unless */
    /*       the limit on the number of calculations of F has been reached. */

L360:
    for (long i__ = 1; i__ <= n; ++i__) {
        x[i__] = std::min(std::max(xl[i__], xbase[i__] + xnew[i__]), xu[i__]);
        if (xnew[i__] == sl[i__]) {
            x[i__] = xl[i__];
        }
        if (xnew[i__] == su[i__]) {
            x[i__] = xu[i__];
        }
    }
    if (nf >= maxfun) {
        DEBUG_LOG("Return from BOBYQA because the objective function has been called "
                  "max_f_evals times.\n");
        goto L720;
    }
    ++nf;
    f = function(n, x + 1);
    if (ntrits == -1) {
        fsave = f;
        goto L720;
    }

    /*     Use the quadratic model to predict the change in F due to the step D, */
    /*       and set DIFF to the error of this prediction. */

    fopt = fval[kopt];
    vquad = 0.0;
    ih = 0;
    for (long j = 1; j <= n; ++j) {
        vquad += d__[j] * gopt[j];
        for (long i__ = 1; i__ <= j; ++i__) {
            ++ih;
            temp = d__[i__] * d__[j];
            if (i__ == j) {
                temp = 0.5 * temp;
            }
            vquad += hq[ih] * temp;
        }
    }
    for (long k = 1; k <= npt; ++k) {
        vquad += 0.5 * pq[k] * square(w[npt + k]);
    }
    diff = f - fopt - vquad;
    diffc = diffb;
    diffb = diffa;
    diffa = std::abs(diff);
    if (dnorm > rho) {
        nfsav = nf;
    }

    /*     Pick the next value of DELTA after a trust region step. */

    if (ntrits > 0) {
        if (vquad >= 0.0) {
            DEBUG_LOG("Return from BOBYQA because a trust region step has failed to reduce Q.\n");
            goto L720;
        }
        ratio = (f - fopt) / vquad;
        if (ratio <= 0.1) {
            delta = std::min(0.5 * delta, dnorm);
        } else if (ratio <= .7) {
            delta = std::max(0.5 * delta, dnorm);
        } else {
            delta = std::max(0.5 * delta, dnorm + dnorm);
        }
        if (delta <= rho * 1.5) {
            delta = rho;
        }

        /*     Recalculate KNEW and DENOM if the new F is less than FOPT. */

        if (f < fopt) {
            ksav = knew;
            densav = denom;
            delsq = delta * delta;
            scaden = 0.0;
            biglsq = 0.0;
            knew = 0;
            for (long k = 1; k <= npt; ++k) {
                hdiag = 0.0;
                for (long jj = 1; jj <= nptm; ++jj) {
                    hdiag += square(zmat[k + jj * zmat_dim1]);
                }
                den = beta * hdiag + square(vlag[k]);
                distsq = 0.0;
                for (long j = 1; j <= n; ++j) {
                    distsq += square(xpt[k + j * xpt_dim1] - xnew[j]);
                }
                temp = std::max(1.0, square(distsq / delsq));
                if (temp * den > scaden) {
                    scaden = temp * den;
                    knew = k;
                    denom = den;
                }
                biglsq = std::max(biglsq, temp * square(vlag[k]));
            }
            if (scaden <= 0.5 * biglsq) {
                knew = ksav;
                denom = densav;
            }
        }
    }

    /*     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be */
    /*     moved. Also update the second derivative terms of the model. */

    update(n, npt, bmat + bmat_offset, zmat + zmat_offset, ndim, vlag,
            beta, denom, knew, w);
    ih = 0;
    pqold = pq[knew];
    pq[knew] = 0.0;
    for (long i__ = 1; i__ <= n; ++i__) {
        temp = pqold * xpt[knew + i__ * xpt_dim1];
        for (long j = 1; j <= i__; ++j) {
            ++ih;
            hq[ih] += temp * xpt[knew + j * xpt_dim1];
        }
    }
    for (long jj = 1; jj <= nptm; ++jj) {
        temp = diff * zmat[knew + jj * zmat_dim1];
        for (long k = 1; k <= npt; ++k) {
            pq[k] += temp * zmat[k + jj * zmat_dim1];
        }
    }

    /*     Include the new interpolation point, and make the changes to GOPT at */
    /*     the old XOPT that are caused by the updating of the quadratic model. */

    fval[knew] = f;
    for (long i__ = 1; i__ <= n; ++i__) {
        xpt[knew + i__ * xpt_dim1] = xnew[i__];
        w[i__] = bmat[knew + i__ * bmat_dim1];
    }
    for (long k = 1; k <= npt; ++k) {
        suma = 0.0;
        for (long jj = 1; jj <= nptm; ++jj) {
            suma += zmat[knew + jj * zmat_dim1] * zmat[k + jj * zmat_dim1];
        }
        sumb = 0.0;
        for (long j = 1; j <= n; ++j) {
            sumb += xpt[k + j * xpt_dim1] * xopt[j];
        }
        temp = suma * sumb;
        for (long i__ = 1; i__ <= n; ++i__) {
            w[i__] += temp * xpt[k + i__ * xpt_dim1];
        }
    }
    for (long i__ = 1; i__ <= n; ++i__) {
        gopt[i__] += diff * w[i__];
    }

    /*     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT. */

    if (f < fopt) {
        kopt = knew;
        xoptsq = 0.0;
        ih = 0;
        for (long j = 1; j <= n; ++j) {
            xopt[j] = xnew[j];
            xoptsq += square(xopt[j]);
            for (long i__ = 1; i__ <= j; ++i__) {
                ++ih;
                if (i__ < j) {
                    gopt[j] += hq[ih] * d__[i__];
                }
                gopt[i__] += hq[ih] * d__[j];
            }
        }
        for (long k = 1; k <= npt; ++k) {
            temp = 0.0;
            for (long j = 1; j <= n; ++j) {
                temp += xpt[k + j * xpt_dim1] * d__[j];
            }
            temp = pq[k] * temp;
            for (long i__ = 1; i__ <= n; ++i__) {
                gopt[i__] += temp * xpt[k + i__ * xpt_dim1];
            }
        }
    }

    /*     Calculate the parameters of the least Frobenius norm interpolant to */
    /*     the current data, the gradient of this interpolant at XOPT being put */
    /*     into VLAG(NPT+I), I=1,2,...,N. */

    if (ntrits > 0) {
        for (long k = 1; k <= npt; ++k) {
            vlag[k] = fval[k] - fval[kopt];
            w[k] = 0.0;
        }
        for (long j = 1; j <= nptm; ++j) {
            sum = 0.0;
            for (long k = 1; k <= npt; ++k) {
                sum += zmat[k + j * zmat_dim1] * vlag[k];
            }
            for (long k = 1; k <= npt; ++k) {
                w[k] += sum * zmat[k + j * zmat_dim1];
            }
        }
        for (long k = 1; k <= npt; ++k) {
            sum = 0.0;
            for (long j = 1; j <= n; ++j) {
                sum += xpt[k + j * xpt_dim1] * xopt[j];
            }
            w[k + npt] = w[k];
            w[k] = sum * w[k];
        }
        gqsq = 0.0;
        gisq = 0.0;
        for (long i__ = 1; i__ <= n; ++i__) {
            sum = 0.0;
            for (long k = 1; k <= npt; ++k) {
                sum = sum + bmat[k + i__ * bmat_dim1] * vlag[k] + xpt[k + i__
                    * xpt_dim1] * w[k];
            }
            if (xopt[i__] == sl[i__]) {
                gqsq += square(std::min(0.0, gopt[i__]));
                gisq += square(std::min(0.0, sum));
            } else if (xopt[i__] == su[i__]) {
                gqsq += square(std::max(0.0, gopt[i__]));
                gisq += square(std::max(0.0, sum));
            } else {
                gqsq += square(gopt[i__]);
                gisq += sum * sum;
            }
            vlag[npt + i__] = sum;
        }

        /*     Test whether to replace the new quadratic model by the least Frobenius */
        /*     norm interpolant, making the replacement if the test is satisfied. */

        ++itest;
        if (gqsq < 10.0 * gisq) {
            itest = 0;
        }
        if (itest >= 3) {
            const long i_n = std::max(npt, nh);
            for (long i__ = 1; i__ <= i_n; ++i__) {
                if (i__ <= n) {
                    gopt[i__] = vlag[npt + i__];
                }
                if (i__ <= npt) {
                    pq[i__] = w[npt + i__];
                }
                if (i__ <= nh) {
                    hq[i__] = 0.0;
                }
                itest = 0;
            }
        }
    }

    /*     If a trust region step has provided a sufficient decrease in F, then */
    /*     branch for another trust region calculation. The case NTRITS=0 occurs */
    /*     when the new interpolation point was reached by an alternative step. */

    if (ntrits == 0) {
        goto L60;
    }
    if (f <= fopt + 0.1 * vquad) {
        goto L60;
    }

    /*     Alternatively, find out if the interpolation points are close enough */
    /*       to the best point so far. */

    distsq = std::max(square(2.0 * delta), square(10.0 * rho));
L650:
    knew = 0;
    for (long k = 1; k <= npt; ++k) {
        sum = 0.0;
        for (long j = 1; j <= n; ++j) {
            sum += square(xpt[k + j * xpt_dim1] - xopt[j]);
        }
        if (sum > distsq) {
            knew = k;
            distsq = sum;
        }
    }

    /*     If KNEW is positive, then ALTMOV finds alternative new positions for */
    /*     the KNEW-th interpolation point within distance ADELT of XOPT. It is */
    /*     reached via label 90. Otherwise, there is a branch to label 60 for */
    /*     another trust region iteration, unless the calculations with the */
    /*     current RHO are complete. */

    if (knew > 0) {
        dist = std::sqrt(distsq);
        if (ntrits == -1) {
            delta = std::min(0.1 * delta, 0.5 * dist);
            if (delta <= rho * 1.5) {
                delta = rho;
            }
        }
        ntrits = 0;
        adelt = std::max(std::min(0.1 * dist, delta), rho);
        dsq = adelt * adelt;
        goto L90;
    }
    if (ntrits == -1) {
        goto L680;
    }
    if (ratio > 0.0) {
        goto L60;
    }
    if (std::max(delta,dnorm) > rho) {
        goto L60;
    }

    /*     The calculations with the current value of RHO are complete. Pick the */
    /*       next values of RHO and DELTA. */

L680:
    if (rho > rhoend) {
        delta = 0.5 * rho;
        ratio = rho / rhoend;
        if (ratio <= 16.) {
            rho = rhoend;
        } else if (ratio <= 250.) {
            rho = std::sqrt(ratio) * rhoend;
        } else {
            rho = 0.1 * rho;
        }
        delta = std::max(delta,rho);
        ntrits = 0;
        nfsav = nf;
        goto L60;
    }

    /*     Return from the calculation, after another Newton-Raphson step, if */
    /*       it is too short to have been tried before. */

    if (ntrits == -1) {
        goto L360;
    }
L720:
    if (fval[kopt] <= fsave) {
        for (long i__ = 1; i__ <= n; ++i__) {
            x[i__] = std::min(std::max(xl[i__], xbase[i__] + xopt[i__]), xu[i__]);
            if (xopt[i__] == sl[i__]) {
                x[i__] = xl[i__];
            }
            if (xopt[i__] == su[i__]) {
                x[i__] = xu[i__];
            }
        }
        f = fval[kopt];
    }

    return f;
}

} // namespace bobyqa_detail
