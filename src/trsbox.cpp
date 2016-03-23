#include "utils.hpp"

#include <algorithm>

namespace bobyqa_detail {

void trsbox(
    const long n,
    const long npt,
    const double *xpt,
    const double *const xopt,
    const double *const gopt,
    const double *const hq,
    const double *const pq,
    const double *const sl,
    const double *const su,
    const double delta,
    double *const xnew,
    double *const d,
    double *const gnew,
    double *const xbdi,
    double *const s,
    double *const hs,
    double *const hred,
    double *const dsq,
    double *const crvmin
)
{
    /* Local variables */
    double ds;
    long iu;
    double dhd, dhs, cth, shs, sth, ssq, beta, sdec, blen;
    long iact = 0, nact = 0;
    double angt, qred;
    long isav;
    double temp = 0, xsav = 0, xsum = 0, angbd = 0, dredg = 0, sredg = 0;
    long iterc;
    double resid = 0, delsq = 0, ggsav = 0, tempa = 0, tempb = 0,
               redmax = 0, dredsq = 0, redsav = 0, gredsq = 0, rednew = 0;
    long itcsav = 0;
    double rdprev = 0, rdnext = 0, stplen = 0, stepsq = 0;
    long itermax = 0;


    /*     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same */
    /*       meanings as the corresponding arguments of BOBYQB. */
    /*     DELTA is the trust region radius for the present calculation, which */
    /*       seeks a small value of the quadratic model within distance DELTA of */
    /*       XOPT subject to the bounds on the variables. */
    /*     XNEW will be set to a new vector of variables that is approximately */
    /*       the one that minimizes the quadratic model within the trust region */
    /*       subject to the SL and SU constraints on the variables. It satisfies */
    /*       as equations the bounds that become active during the calculation. */
    /*     D is the calculated trial step from XOPT, generated iteratively from an */
    /*       initial value of zero. Thus XNEW is XOPT+D after the final iteration. */
    /*     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated */
    /*       when D is updated. */
    /*     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is */
    /*       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the */
    /*       I-th variable has become fixed at a bound, the bound being SL(I) or */
    /*       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This */
    /*       information is accumulated during the construction of XNEW. */
    /*     The arrays S, HS and HRED are also used for working space. They hold the */
    /*       current search direction, and the changes in the gradient of Q along S */
    /*       and the reduced D, respectively, where the reduced D is the same as D, */
    /*       except that the components of the fixed variables are zero. */
    /*     DSQ will be set to the square of the length of XNEW-XOPT. */
    /*     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise */
    /*       it is set to the least curvature of H that occurs in the conjugate */
    /*       gradient searches that are not restricted by any constraints. The */
    /*       value CRVMIN=-1.0D0 is set, however, if all of these searches are */
    /*       constrained. */

    /*     A version of the truncated conjugate gradient is applied. If a line */
    /*     search is restricted by a constraint, then the procedure is restarted, */
    /*     the values of the variables that are at their bounds being fixed. If */
    /*     the trust region boundary is reached, then further changes may be made */
    /*     to D, each one being in the two dimensional space that is spanned */
    /*     by the current D and the gradient of Q at XOPT+D, staying on the trust */
    /*     region boundary. Termination occurs when the reduction in Q seems to */
    /*     be close to the greatest reduction that can be achieved. */

    /*     Set some constants. */

    /* Parameter adjustments */
    const long xpt_dim1 = npt;
    const long xpt_offset = 1 + xpt_dim1;
    xpt -= xpt_offset;

    /* Function Body */

    /*     The sign of GOPT(I) gives the sign of the change to the I-th variable */
    /*     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether */
    /*     or not to fix the I-th variable at one of its bounds initially, with */
    /*     NACT being set to the number of fixed variables. D and GNEW are also */
    /*     set for the first iteration. DELSQ is the upper bound on the sum of */
    /*     squares of the free variables. QRED is the reduction in Q so far. */

    iterc = 0;
    nact = 0;
    for (long i = 1; i <= n; ++i) {
        xbdi[i] = 0.0;
        if (xopt[i] <= sl[i]) {
            if (gopt[i] >= 0.0) {
                xbdi[i] = -1.0;
            }
        } else if (xopt[i] >= su[i]) {
            if (gopt[i] <= 0.0) {
                xbdi[i] = 1.0;
            }
        }
        if (xbdi[i] != 0.0) {
            ++nact;
        }
        d[i] = 0.0;
        gnew[i] = gopt[i];
    }
    delsq = delta * delta;
    qred = 0.0;
    *crvmin = -1.0;

    /*     Set the next search direction of the conjugate gradient method. It is */
    /*     the steepest descent direction initially and when the iterations are */
    /*     restarted because a variable has just been fixed by a bound, and of */
    /*     course the components of the fixed variables are zero. ITERMAX is an */
    /*     upper bound on the indices of the conjugate gradient iterations. */

L20:
    beta = 0.0;
L30:
    stepsq = 0.0;
    for (long i = 1; i <= n; ++i) {
        if (xbdi[i] != 0.0) {
            s[i] = 0.0;
        } else if (beta == 0.0) {
            s[i] = -gnew[i];
        } else {
            s[i] = beta * s[i] - gnew[i];
        }
        stepsq += square(s[i]);
    }
    if (stepsq == 0.0) {
        goto L190;
    }
    if (beta == 0.0) {
        gredsq = stepsq;
        itermax = iterc + n - nact;
    }
    if (gredsq * delsq <= qred * 1e-4 * qred) {
        goto L190;
    }

    /*     Multiply the search direction by the second derivative matrix of Q and */
    /*     calculate some scalars for the choice of steplength. Then set BLEN to */
    /*     the length of the the step to the trust region boundary and STPLEN to */
    /*     the steplength, ignoring the simple bounds. */

    goto L210;
L50:
    resid = delsq;
    ds = 0.0;
    shs = 0.0;
    for (long i = 1; i <= n; ++i) {
        if (xbdi[i] == 0.0) {
            resid -= square(d[i]);
            ds += s[i] * d[i];
            shs += s[i] * hs[i];
        }
    }
    if (resid <= 0.0) {
        goto L90;
    }
    temp = std::sqrt(stepsq * resid + ds * ds);
    if (ds < 0.0) {
        blen = (temp - ds) / stepsq;
    } else {
        blen = resid / (temp + ds);
    }
    stplen = blen;
    if (shs > 0.0) {
        stplen = std::min(blen, gredsq / shs);
    }

    /*     Reduce STPLEN if necessary in order to preserve the simple bounds, */
    /*     letting IACT be the index of the new constrained variable. */

    iact = 0;
    for (long i = 1; i <= n; ++i) {
        if (s[i] != 0.0) {
            xsum = xopt[i] + d[i];
            if (s[i] > 0.0) {
                temp = (su[i] - xsum) / s[i];
            } else {
                temp = (sl[i] - xsum) / s[i];
            }
            if (temp < stplen) {
                stplen = temp;
                iact = i;
            }
        }
    }

    /*     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q. */

    sdec = 0.0;
    if (stplen > 0.0) {
        ++iterc;
        temp = shs / stepsq;
        if (iact == 0 && temp > 0.0) {
            *crvmin = std::min(*crvmin,temp);
            if (*crvmin == -1.0) {
                *crvmin = temp;
            }
        }
        ggsav = gredsq;
        gredsq = 0.0;
        for (long i = 1; i <= n; ++i) {
            gnew[i] += stplen * hs[i];
            if (xbdi[i] == 0.0) {
                gredsq += square(gnew[i]);
            }
            d[i] += stplen * s[i];
        }
        sdec = std::max(stplen * (ggsav - 0.5 * stplen * shs), 0.0);
        qred += sdec;
    }

    /*     Restart the conjugate gradient method if it has hit a new bound. */

    if (iact > 0) {
        ++nact;
        xbdi[iact] = 1.0;
        if (s[iact] < 0.0) {
            xbdi[iact] = -1.0;
        }
        delsq -= square(d[iact]);
        if (delsq <= 0.0) {
            goto L90;
        }
        goto L20;
    }

    /*     If STPLEN is less than BLEN, then either apply another conjugate */
    /*     gradient iteration or RETURN. */

    if (stplen < blen) {
        if (iterc == itermax) {
            goto L190;
        }
        if (sdec <= qred * .01) {
            goto L190;
        }
        beta = gredsq / ggsav;
        goto L30;
    }
L90:
    *crvmin = 0.0;

    /*     Prepare for the alternative iteration by calculating some scalars */
    /*     and by multiplying the reduced D by the second derivative matrix of */
    /*     Q, where S holds the reduced D in the call of GGMULT. */

L100:
    if (nact >= n - 1) {
        goto L190;
    }
    dredsq = 0.0;
    dredg = 0.0;
    gredsq = 0.0;
    for (long i = 1; i <= n; ++i) {
        if (xbdi[i] == 0.0) {
            dredsq += square(d[i]);
            dredg += d[i] * gnew[i];
            gredsq += square(gnew[i]);
            s[i] = d[i];
        } else {
            s[i] = 0.0;
        }
    }
    itcsav = iterc;
    goto L210;

    /*     Let the search direction S be a linear combination of the reduced D */
    /*     and the reduced G that is orthogonal to the reduced D. */

L120:
    ++iterc;
    temp = gredsq * dredsq - dredg * dredg;
    if (temp <= qred * 1e-4 * qred) {
        goto L190;
    }
    temp = std::sqrt(temp);
    for (long i = 1; i <= n; ++i) {
        if (xbdi[i] == 0.0) {
            s[i] = (dredg * d[i] - dredsq * gnew[i]) / temp;
        } else {
            s[i] = 0.0;
        }
    }
    sredg = -temp;

    /*     By considering the simple bounds on the variables, calculate an upper */
    /*     bound on the tangent of half the angle of the alternative iteration, */
    /*     namely ANGBD, except that, if already a free variable has reached a */
    /*     bound, there is a branch back to label 100 after fixing that variable. */

    angbd = 1.0;
    iact = 0;
    for (long i = 1; i <= n; ++i) {
        if (xbdi[i] == 0.0) {
            tempa = xopt[i] + d[i] - sl[i];
            tempb = su[i] - xopt[i] - d[i];
            if (tempa <= 0.0) {
                ++nact;
                xbdi[i] = -1.0;
                goto L100;
            } else if (tempb <= 0.0) {
                ++nact;
                xbdi[i] = 1.0;
                goto L100;
            }
            ssq = square(d[i]) + square(s[i]);
            temp = ssq - square(xopt[i] - sl[i]);
            if (temp > 0.0) {
                temp = std::sqrt(temp) - s[i];
                if (angbd * temp > tempa) {
                    angbd = tempa / temp;
                    iact = i;
                    xsav = -1.0;
                }
            }
            temp = ssq - square(su[i] - xopt[i]);
            if (temp > 0.0) {
                temp = std::sqrt(temp) + s[i];
                if (angbd * temp > tempb) {
                    angbd = tempb / temp;
                    iact = i;
                    xsav = 1.0;
                }
            }
        }
    }

    /*     Calculate HHD and some curvatures for the alternative iteration. */

    goto L210;
L150:
    shs = 0.0;
    dhs = 0.0;
    dhd = 0.0;
    for (long i = 1; i <= n; ++i) {
        if (xbdi[i] == 0.0) {
            shs += s[i] * hs[i];
            dhs += d[i] * hs[i];
            dhd += d[i] * hred[i];
        }
    }

    /*     Seek the greatest reduction in Q for a range of equally spaced values */
    /*     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of */
    /*     the alternative iteration. */

    redmax = 0.0;
    isav = 0;
    redsav = 0.0;
    iu = long(angbd * 17. + 3.1);
    for (long i = 1; i <= iu; ++i) {
        angt = angbd * double(i) / double(iu);
        sth = (angt + angt) / (1.0 + angt * angt);
        temp = shs + angt * (angt * dhd - dhs - dhs);
        rednew = sth * (angt * dredg - sredg - 0.5 * sth * temp);
        if (rednew > redmax) {
            redmax = rednew;
            isav = i;
            rdprev = redsav;
        } else if (i == isav + 1) {
            rdnext = rednew;
        }
        redsav = rednew;
    }

    /*     Return if the reduction is zero. Otherwise, set the sine and cosine */
    /*     of the angle of the alternative iteration, and calculate SDEC. */

    if (isav == 0) {
        goto L190;
    }
    if (isav < iu) {
        temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext);
        angt = angbd * (double(isav) + 0.5 * temp) / double(iu);
    }
    cth = (1.0 - angt * angt) / (1.0 + angt * angt);
    sth = (angt + angt) / (1.0 + angt * angt);
    temp = shs + angt * (angt * dhd - dhs - dhs);
    sdec = sth * (angt * dredg - sredg - 0.5 * sth * temp);
    if (sdec <= 0.0) {
        goto L190;
    }

    /*     Update GNEW, D and HRED. If the angle of the alternative iteration */
    /*     is restricted by a bound on a free variable, that variable is fixed */
    /*     at the bound. */

    dredg = 0.0;
    gredsq = 0.0;
    for (long i = 1; i <= n; ++i) {
        gnew[i] = gnew[i] + (cth - 1.0) * hred[i] + sth * hs[i];
        if (xbdi[i] == 0.0) {
            d[i] = cth * d[i] + sth * s[i];
            dredg += d[i] * gnew[i];
            gredsq += square(gnew[i]);
        }
        hred[i] = cth * hred[i] + sth * hs[i];
    }
    qred += sdec;
    if (iact > 0 && isav == iu) {
        ++nact;
        xbdi[iact] = xsav;
        goto L100;
    }

    /*     If SDEC is sufficiently small, then RETURN after setting XNEW to */
    /*     XOPT+D, giving careful attention to the bounds. */

    if (sdec > qred * .01) {
        goto L120;
    }
L190:
    *dsq = 0.0;
    for (long i = 1; i <= n; ++i) {
        xnew[i] = std::max(std::min(xopt[i] + d[i], su[i]), sl[i]);
        if (xbdi[i] == -1.0) {
            xnew[i] = sl[i];
        }
        if (xbdi[i] == 1.0) {
            xnew[i] = su[i];
        }
        d[i] = xnew[i] - xopt[i];
        *dsq += square(d[i]);
    }
    return;
    /*     The following instructions multiply the current S-vector by the second */
    /*     derivative matrix of the quadratic model, putting the product in HS. */
    /*     They are reached from three different parts of the software above and */
    /*     they can be regarded as an external subroutine. */

L210:
    long ih = 0;
    for (long j = 1; j <= n; ++j) {
        hs[j] = 0.0;
        for (long i = 1; i <= j; ++i) {
            ++ih;
            if (i < j) {
                hs[j] += hq[ih] * s[i];
            }
            hs[i] += hq[ih] * s[j];
        }
    }
    for (long k = 1; k <= npt; ++k) {
        if (pq[k] != 0.0) {
            temp = 0.0;
            for (long j = 1; j <= n; ++j) {
                temp += xpt[k + j * xpt_dim1] * s[j];
            }
            temp *= pq[k];
            for (long i = 1; i <= n; ++i) {
                hs[i] += temp * xpt[k + i * xpt_dim1];
            }
        }
    }
    if (*crvmin != 0.0) {
        goto L50;
    }
    if (iterc > itcsav) {
        goto L150;
    }
    for (long i = 1; i <= n; ++i) {
        hred[i] = hs[i];
    }
    goto L120;
}

} // namespace bobyqa_detail
