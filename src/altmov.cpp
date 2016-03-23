#include "utils.hpp"

#include <algorithm>

namespace bobyqa_detail {

const double one_plus_sqrt_2 = 1.0 + sqrt_2;

void altmov(
    const long n,
    const long npt,
    const double *xpt,
    const double *const xopt,
    const double *bmat,
    const double *zmat,
    const long ndim,
    const double *const sl,
    const double *const su,
    const long kopt,
    const long knew,
    const double adelt,
    double *const xnew,
    double *const xalt,
    double& alpha,
    double& cauchy,
    double *const glag,
    double *const hcol,
    double *const w
)
{
    /* Local variables */
    double gw, diff;
    long ilbd, isbd;
    double slbd;
    long iubd;
    double vlag, subd, temp;
    long ksav = 0;
    double step = 0, curv = 0;
    long iflag;
    double scale = 0, csave = 0, tempa = 0, tempb = 0, tempd = 0, sumin = 0, ggfree = 0;
    long ibdsav = 0;
    double dderiv = 0, bigstp = 0, predsq = 0, presav = 0, distsq = 0, stpsav = 0, wfixsq = 0, wsqsav = 0;

    /*     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have */
    /*       the same meanings as the corresponding arguments of BOBYQB. */
    /*     KOPT is the index of the optimal interpolation point. */
    /*     KNEW is the index of the interpolation point that is going to be moved. */
    /*     ADELT is the current trust region bound. */
    /*     XNEW will be set to a suitable new position for the interpolation point */
    /*       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region */
    /*       bounds and it should provide a large denominator in the next call of */
    /*       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the */
    /*       straight lines through XOPT and another interpolation point. */
    /*     XALT also provides a large value of the modulus of the KNEW-th Lagrange */
    /*       function subject to the constraints that have been mentioned, its main */
    /*       difference from XNEW being that XALT-XOPT is a constrained version of */
    /*       the Cauchy step within the trust region. An exception is that XALT is */
    /*       not calculated if all components of GLAG (see below) are zero. */
    /*     ALPHA will be set to the KNEW-th diagonal element of the H matrix. */
    /*     CAUCHY will be set to the square of the KNEW-th Lagrange function at */
    /*       the step XALT-XOPT from XOPT for the vector XALT that is returned, */
    /*       except that CAUCHY is set to zero if XALT is not calculated. */
    /*     GLAG is a working space vector of length N for the gradient of the */
    /*       KNEW-th Lagrange function at XOPT. */
    /*     HCOL is a working space vector of length NPT for the second derivative */
    /*       coefficients of the KNEW-th Lagrange function. */
    /*     W is a working space vector of length 2N that is going to hold the */
    /*       constrained Cauchy step from XOPT of the Lagrange function, followed */
    /*       by the downhill version of XALT when the uphill step is calculated. */

    /*     Set the first NPT components of W to the leading elements of the */
    /*     KNEW-th column of the H matrix. */

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
    for (long k = 1; k <= npt; ++k) {
        hcol[k] = 0.0;
    }
    const long j_n = npt - n - 1;
    for (long j = 1; j <= j_n; ++j) {
        temp = zmat[knew + j * zmat_dim1];
        for (long k = 1; k <= npt; ++k) {
            hcol[k] += temp * zmat[k + j * zmat_dim1];
        }
    }
    alpha = hcol[knew];
    const double ha = 0.5 * alpha;

    /*     Calculate the gradient of the KNEW-th Lagrange function at XOPT. */

    for (long i = 1; i <= n; ++i) {
        glag[i] = bmat[knew + i * bmat_dim1];
    }
    for (long k = 1; k <= npt; ++k) {
        temp = 0.0;
        for (long j = 1; j <= n; ++j) {
            temp += xpt[k + j * xpt_dim1] * xopt[j];
        }
        temp = hcol[k] * temp;
        for (long i = 1; i <= n; ++i) {
            glag[i] += temp * xpt[k + i * xpt_dim1];
        }
    }

    /*     Search for a large denominator along the straight lines through XOPT */
    /*     and another interpolation point. SLBD and SUBD will be lower and upper */
    /*     bounds on the step along each of these lines in turn. PREDSQ will be */
    /*     set to the square of the predicted denominator for each line. PRESAV */
    /*     will be set to the largest admissible value of PREDSQ that occurs. */

    presav = 0.0;
    for (long k = 1; k <= npt; ++k) {
        if (k == kopt) {
            goto L80;
        }
        dderiv = 0.0;
        distsq = 0.0;
        for (long i = 1; i <= n; ++i) {
            temp = xpt[k + i * xpt_dim1] - xopt[i];
            dderiv += glag[i] * temp;
            distsq += temp * temp;
        }
        subd = adelt / std::sqrt(distsq);
        slbd = -subd;
        ilbd = 0;
        iubd = 0;
        sumin = std::min(1.0,subd);

        /*     Revise SLBD and SUBD if necessary because of the bounds in SL and SU. */

        for (long i = 1; i <= n; ++i) {
            temp = xpt[k + i * xpt_dim1] - xopt[i];
            if (temp > 0.0) {
                if (slbd * temp < sl[i] - xopt[i]) {
                    slbd = (sl[i] - xopt[i]) / temp;
                    ilbd = -i;
                }
                if (subd * temp > su[i] - xopt[i]) {
                    subd = std::max(sumin, (su[i] - xopt[i]) / temp);
                    iubd = i;
                }
            } else if (temp < 0.0) {
                if (slbd * temp > su[i] - xopt[i]) {
                    slbd = (su[i] - xopt[i]) / temp;
                    ilbd = i;
                }
                if (subd * temp < sl[i] - xopt[i]) {
                    subd = std::max(sumin, (sl[i] - xopt[i]) / temp);
                    iubd = -i;
                }
            }
        }

        /*     Seek a large modulus of the KNEW-th Lagrange function when the index */
        /*     of the other interpolation point on the line through XOPT is KNEW. */

        if (k == knew) {
            diff = dderiv - 1.0;
            step = slbd;
            vlag = slbd * (dderiv - slbd * diff);
            isbd = ilbd;
            temp = subd * (dderiv - subd * diff);
            if (std::abs(temp) > std::abs(vlag)) {
                step = subd;
                vlag = temp;
                isbd = iubd;
            }
            tempd = 0.5 * dderiv;
            tempa = tempd - diff * slbd;
            tempb = tempd - diff * subd;
            if (tempa * tempb < 0.0) {
                temp = tempd * tempd / diff;
                if (std::abs(temp) > std::abs(vlag)) {
                    step = tempd / diff;
                    vlag = temp;
                    isbd = 0;
                }
            }

            /*     Search along each of the other lines through XOPT and another point. */

        } else {
            step = slbd;
            vlag = slbd * (1.0 - slbd);
            isbd = ilbd;
            temp = subd * (1.0 - subd);
            if (std::abs(temp) > std::abs(vlag)) {
                step = subd;
                vlag = temp;
                isbd = iubd;
            }
            if (subd > 0.5) {
                if (std::abs(vlag) < .25) {
                    step = 0.5;
                    vlag = .25;
                    isbd = 0;
                }
            }
            vlag *= dderiv;
        }

        /*     Calculate PREDSQ for the current line search and maintain PRESAV. */

        temp = step * (1.0 - step) * distsq;
        predsq = vlag * vlag * (vlag * vlag + ha * temp * temp);
        if (predsq > presav) {
            presav = predsq;
            ksav = k;
            stpsav = step;
            ibdsav = isbd;
        }
L80:
        ;
    }

    /*     Construct XNEW in a way that satisfies the bound constraints exactly. */

    for (long i = 1; i <= n; ++i) {
        temp = xopt[i] + stpsav * (xpt[ksav + i * xpt_dim1] - xopt[i]);
        xnew[i] = std::max(sl[i], std::min(su[i], temp));
    }
    if (ibdsav < 0) {
        xnew[-ibdsav] = sl[-ibdsav];
    }
    if (ibdsav > 0) {
        xnew[ibdsav] = su[ibdsav];
    }

    /*     Prepare for the iterative method that assembles the constrained Cauchy */
    /*     step in W. The sum of squares of the fixed components of W is formed in */
    /*     WFIXSQ, and the free components of W are set to BIGSTP. */

    bigstp = adelt + adelt;
    iflag = 0;
L100:
    wfixsq = 0.0;
    ggfree = 0.0;
    for (long i = 1; i <= n; ++i) {
        w[i] = 0.0;
        tempa = std::min(xopt[i] - sl[i], glag[i]);
        tempb = std::max(xopt[i] - su[i], glag[i]);
        if (tempa > 0.0 || tempb < 0.0) {
            w[i] = bigstp;
            ggfree += square(glag[i]);
        }
    }
    if (ggfree == 0.0) {
        cauchy = 0.0;
        goto L200;
    }

    /*     Investigate whether more components of W can be fixed. */

L120:
    temp = adelt * adelt - wfixsq;
    if (temp > 0.0) {
        wsqsav = wfixsq;
        step = std::sqrt(temp / ggfree);
        ggfree = 0.0;
        for (long i = 1; i <= n; ++i) {
            if (w[i] == bigstp) {
                temp = xopt[i] - step * glag[i];
                if (temp <= sl[i]) {
                    w[i] = sl[i] - xopt[i];
                    wfixsq += square(w[i]);
                } else if (temp >= su[i]) {
                    w[i] = su[i] - xopt[i];
                    wfixsq += square(w[i]);
                } else {
                    ggfree += square(glag[i]);
                }
            }
        }
        if (wfixsq > wsqsav && ggfree > 0.0) {
            goto L120;
        }
    }

    /*     Set the remaining free components of W and all components of XALT, */
    /*     except that W may be scaled later. */

    gw = 0.0;
    for (long i = 1; i <= n; ++i) {
        if (w[i] == bigstp) {
            w[i] = -step * glag[i];
            xalt[i] = std::max(sl[i], std::min(su[i], xopt[i] + w[i]));
        } else if (w[i] == 0.0) {
            xalt[i] = xopt[i];
        } else if (glag[i] > 0.0) {
            xalt[i] = sl[i];
        } else {
            xalt[i] = su[i];
        }
        gw += glag[i] * w[i];
    }

    /*     Set CURV to the curvature of the KNEW-th Lagrange function along W. */
    /*     Scale W by a factor less than one if that can reduce the modulus of */
    /*     the Lagrange function at XOPT+W. Set CAUCHY to the final value of */
    /*     the square of this function. */

    curv = 0.0;
    for (long k = 1; k <= npt; ++k) {
        temp = 0.0;
        for (long j = 1; j <= n; ++j) {
            temp += xpt[k + j * xpt_dim1] * w[j];
        }
        curv += hcol[k] * temp * temp;
    }
    if (iflag == 1) {
        curv = -curv;
    }
    if (curv > -gw && curv < -one_plus_sqrt_2 * gw) {
        scale = -gw / curv;
        for (long i = 1; i <= n; ++i) {
            temp = xopt[i] + scale * w[i];
            xalt[i] = std::max(sl[i], std::min(su[i], temp));
        }
        cauchy = square(0.5 * gw * scale);
    } else {
        cauchy = square(gw + 0.5 * curv);
    }

    /*     If IFLAG is zero, then XALT is calculated as before after reversing */
    /*     the sign of GLAG. Thus two XALT vectors become available. The one that */
    /*     is chosen is the one that gives the larger value of CAUCHY. */

    if (iflag == 0) {
        for (long i = 1; i <= n; ++i) {
            glag[i] = -glag[i];
            w[n + i] = xalt[i];
        }
        csave = cauchy;
        iflag = 1;
        goto L100;
    }
    if (csave > cauchy) {
        for (long i = 1; i <= n; ++i) {
            xalt[i] = w[n + i];
        }
        cauchy = csave;
    }
L200:
    ;
}

} // namespace bobyqa_detail
