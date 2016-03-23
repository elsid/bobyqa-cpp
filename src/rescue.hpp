#pragma once

#include "update.hpp"
#include "utils.hpp"

namespace bobyqa_detail {

template <class Function>
void rescue(
    const Function &function,
    const long n,
    const long npt,
    const double *const xl,
    const double *const xu,
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
    long& nf,
    const double delta,
    long& kopt,
    double *const vlag,
    double *const ptsaux,
    double *const ptsid,
    double *const w
)
{
    /* Local variables */
    long ih, jp, iq, iw;
    double xp = 0, xq = 0, den = 0;
    long ihp = 0;
    long ihq, jpn;
    double sum = 0, diff = 0, beta = 0;
    long kold;
    double winc;
    long nrem, knew;
    double temp, bsum;
    double hdiag = 0, fbase = 0, denom = 0, vquad = 0, sumpq = 0;
    double dsqmin, distsq, vlmxsq;



    /*     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT, */
    /*       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as */
    /*       the corresponding arguments of BOBYQB on the entry to RESCUE. */
    /*     NF is maintained as the number of calls of CALFUN so far, except that */
    /*       NF is set to -1 if the value of MAXFUN prevents further progress. */
    /*     KOPT is maintained so that FVAL(KOPT) is the least calculated function */
    /*       value. Its correct value must be given on entry. It is updated if a */
    /*       new least function value is found, but the corresponding changes to */
    /*       XOPT and GOPT have to be made later by the calling program. */
    /*     DELTA is the current trust region radius. */
    /*     VLAG is a working space vector that will be used for the values of the */
    /*       provisional Lagrange functions at each of the interpolation points. */
    /*       They are part of a product that requires VLAG to be of length NDIM. */
    /*     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and */
    /*       PTSAUX(2,J) specify the two positions of provisional interpolation */
    /*       points when a nonzero step is taken along e_J (the J-th coordinate */
    /*       direction) through XBASE+XOPT, as specified below. Usually these */
    /*       steps have length DELTA, but other lengths are chosen if necessary */
    /*       in order to satisfy the given bounds on the variables. */
    /*     PTSID is also a working space array. It has NPT components that denote */
    /*       provisional new positions of the original interpolation points, in */
    /*       case changes are needed to restore the linear independence of the */
    /*       interpolation conditions. The K-th point is a candidate for change */
    /*       if and only if PTSID(K) is nonzero. In this case let p and q be the */
    /*       long parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p */
    /*       and q are both positive, the step from XBASE+XOPT to the new K-th */
    /*       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise */
    /*       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or */
    /*       p=0, respectively. */
    /*     The first NDIM+NPT elements of the array W are used for working space. */
    /*     The final elements of BMAT and ZMAT are set in a well-conditioned way */
    /*       to the values that are appropriate for the new interpolation points. */
    /*     The elements of GOPT, HQ and PQ are also revised to the values that are */
    /*       appropriate to the final quadratic model. */

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
    const double sfrac = 0.5 / double(np);
    const long nptm = npt - np;

    /*     Shift the interpolation points so that XOPT becomes the origin, and set */
    /*     the elements of ZMAT to zero. The value of SUMPQ is required in the */
    /*     updating of HQ below. The squares of the distances from XOPT to the */
    /*     other interpolation points are set at the end of W. Increments of WINC */
    /*     may be added later to these squares to balance the consideration of */
    /*     the choice of point that is going to become current. */

    sumpq = 0.0;
    winc = 0.0;
    for (long k = 1; k <= npt; ++k) {
        distsq = 0.0;
        for (long j = 1; j <= n; ++j) {
            xpt[k + j * xpt_dim1] -= xopt[j];
            distsq += square(xpt[k + j * xpt_dim1]);
        }
        sumpq += pq[k];
        w[ndim + k] = distsq;
        winc = std::max(winc,distsq);
        for (long j = 1; j <= nptm; ++j) {
            zmat[k + j * zmat_dim1] = 0.0;
        }
    }

    /*     Update HQ so that HQ and PQ define the second derivatives of the model */
    /*     after XBASE has been shifted to the trust region centre. */

    ih = 0;
    for (long j = 1; j <= n; ++j) {
        w[j] = 0.5 * sumpq * xopt[j];
        for (long k = 1; k <= npt; ++k) {
            w[j] += pq[k] * xpt[k + j * xpt_dim1];
        }
        for (long i__ = 1; i__ <= j; ++i__) {
            ++ih;
            hq[ih] = hq[ih] + w[i__] * xopt[j] + w[j] * xopt[i__];
        }
    }

    /*     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and */
    /*     also set the elements of PTSAUX. */

    for (long j = 1; j <= n; ++j) {
        xbase[j] += xopt[j];
        sl[j] -= xopt[j];
        su[j] -= xopt[j];
        xopt[j] = 0.0;
        ptsaux[(j << 1) + 1] = std::min(delta, su[j]);
        ptsaux[(j << 1) + 2] = std::max(-delta, sl[j]);
        if (ptsaux[(j << 1) + 1] + ptsaux[(j << 1) + 2] < 0.0) {
            temp = ptsaux[(j << 1) + 1];
            ptsaux[(j << 1) + 1] = ptsaux[(j << 1) + 2];
            ptsaux[(j << 1) + 2] = temp;
        }
        if (std::abs(ptsaux[(j << 1) + 2]) < 0.5 * std::abs(ptsaux[(j << 1) + 1])) {
            ptsaux[(j << 1) + 2] = 0.5 * ptsaux[(j << 1) + 1];
        }
        for (long i__ = 1; i__ <= ndim; ++i__) {
            bmat[i__ + j * bmat_dim1] = 0.0;
        }
    }
    fbase = fval[kopt];

    /*     Set the identifiers of the artificial interpolation points that are */
    /*     along a coordinate direction from XOPT, and set the corresponding */
    /*     nonzero elements of BMAT and ZMAT. */

    ptsid[1] = sfrac;
    for (long j = 1; j <= n; ++j) {
        jp = j + 1;
        jpn = jp + n;
        ptsid[jp] = double(j + sfrac);
        if (jpn <= npt) {
            ptsid[jpn] = double(j) / double(np) + sfrac;
            temp = 1.0 / (ptsaux[(j << 1) + 1] - ptsaux[(j << 1) + 2]);
            bmat[jp + j * bmat_dim1] = -temp + 1.0 / ptsaux[(j << 1) + 1];
            bmat[jpn + j * bmat_dim1] = temp + 1.0 / ptsaux[(j << 1) + 2];
            bmat[j * bmat_dim1 + 1] = -bmat[jp + j * bmat_dim1] - bmat[jpn + j * bmat_dim1];
            zmat[j * zmat_dim1 + 1] = sqrt_2 / std::abs(ptsaux[(j << 1) + 1] * ptsaux[(j << 1) + 2]);
            zmat[jp + j * zmat_dim1] = zmat[j * zmat_dim1 + 1] * ptsaux[(j << 1) + 2] * temp;
            zmat[jpn + j * zmat_dim1] = -zmat[j * zmat_dim1 + 1] * ptsaux[(j << 1) + 1] * temp;
        } else {
            bmat[j * bmat_dim1 + 1] = -1.0 / ptsaux[(j << 1) + 1];
            bmat[jp + j * bmat_dim1] = 1.0 / ptsaux[(j << 1) + 1];
            bmat[j + npt + j * bmat_dim1] = -0.5 * square(ptsaux[(j << 1) + 1]);
        }
    }

    /*     Set any remaining identifiers with their nonzero elements of ZMAT. */

    if (npt >= n + np) {
        for (long k = np << 1; k <= npt; ++k) {
            iw = long((double(k - np) - 0.5) / double(n)
            );
            const long ip = k - np - iw * n;
            iq = ip + iw;
            if (iq > n) {
                iq -= n;
            }
            ptsid[k] = double(ip) + double(iq) / double(np) + sfrac;
            temp = 1.0 / (ptsaux[(ip << 1) + 1] * ptsaux[(iq << 1) + 1]);
            zmat[(k - np) * zmat_dim1 + 1] = temp;
            zmat[ip + 1 + (k - np) * zmat_dim1] = -temp;
            zmat[iq + 1 + (k - np) * zmat_dim1] = -temp;
            zmat[k + (k - np) * zmat_dim1] = temp;
        }
    }
    nrem = npt;
    kold = 1;
    knew = kopt;

    /*     Reorder the provisional points in the way that exchanges PTSID(KOLD) */
    /*     with PTSID(KNEW). */

L80:
    for (long j = 1; j <= n; ++j) {
        temp = bmat[kold + j * bmat_dim1];
        bmat[kold + j * bmat_dim1] = bmat[knew + j * bmat_dim1];
        bmat[knew + j * bmat_dim1] = temp;
    }
    for (long j = 1; j <= nptm; ++j) {
        temp = zmat[kold + j * zmat_dim1];
        zmat[kold + j * zmat_dim1] = zmat[knew + j * zmat_dim1];
        zmat[knew + j * zmat_dim1] = temp;
    }
    ptsid[kold] = ptsid[knew];
    ptsid[knew] = 0.0;
    w[ndim + knew] = 0.0;
    --nrem;
    if (knew != kopt) {
        temp = vlag[kold];
        vlag[kold] = vlag[knew];
        vlag[knew] = temp;

        /*     Update the BMAT and ZMAT matrices so that the status of the KNEW-th */
        /*     interpolation point can be changed from provisional to original. The */
        /*     branch to label 350 occurs if all the original points are reinstated. */
        /*     The nonnegative values of W(NDIM+K) are required in the search below. */

        update(n, npt, bmat + bmat_offset, zmat + zmat_offset, ndim, vlag, beta, denom, knew, w);
        if (nrem == 0) {
            goto L350;
        }
        for (long k = 1; k <= npt; ++k) {
            w[ndim + k] = std::abs(w[ndim + k]);
        }
    }

    /*     Pick the index KNEW of an original interpolation point that has not */
    /*     yet replaced one of the provisional interpolation points, giving */
    /*     attention to the closeness to XOPT and to previous tries with KNEW. */

L120:
    dsqmin = 0.0;
    for (long k = 1; k <= npt; ++k) {
        if (w[ndim + k] > 0.0) {
            if (dsqmin == 0.0 || w[ndim + k] < dsqmin) {
                knew = k;
                dsqmin = w[ndim + k];
            }
        }
    }
    if (dsqmin == 0.0) {
        goto L260;
    }

    /*     Form the W-vector of the chosen original interpolation point. */

    for (long j = 1; j <= n; ++j) {
        w[npt + j] = xpt[knew + j * xpt_dim1];
    }
    for (long k = 1; k <= npt; ++k) {
        sum = 0.0;
        if (k == kopt) {
        } else if (ptsid[k] == 0.0) {
            for (long j = 1; j <= n; ++j) {
                sum += w[npt + j] * xpt[k + j * xpt_dim1];
            }
        } else {
            const long ip = long(ptsid[k]);
            if (ip > 0) {
                sum = w[npt + ip] * ptsaux[(ip << 1) + 1];
            }
            iq = long(double(np) * ptsid[k] - double(ip * np));
            if (iq > 0) {
                iw = 1;
                if (ip == 0) {
                    iw = 2;
                }
                sum += w[npt + iq] * ptsaux[iw + (iq << 1)];
            }
        }
        w[k] = 0.5 * sum * sum;
    }

    /*     Calculate VLAG and BETA for the required updating of the H matrix if */
    /*     XPT(KNEW,.) is reinstated in the set of interpolation points. */

    for (long k = 1; k <= npt; ++k) {
        sum = 0.0;
        for (long j = 1; j <= n; ++j) {
            sum += bmat[k + j * bmat_dim1] * w[npt + j];
        }
        vlag[k] = sum;
    }
    beta = 0.0;
    for (long j = 1; j <= nptm; ++j) {
        sum = 0.0;
        for (long k = 1; k <= npt; ++k) {
            sum += zmat[k + j * zmat_dim1] * w[k];
        }
        beta -= sum * sum;
        for (long k = 1; k <= npt; ++k) {
            vlag[k] += sum * zmat[k + j * zmat_dim1];
        }
    }
    bsum = 0.0;
    distsq = 0.0;
    for (long j = 1; j <= n; ++j) {
        sum = 0.0;
        for (long k = 1; k <= npt; ++k) {
            sum += bmat[k + j * bmat_dim1] * w[k];
        }
        jp = j + npt;
        bsum += sum * w[jp];
        for (long ip = npt + 1; ip <= ndim; ++ip) {
            sum += bmat[ip + j * bmat_dim1] * w[ip];
        }
        bsum += sum * w[jp];
        vlag[jp] = sum;
        distsq += square(xpt[knew + j * xpt_dim1]);
    }
    beta = 0.5 * distsq * distsq + beta - bsum;
    vlag[kopt] += 1.0;

    /*     KOLD is set to the index of the provisional interpolation point that is */
    /*     going to be deleted to make way for the KNEW-th original interpolation */
    /*     point. The choice of KOLD is governed by the avoidance of a small value */
    /*     of the denominator in the updating calculation of UPDATE. */

    denom = 0.0;
    vlmxsq = 0.0;
    for (long k = 1; k <= npt; ++k) {
        if (ptsid[k] != 0.0) {
            hdiag = 0.0;
            for (long j = 1; j <= nptm; ++j) {
                hdiag += square(zmat[k + j * zmat_dim1]);
            }
            den = beta * hdiag + square(vlag[k]);
            if (den > denom) {
                kold = k;
                denom = den;
            }
        }
        vlmxsq = std::max(vlmxsq, square(vlag[k]));
    }
    if (denom <= vlmxsq * .01) {
        w[ndim + knew] = -w[ndim + knew] - winc;
        goto L120;
    }
    goto L80;

    /*     When label 260 is reached, all the final positions of the interpolation */
    /*     points have been chosen although any changes have not been included yet */
    /*     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart */
    /*     from the shift of XBASE, the updating of the quadratic model remains to */
    /*     be done. The following cycle through the new interpolation points begins */
    /*     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero, */
    /*     except that a RETURN occurs if MAXFUN prohibits another value of F. */

L260:
    for (long kpt = 1; kpt <= npt; ++kpt) {
        if (ptsid[kpt] == 0.0) {
            goto L340;
        }
        if (nf >= maxfun) {
            nf = -1;
            goto L350;
        }
        ih = 0;
        for (long j = 1; j <= n; ++j) {
            w[j] = xpt[kpt + j * xpt_dim1];
            xpt[kpt + j * xpt_dim1] = 0.0;
            temp = pq[kpt] * w[j];
            for (long i__ = 1; i__ <= j; ++i__) {
                ++ih;
                hq[ih] += temp * w[i__];
            }
        }
        pq[kpt] = 0.0;

        {

        const long ip = long(ptsid[kpt]);
        iq = long(double(np) * ptsid[kpt] - double(ip * np));
        if (ip > 0) {
            xp = ptsaux[(ip << 1) + 1];
            xpt[kpt + ip * xpt_dim1] = xp;
        }
        if (iq > 0) {
            xq = ptsaux[(iq << 1) + 1];
            if (ip == 0) {
                xq = ptsaux[(iq << 1) + 2];
            }
            xpt[kpt + iq * xpt_dim1] = xq;
        }

        /*     Set VQUAD to the value of the current model at the new point. */

        vquad = fbase;
        if (ip > 0) {
            ihp = (ip + ip * ip) / 2;
            vquad += xp * (gopt[ip] + 0.5 * xp * hq[ihp]);
        }
        if (iq > 0) {
            ihq = (iq + iq * iq) / 2;
            vquad += xq * (gopt[iq] + 0.5 * xq * hq[ihq]);
            if (ip > 0) {
                iw = std::max(ihp, ihq) - std::abs(ip - iq);
                vquad += xp * xq * hq[iw];
            }
        }
        for (long k = 1; k <= npt; ++k) {
            temp = 0.0;
            if (ip > 0) {
                temp += xp * xpt[k + ip * xpt_dim1];
            }
            if (iq > 0) {
                temp += xq * xpt[k + iq * xpt_dim1];
            }
            vquad += 0.5 * pq[k] * temp * temp;
        }

        }

        /*     Calculate F at the new interpolation point, and set DIFF to the factor */
        /*     that is going to multiply the KPT-th Lagrange function when the model */
        /*     is updated to provide interpolation to the new function value. */

        for (long i__ = 1; i__ <= n; ++i__) {
            w[i__] = std::min(std::max(xl[i__], xbase[i__] + xpt[kpt + i__ * xpt_dim1]), xu[i__]);
            if (xpt[kpt + i__ * xpt_dim1] == sl[i__]) {
                w[i__] = xl[i__];
            }
            if (xpt[kpt + i__ * xpt_dim1] == su[i__]) {
                w[i__] = xu[i__];
            }
        }
        ++(nf);
        {
            const double f = function(n, w + 1);
            fval[kpt] = f;
            if (f < fval[kopt]) {
                kopt = kpt;
            }
            diff = f - vquad;
        }

        /*     Update the quadratic model. The RETURN from the subroutine occurs when */
        /*     all the new interpolation points are included in the model. */

        for (long i__ = 1; i__ <= n; ++i__) {
            gopt[i__] += diff * bmat[kpt + i__ * bmat_dim1];
        }
        for (long k = 1; k <= npt; ++k) {
            sum = 0.0;
            for (long j = 1; j <= nptm; ++j) {
                sum += zmat[k + j * zmat_dim1] * zmat[kpt + j * zmat_dim1];
            }
            temp = diff * sum;
            if (ptsid[k] == 0.0) {
                pq[k] += temp;
            } else {
                const long ip = long(ptsid[k]);
                iq = long(double(np) * ptsid[k] - double(ip * np));
                ihq = (iq * iq + iq) / 2;
                if (ip == 0) {
                    hq[ihq] += temp * square(ptsaux[(iq << 1) + 2]);
                } else {
                    ihp = (ip * ip + ip) / 2;
                    hq[ihp] += temp * square(ptsaux[(ip << 1) + 1]);
                    if (iq > 0) {
                        hq[ihq] += temp * square(ptsaux[(iq << 1) + 1]);
                        iw = std::max(ihp, ihq) - std::abs(iq - ip);
                        hq[iw] += temp * ptsaux[(ip << 1) + 1] * ptsaux[(iq << 1) + 1];
                    }
                }
            }
        }
        ptsid[kpt] = 0.0;
L340:
        ;
    }
L350:
    ;
}

} // namespace bobyqa_detail
