#include <algorithm>
#include <cmath>

namespace bobyqa_detail {

double less_abs(double lhs, double rhs) {
    return std::abs(lhs) < std::abs(rhs);
}

void update(
    const long n,
    const long npt,
    double *bmat,
    double *zmat,
    const long ndim,
    double *const vlag,
    const double beta,
    const double denom,
    const long knew,
    double *const w
) {
    /*     The arrays BMAT and ZMAT are updated, as required by the new position */
    /*     of the interpolation point that has the index KNEW. The vector VLAG has */
    /*     N+NPT components, set on entry to the first NPT and last N components */
    /*     of the product Hw in equation (4.11) of the Powell (2006) paper on */
    /*     NEWUOA. Further, BETA is set on entry to the value of the parameter */
    /*     with that name, and DENOM is set to the denominator of the updating */
    /*     formula. Elements of ZMAT may be treated as zero if their moduli are */
    /*     at most ZTEST. The first NDIM elements of W are used for working space. */

    /*     Set some constants. */

    /* Parameter adjustments */
    const long zmat_dim1 = npt;
    const long zmat_offset = 1 + zmat_dim1;
    zmat -= zmat_offset;
    const long bmat_dim1 = ndim;
    const long bmat_offset = 1 + bmat_dim1;
    bmat -= bmat_offset;

    /* Function Body */
    const long nptm = npt - n - 1;
    const auto zmat_end = zmat + zmat_offset + nptm * npt;
    const auto zmat_max = std::max_element(zmat + zmat_offset, zmat_end, less_abs);
    const double ztest = zmat_max == zmat_end ? 0 : *zmat_max * 1e-20;

    /*     Apply the rotations that put zeros in the KNEW-th row of ZMAT. */

    for (long j = 2; j <= nptm; ++j) {
        if (std::abs(zmat[knew + j * zmat_dim1]) > ztest) {
            const double temp = std::hypot(zmat[knew + zmat_dim1], zmat[knew + j * zmat_dim1]);
            const double tempa = zmat[knew + zmat_dim1] / temp;
            const double tempb = zmat[knew + j * zmat_dim1] / temp;
            for (long i = 1; i <= npt; ++i) {
                const double temp = tempa * zmat[i + zmat_dim1] + tempb * zmat[i + j * zmat_dim1];
                zmat[i + j * zmat_dim1] = tempa * zmat[i + j * zmat_dim1] - tempb * zmat[i + zmat_dim1];
                zmat[i + zmat_dim1] = temp;
            }
        }
        zmat[knew + j * zmat_dim1] = 0.0;
    }

    /*     Put the first NPT components of the KNEW-th column of HLAG into W, */
    /*     and calculate the parameters of the updating formula. */

    for (long i = 1; i <= npt; ++i) {
        w[i] = zmat[knew + zmat_dim1] * zmat[i + zmat_dim1];
    }
    const double alpha = w[knew];
    const double tau = vlag[knew];
    vlag[knew] -= 1.0;

    /*     Complete the updating of ZMAT. */
    const double temp = std::sqrt(denom);
    const double tempb = zmat[knew + zmat_dim1] / temp;
    const double tempa = tau / temp;
    for (long i = 1; i <= npt; ++i) {
        zmat[i + zmat_dim1] = tempa * zmat[i + zmat_dim1] - tempb * vlag[i];
    }

    /*     Finally, update the matrix BMAT. */

    for (long j = 1; j <= n; ++j) {
        const long jp = npt + j;
        w[jp] = bmat[knew + j * bmat_dim1];
        const double tempa = (alpha * vlag[jp] - tau * w[jp]) / denom;
        const double tempb = (-(beta) * w[jp] - tau * vlag[jp]) / denom;
        for (long i = 1; i <= jp; ++i) {
            bmat[i + j * bmat_dim1] = bmat[i + j * bmat_dim1] + tempa * vlag[i] + tempb * w[i];
            if (i > npt) {
                bmat[jp + (i - npt) * bmat_dim1] = bmat[i + j * bmat_dim1];
            }
        }
    }
}

} // namespace bobyqa_detail
