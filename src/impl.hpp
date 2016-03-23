#pragma once

#include "bobyqb.hpp"

namespace bobyqa_detail {

template <class Function>
double impl(const Function &function, long n, long npt, double *x,
        const double *xl, const double *xu, double rhobeg, double rhoend,
        long maxfun, double *w) {
    /*     This subroutine seeks the least value of a function of many variables, */
    /*     by applying a trust region method that forms quadratic models by */
    /*     interpolation. There is usually some freedom in the interpolation */
    /*     conditions, which is taken up by minimizing the Frobenius norm of */
    /*     the change to the second derivative of the model, beginning with the */
    /*     zero matrix. The values of the variables are constrained by upper and */
    /*     lower bounds. The arguments of the subroutine are as follows. */

    /*     N must be set to the number of variables and must be at least two. */
    /*     NPT is the number of interpolation conditions. Its value must be in */
    /*       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not */
    /*       recommended. */
    /*     Initial values of the variables must be set in X(1),X(2),...,X(N). They */
    /*       will be changed to the values that give the least calculated F. */
    /*     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper */
    /*       bounds, respectively, on X(I). The construction of quadratic models */
    /*       requires XL(I) to be strictly less than XU(I) for each I. Further, */
    /*       the contribution to a model from changes to the I-th variable is */
    /*       damaged severely by rounding errors if XU(I)-XL(I) is too small. */
    /*     RHOBEG and RHOEND must be set to the initial and final values of a trust */
    /*       region radius, so both must be positive with RHOEND no greater than */
    /*       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest */
    /*       expected change to a variable, while RHOEND should indicate the */
    /*       accuracy that is required in the final values of the variables. An */
    /*       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N, */
    /*       is less than 2*RHOBEG. */
    /*     MAXFUN must be set to an upper bound on the number of calls of CALFUN. */
    /*     The array W will be used for working space. Its length must be at least */
    /*       (NPT+5)*(NPT+N)+3*N*(N+5)/2. */

    /* Parameter adjustments */
    --w;
    --xu;
    --xl;
    --x;

    /* Function Body */
    const long np = n + 1;

    /*     Return if the value of NPT is unacceptable. */
    if (npt < n + 2 || npt > (n + 2) * np / 2) {
        DEBUG_LOG("Return from BOBYQA because NPT is not in the required interval.\n");
        return 0.0;
    }

    /*     Partition the working space array, so that different parts of it can */
    /*     be treated separately during the calculation of BOBYQB. The partition */
    /*     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the */
    /*     space that is taken by the last array in the argument list of BOBYQB. */

    const long ndim = npt + n;
    const long ixp = 1 + n;
    const long ifv = ixp + n * npt;
    const long ixo = ifv + npt;
    const long igo = ixo + n;
    const long ihq = igo + n;
    const long ipq = ihq + n * np / 2;
    const long ibmat = ipq + npt;
    const long izmat = ibmat + ndim * n;
    const long isl = izmat + npt * (npt - np);
    const long isu = isl + n;
    const long ixn = isu + n;
    const long ixa = ixn + n;
    const long id_ = ixa + n;
    const long ivl = id_ + n;
    const long iw = ivl + ndim;

    /*     Return if there is insufficient space between the bounds. Modify the */
    /*     initial X if necessary in order to avoid conflicts between the bounds */
    /*     and the construction of the first quadratic model. The lower and upper */
    /*     bounds on moves from the updated X are set now, in the ISL and ISU */
    /*     partitions of W, in order to provide useful and exact information about */
    /*     components of X that become within distance RHOBEG from their bounds. */

    for (long j = 1; j <= n; ++j) {
        const double temp = xu[j] - xl[j];
        if (temp < rhobeg + rhobeg) {
            DEBUG_LOG("Return from BOBYQA because one of the differences in x_lower and x_upper "
                      "is less than 2*rho_begin.\n");
            return 0.0;
        }
        const long jsl = isl + j - 1;
        const long jsu = jsl + n;
        w[jsl] = xl[j] - x[j];
        w[jsu] = xu[j] - x[j];
        if (w[jsl] >= -(rhobeg)) {
            if (w[jsl] >= 0.0) {
                x[j] = xl[j];
                w[jsl] = 0.0;
                w[jsu] = temp;
            } else {
                x[j] = xl[j] + rhobeg;
                w[jsl] = -(rhobeg);
                w[jsu] = std::max(xu[j] - x[j], rhobeg);
            }
        } else if (w[jsu] <= rhobeg) {
            if (w[jsu] <= 0.0) {
                x[j] = xu[j];
                w[jsl] = -temp;
                w[jsu] = 0.0;
            } else {
                x[j] = xu[j] - rhobeg;
                w[jsl] = std::min(xl[j] - x[j], -rhobeg);
                w[jsu] = rhobeg;
            }
        }
    }

    /*     Make the call of BOBYQB. */

    return bobyqb(function, n, npt, x, xl, xu, rhobeg, rhoend, maxfun, w,
            w + ixp, w + ifv - 1, w + ixo - 1, w + igo - 1, w + ihq - 1, w + ipq - 1, w + ibmat,
            w + izmat, ndim, w + isl - 1, w + isu - 1, w + ixn - 1, w + ixa - 1, w + id_ - 1,
            w + ivl - 1, w + iw - 1);
}

} // namespace bobyqa_detail
