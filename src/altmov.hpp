#pragma once

namespace bobyqa_detail {

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
);

} // namespace bobyqa_detail
