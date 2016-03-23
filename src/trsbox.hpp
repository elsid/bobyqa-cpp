#pragma once

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
    double *const d__,
    double *const gnew,
    double *const xbdi,
    double *const s,
    double *const hs,
    double *const hred,
    double *const dsq,
    double *const crvmin
);

} // namespace bobyqa_detail
