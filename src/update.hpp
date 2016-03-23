#pragma once

namespace bobyqa_detail {

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
);

} // namespace bobyqa_detail
