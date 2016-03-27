#include "impl.hpp"

#include <bobyqa.h>

double bobyqa(
    BobyqaFunction function,
    long n,
    long npt,
    double *x,
    const double *xl,
    const double *xu,
    double rhobeg,
    double rhoend,
    long maxfun,
    double *w
) {
    return bobyqa_detail::impl([=] (long n, const double *x) -> double {
        return function(n, x);
    }, n, npt, x, xl, xu, rhobeg, rhoend, maxfun, w);
}

double bobyqa_closure(
    BobyqaClosure *closure,
    long n,
    long npt,
    double *x,
    const double *xl,
    const double *xu,
    double rhobeg,
    double rhoend,
    long maxfun,
    double *w
) {
    return bobyqa_detail::impl([=] (long n, const double *x) -> double {
        return closure->function(closure->data, n, x);
    }, n, npt, x, xl, xu, rhobeg, rhoend, maxfun, w);
}

double bobyqa_closure_const(
    BobyqaClosureConst *closure,
    long n,
    long npt,
    double *x,
    const double *xl,
    const double *xu,
    double rhobeg,
    double rhoend,
    long maxfun,
    double *w
) {
    return bobyqa_detail::impl([=] (long n, const double *x) -> double {
        return closure->function(closure->data, n, x);
    }, n, npt, x, xl, xu, rhobeg, rhoend, maxfun, w);
}
