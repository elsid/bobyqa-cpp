#include "impl.hpp"

#include <bobyqa.h>

double bobyqa(
    const BobyqaFunction function,
    const long n,
    const long npt,
    double *x,
    const double *xl,
    const double *xu,
    const double rhobeg,
    const double rhoend,
    const long maxfun,
    double *w
) {
    return bobyqa_detail::impl([=] (long n, const double *x) -> double {
        return function(n, x);
    }, n, npt, x, xl, xu, rhobeg, rhoend, maxfun, w);
}

double bobyqa_closure(
    BobyqaClosure *const closure,
    const long n,
    const long npt,
    double *x,
    const double *xl,
    const double *xu,
    const double rhobeg,
    const double rhoend,
    const long maxfun,
    double *w
) {
    return bobyqa_detail::impl([&] (long n, const double *x) -> double {
        return closure->function(closure->data, n, x);
    }, n, npt, x, xl, xu, rhobeg, rhoend, maxfun, w);
}

double bobyqa_closure_const(
    const BobyqaClosureConst *const closure,
    const long n,
    const long npt,
    double *x,
    const double *xl,
    const double *xu,
    const double rhobeg,
    const double rhoend,
    const long maxfun,
    double *w
) {
    return bobyqa_detail::impl([&] (long n, const double *x) -> double {
        return closure->function(closure->data, n, x);
    }, n, npt, x, xl, xu, rhobeg, rhoend, maxfun, w);
}
