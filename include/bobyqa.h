#ifndef BOBYQA_H
#define BOBYQA_H

typedef double (*BobyqaFunction)(long n, const double *x);
typedef double (*BobyqaClosureFunction)(void *data, long n, const double *x);

typedef struct {
    void *data;
    BobyqaClosureFunction function;
} BobyqaClosure;

#ifdef __cplusplus
extern "C" {
#endif

/* This subroutine seeks the least value of a function of many variables,
 * by applying a trust region method that forms quadratic models by
 * interpolation. There is usually some freedom in the interpolation
 * conditions, which is taken up by minimizing the Frobenius norm of
 * the change to the second derivative of the model, beginning with the
 * zero matrix. The values of the variables are constrained by upper and
 * lower bounds. The arguments of the subroutine are as follows. */

/* N must be set to the number of variables and must be at least two.
 * NPT is the number of interpolation conditions. Its value must be in
 *   the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
 *   recommended.
 * Initial values of the variables must be set in X(1),X(2),...,X(N). They
 *   will be changed to the values that give the least calculated F.
 * For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
 *   bounds, respectively, on X(I). The construction of quadratic models
 *   requires XL(I) to be strictly less than XU(I) for each I. Further,
 *   the contribution to a model from changes to the I-th variable is
 *   damaged severely by rounding errors if XU(I)-XL(I) is too small.
 * RHOBEG and RHOEND must be set to the initial and final values of a trust
 *   region radius, so both must be positive with RHOEND no greater than
 *   RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
 *   expected change to a variable, while RHOEND should indicate the
 *   accuracy that is required in the final values of the variables. An
 *   error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
 *   is less than 2*RHOBEG.
 * MAXFUN must be set to an upper bound on the number of calls of CALFUN.
 * The array W will be used for working space. Its length must be at least
 *   (NPT+5)*(NPT+N)+3*N*(N+5)/2. */

double bobyqa(BobyqaFunction function, long n, long npt, double *x,
    const double *xl, const double *xu, double rhobeg, double rhoend,
    long maxfun, double *w);

double bobyqa_closure(BobyqaClosure *closure, long n, long npt, double *x,
    const double *xl, const double *xu, double rhobeg, double rhoend,
    long maxfun, double *w);

#ifdef __cplusplus
}
#endif

#define BOBYQA_WORKING_SPACE_SIZE(n, npt) \
    (npt + 5)*(npt + n) + 3*n*(n + 5)/2

#endif // BOBYQA_H
