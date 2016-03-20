#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <bobyqa.h>

static void print_array(const size_t size, const double *values) {
    printf("[%lu]", size);
    for (size_t i = 0; i < size; ++i) {
        printf(" %lf", values[i]);
    }
}

static size_t function_calls_count = 0;

static double function(long n, const double *x) {
    assert(n == 2);
    (void) n;
    ++function_calls_count;
    return -4*x[0]*x[1] + 5*x[0]*x[0] + 8*x[1]*x[1] + 16*sqrt(5.0)*x[0] + 8*sqrt(5.0)*x[1] - 44.0;
}

int main() {
    const int variables_count = 2;
//    const int number_of_interpolation_conditions = variables_count + 2;
    const int number_of_interpolation_conditions = (variables_count + 1)*(variables_count + 2)/2;
    double variables_values[] = {0.0, -sqrt(5.0)};
    const double lower_bound[] = {-4.0, -3.0};
    const double upper_bound[] = {5.0, 5.0};
    const double initial_trust_region_radius = 1e-3;
    const double final_trust_region_radius = 1e3;
    const int max_function_calls_count = 100;
    const size_t working_space_size = BOBYQA_WORKING_SPACE_SIZE(variables_count,
                                                                number_of_interpolation_conditions);
    double working_space[working_space_size];

    printf("initial: ");
    print_array(variables_count, variables_values);
    printf("\n");

    const double result = bobyqa(
            function,
            variables_count,
            number_of_interpolation_conditions,
            variables_values,
            lower_bound,
            upper_bound,
            initial_trust_region_radius,
            final_trust_region_radius,
            max_function_calls_count,
            working_space);

    printf("final: ");
    print_array(variables_count, variables_values);
    printf("\n");

    printf("function_calls_count: %ld\n", function_calls_count);

    printf("result: %lf\n", result);

    return 0;
}
