#define CATCH_CONFIG_MAIN

#include <array>
#include <random>

#include "catch.hpp"

#include <bobyqa.h>

namespace std {

template <std::size_t N>
std::ostream& operator << (std::ostream& stream, const std::array<double, N>& values) {
    const auto precision = stream.precision();
    stream.precision(100);
    stream << "{";
    for (auto& value : values) {
        stream << value << ", ";
    }
    stream << "}";
    stream.precision(precision);
    return stream;
}

}

TEST_CASE("Const function", "[bobyqa]") {
    struct Function {
        static double call(long n, const double *) {
            REQUIRE(n == 2);
            return 0;
        }
    };
    const long variables_count = 2;
    const long number_of_interpolation_conditions = variables_count + 2;
    using Values = std::array<double, variables_count>;
    Values variables_values = {{0, 0}};
    const Values lower_bound = {{-1, -1}};
    const Values upper_bound = {{1, 1}};
    const double initial_trust_region_radius = 1e-3;
    const double final_trust_region_radius = 1e3;
    const long max_function_calls_count = 5;
    const std::size_t working_space_size = BOBYQA_WORKING_SPACE_SIZE(
        variables_count, number_of_interpolation_conditions);
    using WorkingSpace = std::array<double, working_space_size>;
    WorkingSpace working_space;
    const double result = bobyqa(
        Function::call,
        variables_count,
        number_of_interpolation_conditions,
        variables_values.data(),
        lower_bound.data(),
        upper_bound.data(),
        initial_trust_region_radius,
        final_trust_region_radius,
        max_function_calls_count,
        working_space.data()
    );
    CHECK(result == 0.0);
    CHECK(variables_values == Values({{0.0, 0.0}}));
}

TEST_CASE("Plane function", "[bobyqa]") {
    struct Function {
        static double call(long n, const double *x) {
            REQUIRE(n == 2);
            return x[0] + x[1];
        }
    };
    const long variables_count = 2;
    const long number_of_interpolation_conditions = variables_count + 2;
    using Values = std::array<double, variables_count>;
    Values variables_values = {{0.5, 0.5}};
    const Values lower_bound = {{-1, -1}};
    const Values upper_bound = {{1, 1}};
    const double initial_trust_region_radius = 1e-3;
    const double final_trust_region_radius = 1e3;
    const long max_function_calls_count = 20;
    const std::size_t working_space_size = BOBYQA_WORKING_SPACE_SIZE(
        variables_count, number_of_interpolation_conditions);
    using WorkingSpace = std::array<double, working_space_size>;
    WorkingSpace working_space;
    const double result = bobyqa(
        Function::call,
        variables_count,
        number_of_interpolation_conditions,
        variables_values.data(),
        lower_bound.data(),
        upper_bound.data(),
        initial_trust_region_radius,
        final_trust_region_radius,
        max_function_calls_count,
        working_space.data()
    );
    CHECK(result == -2);
    CHECK(variables_values == Values({{-1, -1}}));
}

TEST_CASE("Simple quadratic function", "[bobyqa]") {
    struct Function {
        static double call(long n, const double *x) {
            REQUIRE(n == 2);
            return x[0] * x[0] + x[1] * x[1];
        }
    };
    const long variables_count = 2;
    const long number_of_interpolation_conditions = variables_count + 2;
    using Values = std::array<double, variables_count>;
    Values variables_values = {{0.5, 0.5}};
    const Values lower_bound = {{-1, -1}};
    const Values upper_bound = {{1, 1}};
    const double initial_trust_region_radius = 1e-3;
    const double final_trust_region_radius = 1e3;
    const long max_function_calls_count = 20;
    const std::size_t working_space_size = BOBYQA_WORKING_SPACE_SIZE(
        variables_count, number_of_interpolation_conditions);
    using WorkingSpace = std::array<double, working_space_size>;
    WorkingSpace working_space;
    const double result = bobyqa(
        Function::call,
        variables_count,
        number_of_interpolation_conditions,
        variables_values.data(),
        lower_bound.data(),
        upper_bound.data(),
        initial_trust_region_radius,
        final_trust_region_radius,
        max_function_calls_count,
        working_space.data()
    );
    CHECK(result == 5.074286393753844992191614549949252310767633389332331717014312744140625e-09);
    CHECK(variables_values == Values({{
        -1.590456929081670613135290892614648328162729740142822265625e-05,
        6.9435805384739763233825637911422745673917233943939208984375e-05
    }}));
}

TEST_CASE("Complex quadratic function", "[bobyqa]") {
    struct Function {
        static double call(long n, const double *x) {
            REQUIRE(n == 2);
            return -4*x[0]*x[1] + 5*x[0]*x[0] + 8*x[1]*x[1] + 16*sqrt(5.0)*x[0] + 8*sqrt(5.0)*x[1] - 44.0;
        }
    };
    const long variables_count = 2;
    const long number_of_interpolation_conditions = (variables_count + 1)*(variables_count + 2)/2;
    using Values = std::array<double, variables_count>;
    Values variables_values = {{0, -sqrt(5.0)}};
    const Values lower_bound = {{-4, -3}};
    const Values upper_bound = {{5, 5}};
    const double initial_trust_region_radius = 1e-3;
    const double final_trust_region_radius = 1e3;
    const long max_function_calls_count = 22;
    const std::size_t working_space_size = BOBYQA_WORKING_SPACE_SIZE(
        variables_count, number_of_interpolation_conditions);
    using WorkingSpace = std::array<double, working_space_size>;
    WorkingSpace working_space;
    const double result = bobyqa(
        Function::call,
        variables_count,
        number_of_interpolation_conditions,
        variables_values.data(),
        lower_bound.data(),
        upper_bound.data(),
        initial_trust_region_radius,
        final_trust_region_radius,
        max_function_calls_count,
        working_space.data()
    );
    CHECK(result == -142.99689437998489438541582785546779632568359375);
    CHECK(variables_values == Values({{-4, -2.118033988750505525189282707287929952144622802734375}}));
}

template <class F>
BobyqaClosureConst make_closure_const(const F &function) {
    struct Wrap {
        static double call(const void *data, long n, const double *values) {
            return reinterpret_cast<const F *>(data)->operator()(n, values);
        }
    };
    return BobyqaClosureConst {&function, &Wrap::call};
}

TEST_CASE("Quadratic function with jump discontinuity", "[bobyqa]") {
    const struct Function {
        Function() {}
        double operator ()(long n, const double *x) const {
            REQUIRE(n == 2);
            return x[1] < 0.5 ? x[0] * x[0] + x[1] * x[1] : x[0] * x[0] + x[1] * x[1] + 10;
        }
    } function;
    const auto closure = make_closure_const(function);
    const long variables_count = 2;
    const long number_of_interpolation_conditions = (variables_count + 1)*(variables_count + 2)/2;
    using Values = std::array<double, variables_count>;
    Values variables_values = {{0.5, 0.5}};
    const Values lower_bound = {{-1, -1}};
    const Values upper_bound = {{1, 1}};
    const double initial_trust_region_radius = 1e-3;
    const double final_trust_region_radius = 1e3;
    const long max_function_calls_count = 8;
    const std::size_t working_space_size = BOBYQA_WORKING_SPACE_SIZE(
        variables_count, number_of_interpolation_conditions);
    using WorkingSpace = std::array<double, working_space_size>;
    WorkingSpace working_space;
    const double result = bobyqa_closure_const(
        &closure,
        variables_count,
        number_of_interpolation_conditions,
        variables_values.data(),
        lower_bound.data(),
        upper_bound.data(),
        initial_trust_region_radius,
        final_trust_region_radius,
        max_function_calls_count,
        working_space.data()
    );
    CHECK(result == 0.49700493360635966677563146731699816882610321044921875);
    CHECK(variables_values == Values({{
        0.498999933471093226611259296987554989755153656005859375,
        0.498000000002213061289779716389602981507778167724609375
    }}));
}

template <class F>
BobyqaClosure make_closure(F &function) {
    struct Wrap {
        static double call(void *data, long n, const double *values) {
            return reinterpret_cast<F *>(data)->operator()(n, values);
        }
    };
    return BobyqaClosure {&function, &Wrap::call};
}

TEST_CASE("Uniform random function with funciton calls calculation", "[bobyqa_closure]") {
    class Function {
    public:
        Function() {
            generator.seed(0);
        }

        double operator ()(long n, const double *) {
            REQUIRE(n == 2);
            return distribution(generator);
        }

    private:
        std::random_device device;
        std::mt19937 generator = std::mt19937(device());
        std::uniform_real_distribution<> distribution = std::uniform_real_distribution<>(-1, 1);
    } function;
    auto closure = make_closure(function);
    const long variables_count = 2;
    const long number_of_interpolation_conditions = variables_count + 2;
    using Values = std::array<double, variables_count>;
    Values variables_values = {{0, 0}};
    const Values lower_bound = {{-1, -1}};
    const Values upper_bound = {{1, 1}};
    const double initial_trust_region_radius = 1e-3;
    const double final_trust_region_radius = 1e3;
    const long max_function_calls_count = 5;
    const std::size_t working_space_size = BOBYQA_WORKING_SPACE_SIZE(
        variables_count, number_of_interpolation_conditions);
    using WorkingSpace = std::array<double, working_space_size>;
    WorkingSpace working_space;
    const double result = bobyqa_closure(
        &closure,
        variables_count,
        number_of_interpolation_conditions,
        variables_values.data(),
        lower_bound.data(),
        upper_bound.data(),
        initial_trust_region_radius,
        final_trust_region_radius,
        max_function_calls_count,
        working_space.data()
    );
    CHECK(result == 0.1856892330333652640916852760710753500461578369140625);
    CHECK(variables_values == Values({{0, 0}}));
}
