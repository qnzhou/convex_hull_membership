#include <convex_hull_membership/contains.h>

#include <array>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <random>

template <typename T, size_t DIM>
void check() {
    constexpr size_t N = 20;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-10, 10);

    std::array<T, N * DIM> pts;
    std::for_each(pts.begin(), pts.end(), [&](auto& pt) { pt = dist(gen); });

    std::array<T, DIM> query;
    std::for_each(query.begin(), query.end(), [&](auto& pt) { pt = 0; });
    REQUIRE(convex_hull_membership::contains<DIM, T>(pts, query));

    std::for_each(query.begin(), query.end(), [&](auto& pt) { pt = 20; });
    REQUIRE(!convex_hull_membership::contains<DIM, T>(pts, query));

    for (size_t i = 0; i < N; ++i) {
        std::span<T> p(pts.data() + i * DIM, DIM);
        REQUIRE(convex_hull_membership::contains<DIM, T>(pts, p));
    }
}

template <typename T>
void check_origin() {
    constexpr size_t N = 35;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(1, 10);

    std::array<T, N * 2> pts;
    std::for_each(pts.begin(), pts.end(), [&](auto& pt) { pt = dist(gen); });

    REQUIRE(!convex_hull_membership::contains_origin<2, N, T>(pts));

    pts[0] = 0;
    pts[1] = 0;
    REQUIRE(convex_hull_membership::contains_origin<2, N, T>(pts));

    pts[0] = -pts[2] - pts[4];
    pts[1] = -pts[3] - pts[5];
    REQUIRE(convex_hull_membership::contains_origin<2, N, T>(pts));
}

TEST_CASE("convex_hull_membership", "[check]") {
    using Scalar = double;
    check<Scalar, 2>();
    check<Scalar, 3>();

    check_origin<Scalar>();
}

TEST_CASE("benchmark", "[convext_hull_membership][.benchmark]") {
    using Scalar = double;
    constexpr size_t N = 35;
    BENCHMARK_ADVANCED("2D")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t DIM = 2;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(-10, 10);

        std::array<Scalar, N * DIM> data;
        std::for_each(data.begin(), data.end(),
                      [&](auto& pt) { pt = dist(gen); });

        std::array<Scalar, DIM> query;
        std::for_each(query.begin(), query.end(),
                      [&](auto& pt) { pt = dist(gen); });

        meter.measure([&] {
            return convex_hull_membership::contains<DIM, Scalar>(data, query);
        });
    };

    BENCHMARK_ADVANCED("2D origin")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t DIM = 2;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0, 20);

        std::array<Scalar, N * DIM> data;
        std::for_each(data.begin(), data.end(),
                      [&](auto& pt) { pt = dist(gen); });

        meter.measure([&] {
            return convex_hull_membership::contains_origin<DIM, N, Scalar>(data);
        });
    };

    BENCHMARK_ADVANCED("3D")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t DIM = 3;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(-10, 10);

        std::array<Scalar, N * DIM> data;
        std::for_each(data.begin(), data.end(),
                      [&](auto& pt) { pt = dist(gen); });

        std::array<Scalar, DIM> query;
        std::for_each(query.begin(), query.end(),
                      [&](auto& pt) { pt = dist(gen); });

        meter.measure([&] {
            return convex_hull_membership::contains<DIM, Scalar>(data, query);
        });
    };
}
