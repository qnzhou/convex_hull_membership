#pragma once

#include <predicates.h>

#include <cassert>
#include <cmath>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace convex_hull_membership {

namespace details {

template <typename T>
bool is_separating_plane(std::span<T> pts, std::span<T> query_point, size_t i,
                         size_t j, size_t k) {
    constexpr size_t DIM = 3;
    assert(query_point.size() == DIM);
    T* pa = pts.data() + DIM * i;
    T* pb = pts.data() + DIM * j;
    T* pc = pts.data() + DIM * k;
    T* pd = query_point.data();
    const auto ori_query = orient3d(pa, pb, pc, pd);
    if (ori_query == 0) {
        // Query point is coplanar with the triangle formed by i, j, k.
        // It is contained in the convex hull by definition.
        return false;
    }

    const size_t num_pts = pts.size() / DIM;
    for (size_t l = 0; l < num_pts; l++) {
        if (l == i || l == j || l == k) {
            continue;
        }
        T* pe = pts.data() + DIM * l;
        const auto ori = orient3d(pa, pb, pc, pe);
        if (ori * ori_query > 0) {
            return false;
        }
    }
    return true;
}

template <typename T>
bool is_separating_plane(std::span<T> pts, std::span<T> query_point, size_t i,
                         size_t j) {
    constexpr size_t DIM = 2;
    assert(query_point.size() == DIM);
    T* pa = pts.data() + DIM * i;
    T* pb = pts.data() + DIM * j;
    T* pc = query_point.data();
    const auto ori_query = orient2d(pa, pb, pc);
    if (ori_query == 0) {
        // Query point is coplanar with the triangle formed by i, j, k.
        // It is contained in the convex hull by definition.
        return false;
    }

    const size_t num_pts = pts.size() / DIM;
    for (size_t l = 0; l < num_pts; l++) {
        if (l == i || l == j) {
            continue;
        }
        T* pd = pts.data() + DIM * l;
        const auto ori = orient2d(pa, pb, pd);
        if (ori * ori_query > 0) {
            return false;
        }
    }
    return true;
}

inline int det2_filtered(double p0, double p1, double q0, double q1) {
    double d1 = p0 * q1;
    double d2 = p1 * q0;
    double m = d1 - d2;

    double _tmp_fabs;

    double max_var = 0.0;
    if ((_tmp_fabs = std::fabs(p0)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = std::fabs(p1)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = std::fabs(q0)) > max_var) max_var = _tmp_fabs;
    if ((_tmp_fabs = std::fabs(q1)) > max_var) max_var = _tmp_fabs;
    double epsilon = max_var;
    epsilon *= epsilon;
    epsilon *= 4.440892098500627e-16;
    if (m > epsilon) return 1;
    if (-m > epsilon) return -1;
    return 0;
}

inline int det2(std::span<double> pa, std::span<double> pb) {
    auto d = pa[0] * pb[1] - pa[1] * pb[0];
    if (d > 0) return 1;
    if (d < 0) return -1;
    return 0;
}

}  // namespace details

template <int DIM, typename T>
bool contains(std::span<T> pts, std::span<T> query_point) {
    static_assert(std::is_same_v<T, IGL_PREDICATES_REAL>);
    exactinit();
    const size_t num_pts = pts.size() / DIM;

    if constexpr (DIM == 3) {
        for (size_t i = 0; i < num_pts; i++) {
            for (size_t j = i + 1; j < num_pts; j++) {
                for (size_t k = j + 1; k < num_pts; k++) {
                    if (details::is_separating_plane(pts, query_point, i, j, k))
                        return false;
                }
            }
        }
    } else if constexpr (DIM == 2) {
        for (size_t i = 0; i < num_pts; i++) {
            for (size_t j = i + 1; j < num_pts; j++) {
                if (details::is_separating_plane(pts, query_point, i, j))
                    return false;
            }
        }
    }
    return true;
}

template <int DIM, int N, typename T>
bool contains_origin(std::span<T, DIM * N> pts) {
    static_assert(DIM == 2);
    static_assert(std::is_same_v<T, IGL_PREDICATES_REAL>);

    std::array<T, N * (N - 1) / 2> queries;
    size_t count = 0;
    for (size_t i = 0; i < N; i++) {
        std::span<T, 2> pi{pts.data() + DIM * i, 2};
        for (size_t j = i + 1; j < N; j++) {
            std::span<T, 2> pj{pts.data() + DIM * j, 2};
            const size_t idx = (N - 1 + N - i) * i / 2 + j - i - 1;
            assert(idx == count);
            queries[idx] = details::det2(pi, pj);
            count++;
        }
    }
    assert(count == N * (N-1) / 2);

    for (size_t i = 0; i < N; i++) {
        for (size_t j = i + 1; j < N; j++) {
            const size_t idx_ij = (N - 1 + N - i) * i / 2 + j - i - 1;
            for (size_t k = j + 1; k < N; k++) {
                const size_t idx_ik = (N - 1 + N - i) * i / 2 + k - i - 1;
                const size_t idx_jk = (N - 1 + N - j) * j / 2 + k - j - 1;
                if (queries[idx_ij] <= 0 && queries[idx_jk] <= 0 && queries[idx_ik] >= 0) return true;
                if (queries[idx_ij] >= 0 && queries[idx_jk] >= 0 && queries[idx_ik] <= 0) return true;
            }
        }
    }
    return false;
}

}  // namespace convex_hull_membership
