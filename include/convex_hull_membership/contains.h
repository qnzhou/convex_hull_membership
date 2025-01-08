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
    int r = det2_filtered(pa[0], pa[1], pb[0], pb[1]);
    if (r != 0) {
        return r;
    } else {
        double origin[2] = {0, 0};
        return orient2d(pa.data(), pb.data(), origin);
    }
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

template <int DIM, typename T>
bool contains_origin(std::span<T> pts) {
    static_assert(std::is_same_v<T, IGL_PREDICATES_REAL>);
    exactinit();
    const size_t num_pts = pts.size() / DIM;

    if constexpr (DIM == 2) {
        for (size_t i = 0; i < num_pts; i++) {
            T* pi = pts.data() + DIM * i;
            for (size_t j = i + 1; j < num_pts; j++) {
                T* pj = pts.data() + DIM * j;

                const auto ori_ij = details::det2({pi, 2}, {pj, 2});
                if (ori_ij == 0) return true;

                for (size_t k = j + 1; k < num_pts; k++) {
                    T* pk = pts.data() + DIM * k;

                    const auto ori_jk = details::det2({pj, 2}, {pk, 2});
                    const auto ori_ki = details::det2({pk, 2}, {pi, 2});
                    if (ori_ij <= 0 && ori_jk <= 0 && ori_ki <= 0) return true;
                    if (ori_ij >= 0 && ori_jk >= 0 && ori_ki >= 0) return true;
                }
            }
        }
    } else {
        throw std::runtime_error("Not implemented");
    }
    return false;
}

}  // namespace convex_hull_membership
