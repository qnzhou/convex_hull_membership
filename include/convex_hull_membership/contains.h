#pragma once

#include <predicates.h>

#include <cassert>
#include <span>
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

}  // namespace convexhull_membership
