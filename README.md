# Brute-force Convex Hull Membership Test

The repo is a header-only library implementing a brute-force convex hull membership test.

## 2D Test

```c++
#include <convex_hull_membership/contains.h>

std::vector<double> pts; // X0, Y0, X1, Y1, ...
std::vector<double> query; // X, Y

bool r = convex_hull_membership::contains<2, double>(pts, query);
```

## 3D Test

```c++
#include <convex_hull_membership/contains.h>

std::vector<double> pts; // X0, Y0, Z0, X1, Y1, Z0, ...
std::vector<double> query; // X, Y, Z

bool r = convex_hull_membership::contains<3, double>(pts, query);
```

## Benchmark

Tested on Macbook Pro with 2.4 GHz 8-Core Intel Core i9 and 64G memory.

```
benchmark name                       samples       iterations    estimated
                                     mean          low mean      high mean
                                     std dev       low std dev   high std dev
-------------------------------------------------------------------------------
2D                                             100            17     4.5883 ms
                                          2.169 us    1.97718 us     2.3594 us
                                        978.984 ns    861.917 ns    1.18144 us

3D                                             100             1     4.3963 ms
                                        29.4295 us    25.1695 us    34.0437 us
                                         22.663 us    20.4476 us    26.5526 us
```
