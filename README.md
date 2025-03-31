gf2poly
===

This rust crate provides a `Gf2Poly` type which implements polynomial arithmetic over GF(2).
It links against the `gf2x` C library and care is taken to be asymptotically efficient.
For example, multiplication is `n log n` (because the `gf2x` implementation is), and as a result, division can also be implemented in `n log n`.
This crate also implements a fast `gcd` in `n log² n` time and a basic implementation of factorization.

Building
--------
This crate requires the `gf2x` library to be installed.
It can either be installed through a package manager if you're lucky (but note that this is going to be a slow version because distro maintainers are generally conservative in what CPU features they enable), or it can be compiled using a release [from INRIA's gitlab](https://gitlab.inria.fr/gf2x/gf2x).

`gf2poly` provides two environment variables for controlling the linking:
 * `GF2POLY_STATIC_LIB`: If this is `1`, `gf2x` will be linked statically.
 * `GF2POLY_LIBRARY_PATH`: An additional path for the linker to search for static libraries at compile time.
