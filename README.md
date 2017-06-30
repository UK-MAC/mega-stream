# MEGA-STREAM

The mega-stream benchmark aims to test the theory that streaming many arrays (with different sizes)
causes memory bandwidth limits not to be reached; resulting in latency becoming a dominant factor.

We run a kernel with a similar form to STREAM Triad, but with more than 3 input arrays.
Additionally we run a small reduction, requiring results of the Triad-style computation.

## The kernel
The main kernel consists of 8 multi-dimensional arrays with the following properties:

 * r and q are large
 * x, y and z are medium
 * a, b, and c are small

The computational kernel is found inside a triple-nested loop, and can be expressed as

```
r(i,j,k,l,m) = q(i,j,k,l,m) + a(i)*x(i,j,k,m) + b(i)*y(i,j,l,m) + c(i)*z(i,k,l,m)
sum(j,k,l,m) = SUM(r(:,j,k,l,m))
```

## Building
The benchmark should build with `make`, and by default uses the Intel compiler.
This can be changed by specifying `CC`, for example `make CC=cc`.
Additional options can be passed to the Makefile as `make OPTIONS=`.

## Notes
The Fortran version does not have any command line argument checking.

