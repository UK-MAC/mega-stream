# MEGA-STREAM

The mega-stream benchmark aims to test the theory that streaming many arrays (with different sizes)
causes memory bandwidth limits not to be reached; resulting in latency becoming a dominant factor.

We run a kernel with a similar form to STREAM Triad, but with more than 3 input arrays.

## The kernel
The main kernel consists of 8 arrays with the following properties:

 * r and q are large
 * x, y and z are medium
 * a, c and c are small

```r[i] = q[i] + a[i%S]*x[i%M] + b[i%S]*y[i%M] + c[i%S]*z[i%M]```

## Building
The benchmark should build with `make`, and by default uses the Intel compiler.
This can be changed by specifying `CC`, for example `make CC=cc`.
Additional options can be passed to the Makefile as `make OPTIONS=`.

