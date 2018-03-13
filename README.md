# mega-stream and mega-sweep

## mega-stream

The mega-stream mini-app initially aimed to test the theory that streaming many arrays (with different sizes)
causes memory bandwidth limits not to be reached; resulting in latency becoming a dominant factor.
We ran a kernel with a similar form to STREAM Triad, but with more than 3 input arrays.
Additionally we run a small reduction, requiring results of the Triad-style computation.

This was then extended to also include a finite difference style update on the medium side arays.

### The kernel

The main kernel consists of 8 multi-dimensional arrays with the following properties:

 * r and q are large
 * x, y and z are medium
 * a, b, and c are small

The computational kernel is found inside a triple-nested loop, and can be expressed as

```
r(i,j,k,l,m) = q(i,j,k,l,m) + a(i)*x(i,j,k,m) + b(i)*y(i,j,l,m) + c(i)*z(i,k,l,m)
x(i,j,k,m) = 0.2*r(i,j,k,l,m) - x(i,j,k,m)
y(i,j,l,m) = 0.2*r(i,j,k,l,m) - y(i,j,l,m)
z(i,k,l,m) = 0.2*r(i,j,k,l,m) - z(i,k,l,m)
sum(j,k,l,m) = SUM(r(:,j,k,l,m))
```

### Building

The benchmark should build with `make`, and by default uses the Intel compiler.
This can be changed by specifying `CC`, for example `make CC=cc`.
Additional options can be passed to the Makefile as `make OPTIONS=`.

### Notes

The Fortran version does not have any command line argument checking.
A baseline and an optimised version are kept in this repository.


## mega-sweep

The mega-sweep mini-app aims to take the mega-stream kernel and fit it within a KBA-style sweep iteration structure.
The compute kernel is similar to mega-stream, but additional looping along with MPI communication is included.

This mini-app is Fortran only. There are 2- and 3-spatial dimensions versions available.

### Building

The benchmark should build with `make`, and by default uses the default MPI compiler.
A compiler toolchain should be specified with `COMPILER`, for example `make COMPILER=INTEL`.
Additional options can be passed to the Makefile as `make OPTIONS=`.


## Publications

- Tom Deakin, John Pennycook, Andrew Mallinson, Wayne Gaudin and Simon McIntosh-Smith. The MEGA-STREAM Benchmark on Intel Xeon Phi Processors (Knights Landing). The Intel Xeon Phi Users Group Spring Meeting, 2017.
- Tom Deakin, Simon McIntosh-Smith and Wayne Gaudin. On the Mitigation of Cache Hostile Memory Access Patterns on Many-core CPU Architectures. The Intel Xeon Phi Users Group Workshop at International Conference on High Performance Computing, 2017.

