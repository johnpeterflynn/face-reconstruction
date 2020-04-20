# Face Reconstruction from Kinect Scan

<p align="center">
  <img height="400" src="/content/basic-anim.gif">
</p>

This is a reproduction of part of the work of Dr. Justus Thies in his paper _Real-time Expression Transfer for Facial Reenactment_ [1] to reconstruct high resolution facial expression and geometry (gray) from a Kinect face scan (red).

### Algorithm

See `src/facesolver.cpp` for main optimization loop.

1. Minimize the sum-squared-distance error for sparse correspondence points between face model and scan.
2. For number of iterations `num-iters`:
   1. Randomly sample a subset `frac-used-vertices` of model vertices.
   2. Find new correspondence points by computing the 1-KNN of each model vertex to scan vertices.
   3. Filter scan vertices at border (if `ignore-borders = true`) and model vertices not in `model/mesh/averageMesh_blackface.off`.
   4. Minimize the sum-squared-distance error for all correspondence points that are below a distance threshold `knn-dist-thresh`.

### Optimization

Significant performance improvements were made possible by only optimizing over the first 40 PCA components for geometry and 30 PCA components for expression (see `include/config.h`).

We achieved further performance improvements by randomly sampling a subset of face scan vertices during optimization (`frac-used-vertices`). Below are the results of using various sampling percentages and their runtimes on a commodity Intel-i5 CPU. By qualitative inspection, even using 20% of available vertices achieves nearly the same reconstruction as when using 100% of available vertices in just a third of the time.

![](/content/frac-vertices-time.png)



## Build and Run

To build, clone the repository and run the following in the root directory.

```
mkdir build & cd build
cmake ..
make
```



To run a reconstruction you must specify a scan file, correspondence file and output file. Example run:

```
face_recon_2019 --scan scan/kd_example.off --corr scan/sparse.corr --out synth.off
```



There are additional options to adjust optimization parameters including fraction of vertices to sample and number of iterations. For a full list of options with their descriptions, run:

```
face_recon_2019 --help
```



## Dependencies

face-reconstruction uses [CMake](https://cmake.org/) as a build system with the following dependencies:

* [ceres-solver](https://github.com/ceres-solver/ceres-solver)
* [eigen](https://gitlab.com/libeigen/eigen)
* [Sophus](https://github.com/strasdat/Sophus)
* [nabo](https://github.com/ethz-asl/libnabo)
* [Boost](https://github.com/boostorg/boost)

In addition, face scan and face landmark .ply files are expected from Kinect HD Face in the scan subfolder (see config.h).

## Notes
This repository does not include the basis for face geometry and expression.


## References
1. Justus Thies, Michael Zollhöfer, Matthias Nießner, Levi Valgaerts, Marc Stamminger,
and Christian Theobalt. 2015. Real-time Expression Transfer for Facial Reenactment.
ACM TOG 34, 6, Article 183 (2015)
