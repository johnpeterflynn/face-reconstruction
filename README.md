# Face Reconstruction from Kinect Scan

This is a reproduction of part of the work of Dr. Justus Thies in his paper _Real-time Expression Transfer for Facial Reenactment_ [1] to reconstruct high resolution facial expression and geometry form a Kinect face scan.

![](/content/example.png)

## Dependencies
face-reconstruction uses [CMake](https://cmake.org/) as a build system with the following dependencies:

* [ceres-solver](https://github.com/ceres-solver/ceres-solver)
* [eigen](https://gitlab.com/libeigen/eigen)
* [Sophus](https://github.com/strasdat/Sophus)
* [nabo](https://github.com/ethz-asl/libnabo)

In addition, a face scan and face landmarks .ply files are expected from Kinect HD Face in the scan subfolder (see config.h).

## Notes
This repository does not include the basis for face geometry and expression.


## References
1. Justus Thies, Michael Zollhöfer, Matthias Nießner, Levi Valgaerts, Marc Stamminger,
and Christian Theobalt. 2015. Real-time Expression Transfer for Facial Reenactment.
ACM TOG 34, 6, Article 183 (2015)
