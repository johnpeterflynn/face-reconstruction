#ifndef FACESOLVER_H
#define FACESOLVER_H

#include <iostream>
#include <set>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <sophus/se3.hpp>

#include "facemodel.h"
#include "rgbdscan.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

class FaceSolver {
public:
    FaceSolver(double geo_regularization, int num_iterations);

    void solve(FaceModel& face_model, RGBDScan face_scan, Eigen::VectorXd& alpha, Eigen::VectorXd& delta, Sophus::SE3d& T_xy);

private:
    void calculate_knn(const Eigen::MatrixXf& M, const Eigen::MatrixXf& q,
                       Eigen::MatrixXi& indices);
    int knn_model_to_scan(const MyMesh& synth_mesh, const MyMesh& scanned_mesh, int K, Eigen::MatrixXi& indices);
    void meshToMatrix(const MyMesh& mesh, Eigen::MatrixXf& M_out) ;

    void runCeres(const MyMesh& avg_face_mesh, const MyMesh& scanned_mesh,
                  const Eigen::MatrixXi& indices,
                  const std::map<int, int>& match_indices,
                  const Eigen::MatrixXf& shapeBasisEigen,
                  const Eigen::MatrixXf& exprBasisEigen,
                  const Eigen::VectorXf& shapeDevEigen,
                  const Eigen::VectorXf& exprDevEigen,
                  bool weigh_separately, int max_num_iterations,
                  Eigen::VectorXd& alpha, Eigen::VectorXd& delta,
                  Sophus::SE3d& T_xy);

private:
    double m_geo_regularization;
    int m_num_iterations;
};

#endif // FACESOLVER_H
