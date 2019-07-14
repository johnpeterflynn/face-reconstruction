#ifndef FACESOLVER_H
#define FACESOLVER_H

#include <iostream>
#include <set>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include "facemodel.h"
#include "rgbdscan.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

class FaceSolver {
public:
    FaceSolver();

    void solve(FaceModel& face_model, RGBDScan face_scan, Eigen::VectorXd& alpha, Eigen::VectorXd& delta);

private:
    void calculate_knn(const Eigen::MatrixXf& M, const Eigen::MatrixXf& q,
                       Eigen::MatrixXi& indices);
    int knn_model_to_scan(const FaceModel& face_model, const MyMesh& scanned_mesh, int K, Eigen::MatrixXi& indices);
    void meshToMatrix(const MyMesh& mesh, Eigen::MatrixXf& M_out) ;

    void runCeres(const MyMesh& avg_face_mesh, const MyMesh& scanned_mesh,
                  const Eigen::MatrixXi& indices,
                  const std::set<int>& matches,
                  const Eigen::MatrixXf& shapeBasisEigen,
                  const Eigen::MatrixXf& exprBasisEigen,
                  const Eigen::VectorXf& shapeDevEigen,
                  const Eigen::VectorXf& exprDevEigen,
                  Eigen::VectorXd& alpha, Eigen::VectorXd& delta);

private:
};

#endif // FACESOLVER_H
