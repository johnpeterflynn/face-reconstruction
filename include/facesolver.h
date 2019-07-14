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
    FaceSolver(int n_eigenvecs, int n_expr, int n_vertices);

    void solve(FaceModel& face_model, RGBDScan face_scan, Eigen::VectorXd& alpha, Eigen::VectorXd& delta);

private:
    void calculate_knn(const Eigen::MatrixXf& M, const Eigen::MatrixXf& q,
                       Eigen::MatrixXi& indices);
    int knn_test(const FaceModel& face_model, const MyMesh& scanned_mesh, int K, Eigen::MatrixXi& indices);
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
    //const int NumberOfEigenvectors;
    //const int NumberOfExpressions;
    //const int nVertices;
};

#endif // FACESOLVER_H
