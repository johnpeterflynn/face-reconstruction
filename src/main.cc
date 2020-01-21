#include "pch.h"
#include <iostream>

#include <string>
#include <Eigen/Dense>
#include <vector>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

#include <chrono>

#include "facemodel.h"
#include "rgbdscan.h"
#include "facesolver.h"
#include"config.h"


int main()
{
    Eigen::VectorXd alpha = Eigen::VectorXd::Zero(NumberOfEigenvectors);
    Eigen::VectorXd delta = Eigen::VectorXd::Zero(NumberOfExpressions);

    Sophus::SE3d T_xy;

    FaceModel face_model(MODEL_PATH);
    RGBDScan face_scan(FILENAME_SCANNED_MESH, FILENAME_SCANNED_LANDMARKS);
    FaceSolver face_solver(0.0001, 0.004, 0.003, 8, 0.8, true);// percent = 0.1 works pretty okay too

    auto start = std::chrono::high_resolution_clock::now();
    face_solver.solve(face_model, face_scan, alpha, delta, T_xy);
    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << "Writing synthesized model to file\n";
    face_model.writeSynthesizedModel(alpha, delta, T_xy);

    // Print alpha (geometry parameter we solved for)
    //std::cout << alpha << std::endl;

    std::cout << "Transform" << std::endl;
    std::cout << T_xy.matrix() << std::endl;

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "Face solver took: "
             << duration.count() / 1000000.0 << " seconds" << std::endl;

}
