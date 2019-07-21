#include "pch.h"
#include <iostream>

#include <string>
#include <Eigen/Dense>
#include <vector>

#include <OpenMesh/Tools/Utils/getopt.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

int main()
{
    MyMesh face_mesh;
    MyMesh landmark_mesh;
    OpenMesh::IO::Options ropt;
    std::string face_path = "/Users/ardakeskiner/Downloads/2019_07_19_16_37_07/FaceRaw/Face_Raw_000000.ply";
    std::string landmark_path = "/Users/ardakeskiner/Downloads/2019_07_19_16_37_07/HDFace/HD_Face_000000.ply";

    // Set input options
    ropt += OpenMesh::IO::Options::VertexColor;

    if (!OpenMesh::IO::read_mesh(face_mesh, face_path, ropt)) {
        std::cerr << "ERROR: Could not load " << face_path << "\n";
        return 1;
    }

    if (!OpenMesh::IO::read_mesh(landmark_mesh, landmark_path, ropt)) {
        std::cerr << "ERROR: Could not load " << landmark_path << "\n";
        return 1;
    }

    // Load face vertices
    Eigen::MatrixXf face_vertices(3, face_mesh.n_vertices());
    for (MyMesh::VertexIter v_it = face_mesh.vertices_begin(); v_it != face_mesh.vertices_end(); ++v_it)
    {
      MyMesh::Point p3 = face_mesh.point(*v_it);
      Eigen::Vector3f v3(p3[0], p3[1], p3[2]);

      face_vertices.col(v_it->idx()) = v3;
    }

    // Load Landmark vertices
    Eigen::MatrixXf landmark_vertices(3, landmark_mesh.n_vertices());
    for (MyMesh::VertexIter v_it = landmark_mesh.vertices_begin(); v_it != landmark_mesh.vertices_end(); ++v_it)
    {
      MyMesh::Point p3 = landmark_mesh.point(*v_it);
      Eigen::Vector3f v3(p3[0], p3[1], p3[2]);
      // std::cout << p3[0] << "\n";

      landmark_vertices.col(v_it->idx()) = v3;
    }

    // Get Closest Neighbor
    const int K = 1;

    // Get Facial 4 Landmarks
    Eigen::MatrixXf q(3, 5);
    q.col(0) = landmark_vertices.col(469); // HighDetailFacePoints_LefteyeOutercorner // 14728 in average face model 44355.4 33340 85867.5 // 470
    q.col(1) = landmark_vertices.col(1117); // HighDetailFacePoints_RighteyeOutercorner // 1700 in average face model -45522.1 33283 85384 // 1118
    q.col(2) = landmark_vertices.col(18); // HighDetailFacePoints_NoseTip // 8320 in average face model 20.1202 2232.26 131944 // 17
    q.col(3) = landmark_vertices.col(4); // HighDetailFacePoints_ChinCenter // 48112 in average face model -1978.04 -58610.7 108977 // 3
    q.col(4) = landmark_vertices.col(28); // HighDetailFacePoints_ForeheadCenter // 40824 in average face model 973.729 66819.7 108082 // 
    

    Nabo::NNSearchF* nns = Nabo::NNSearchF::createKDTreeLinearHeap(face_vertices);
    Eigen::MatrixXi indices(K,5);
    Eigen::MatrixXf dists2(K,5);

    // ALLOW_SELF_MATCH appears to be necessary to match vertices that are
    // exactly the same.
    nns->knn(q, indices, dists2, K, 0, Nabo::NNSearchF::ALLOW_SELF_MATCH);
    std::cout << indices << "\n";
    delete nns;
}