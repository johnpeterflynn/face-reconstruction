#include "rgbdscan.h"

#include <Eigen/Dense>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/getopt.h>
#include "nabo/nabo.h"

#include "config.h"

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;


RGBDScan::RGBDScan(const std::string& face_path,
                   const std::string& landmark_path) {
    loadMesh(face_path);

    // Must be called after loadMesh() because it uses m_scanned_mesh
    loadMatchIndices(landmark_path);
}

void RGBDScan::loadMatchIndices(const std::string& landmark_path) {

    OpenMesh::IO::Options ropt;

    if (OpenMesh::IO::read_mesh(m_landmarks, landmark_path, ropt)) {
        // Load face vertices
        Eigen::MatrixXf face_vertices(3, m_scanned_mesh.n_vertices());
        for (MyMesh::VertexIter v_it = m_scanned_mesh.vertices_begin(); v_it != m_scanned_mesh.vertices_end(); ++v_it)
        {
          MyMesh::Point p3 = m_scanned_mesh.point(*v_it);
          Eigen::Vector3f v3(p3[0], p3[1], p3[2]);

          face_vertices.col(v_it->idx()) = v3;
        }

        // Load Landmark vertices
        Eigen::MatrixXf landmark_vertices(3, m_landmarks.n_vertices());
        for (MyMesh::VertexIter v_it = m_landmarks.vertices_begin(); v_it != m_landmarks.vertices_end(); ++v_it)
        {
          MyMesh::Point p3 = m_landmarks.point(*v_it);
          Eigen::Vector3f v3(p3[0], p3[1], p3[2]);
          // std::cout << p3[0] << "\n";

          landmark_vertices.col(v_it->idx()) = v3;
        }

        // Get Closest Neighbor
        const int K = 1;

        // Kinect HD Face fits a landmark model to a face scan. A subset of those
        // landmarks are manually selected that correspond strongly with pixels
        // on the average mesh model (for example, landmarks at the edges of the
        // lips, eyes, chin and top of forehead. Since the Kinect HD face and
        // average mesh model are always the same, this correspondence only
        // needs to be selected once.

        // Add relevant model indices
        std::vector<int> model_indices;
        model_indices.push_back(8288);
        model_indices.push_back(8320);
        model_indices.push_back(11039);
        model_indices.push_back(12273);
        model_indices.push_back(10214);
        model_indices.push_back(14346);
        model_indices.push_back(4275);
        model_indices.push_back(5620);
        model_indices.push_back(2089);
        model_indices.push_back(4287);
        model_indices.push_back(12156);

        int num_correspondences = model_indices.size();

        // Add relevant landmark indices
        Eigen::MatrixXf q(3, num_correspondences);
        q.col(0) = landmark_vertices.col(11);
        q.col(1) = landmark_vertices.col(21);
        q.col(2) = landmark_vertices.col(145);
        q.col(3) = landmark_vertices.col(244);
        q.col(4) = landmark_vertices.col(329);
        q.col(5) = landmark_vertices.col(470);
        q.col(6) = landmark_vertices.col(731);
        q.col(7) = landmark_vertices.col(763);
        q.col(8) = landmark_vertices.col(1010);
        q.col(9) = landmark_vertices.col(1090);
        q.col(10) = landmark_vertices.col(1105);

        Nabo::NNSearchF* nns = Nabo::NNSearchF::createKDTreeLinearHeap(face_vertices);
        Eigen::MatrixXi indices(K,num_correspondences);
        Eigen::MatrixXf dists2(K,num_correspondences);

        // ALLOW_SELF_MATCH appears to be necessary to match vertices that are
        // exactly the same.
        nns->knn(q, indices, dists2, K, 0, Nabo::NNSearchF::ALLOW_SELF_MATCH);
        //std::cout << indices << "\n";
        delete nns;

        for (int i = 0; i < num_correspondences; i++) {
            m_match_indices[model_indices[i]] = indices(0, i);
        }
    }
    else {
        // NOTE: This defaults to hard-coded demo mode for Justus' face model.

        std::cerr << "ERROR: Could not load " << landmark_path << "\n";

        // WARNING: Indices for the Justus' face model. Not generalized
        m_match_indices[12280] = 9734;
        m_match_indices[4280] = 10045;
        m_match_indices[8319] = 12539;
        m_match_indices[48246] = 19679;
        return;
    }
}

std::map<int, int> RGBDScan::getMatchIndices() {
    return m_match_indices;
}

int RGBDScan::loadMesh(const std::string& path) {
    OpenMesh::IO::Options ropt;

    // Set input options
    ropt += OpenMesh::IO::Options::VertexColor;

    m_scanned_mesh.request_vertex_colors();

    // assure we have vertex normals
    if (!m_scanned_mesh.has_vertex_colors())
    {
      std::cerr << "ERROR: Standard vertex property 'Colors' not available for scanned mesh!\n";
      //return 1;
    }

    if (OpenMesh::IO::read_mesh(m_scanned_mesh, path, ropt)) {
        std::cout << "Loaded " << path << "\n";
    }
    else if (OpenMesh::IO::read_mesh(m_scanned_mesh, PATH_DEFAULT_SCANNED_MESG, ropt)) {
        std::cout << "Loaded " << PATH_DEFAULT_SCANNED_MESG << "\n";
    }

    return 0;
}