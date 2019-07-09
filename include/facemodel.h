#ifndef FACEMODEL_H
#define FACEMODEL_H

#include <Eigen/Dense>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

static constexpr float SCALE_AVG_MESH = 1.0 / 1000000.0;

class FaceModel
{
public:
    typedef OpenMesh::TriMesh_ArrayKernelT<> FaceMesh;

    FaceModel();

    int loadAverageMesh();
    int writeSynthesizedModel(const Eigen::VectorXf& diff_vertices);

public:
    FaceMesh m_avg_mesh;

};

#endif // FACEMODEL_H