#ifndef FACEMODEL_H
#define FACEMODEL_H

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

static constexpr float SCALE_AVG_MESH = 1 / 1000000.0;

class FaceModel
{
public:
    FaceModel();

    int loadAverageMesh();

private:
    typedef OpenMesh::TriMesh_ArrayKernelT<>  AvgMesh;

private:
    AvgMesh m_avg_mesh;

};

#endif // FACEMODEL_H