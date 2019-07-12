#ifndef FACEMODEL_H
#define FACEMODEL_H

#include <Eigen/Dense>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

static constexpr float SCALE_AVG_MESH = 1.0 / 1000000.0;

typedef struct float4 {
        float x;
        float y;
        float z;
        float w;
} float4;

class FaceModel
{
public:
    typedef OpenMesh::TriMesh_ArrayKernelT<> FaceMesh;

    FaceModel(const std::string& path, int n_eigenvec, int n_exp, int n_vert);

    void forwardPass(const Eigen::VectorXd& alpha, const Eigen::VectorXd& delta,
                      Eigen::VectorXf& vertices_out);

    int synthesizeModel(const Eigen::VectorXf& diff_vertices);
    int writeSynthesizedModel();

private:
    void load(const std::string& path);
    int loadAverageMesh();

    void loadVector(const std::string &filename, float *res, unsigned int length);
    void progressBar(const char* str, int num);
    float* loadEigenvectors(const std::string &filename, unsigned int components, unsigned int numberOfEigenvectors);
    int convertBasis(Eigen::MatrixXf& basis, int nVertices, int numCols, float4* basisCPU);
    int convertDeviations(Eigen::VectorXf& devs, int num_dims, float* devCPU);

public:
    FaceMesh m_avg_mesh;
    FaceMesh m_synth_mesh;

    const int NumberOfEigenvectors;
    const int NumberOfExpressions;
    const int nVertices;

    Eigen::MatrixXf shapeBasisEigen;
    Eigen::MatrixXf exprBasisEigen;

    Eigen::VectorXf shapeDevEigen;
    Eigen::VectorXf exprDevEigen;

};

#endif // FACEMODEL_H