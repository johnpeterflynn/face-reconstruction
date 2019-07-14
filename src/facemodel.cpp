#include "facemodel.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/getopt.h>


constexpr const char* FILENAME_AVG_MESH = "../models/averageMesh.off";
constexpr const char* FILENAME_OUT_SYNTH_MESH = "synthesizedMesh.off";

constexpr const char* filenameBasisShape = "ShapeBasis.matrix";
constexpr const char* filenameBasisExpression = "ExpressionBasis.matrix";
constexpr const char* filenameStdDevShape = "StandardDeviationShape.vec";
constexpr const char* filenameStdDevExpression = "StandardDeviationExpression.vec";


FaceModel::FaceModel(const std::string &path, int n_eigenvec, int n_exp, int n_vert)
    : NumberOfEigenvectors(n_eigenvec), NumberOfExpressions(n_exp),
      nVertices(n_vert),
    shapeBasisEigen(3 * nVertices, NumberOfEigenvectors),
    exprBasisEigen(3 * nVertices, NumberOfExpressions),
    shapeDevEigen(NumberOfEigenvectors),
    exprDevEigen(NumberOfExpressions)
{
    load(path);
}

void FaceModel::load(const std::string &path) {
    std::string pathBasisShape = path + "/" + filenameBasisShape;
    std::string pathBasisExpression = path + "/" + filenameBasisExpression;
    std::string pathStdDevShape = path + "/" + filenameStdDevShape;
    std::string pathStdDevExpression = path + "/" + filenameStdDevExpression;

    auto shapeBasisCPU = new float4[nVertices * NumberOfEigenvectors];
    auto expressionBasisCPU = new float4[nVertices * NumberOfExpressions];
    loadVector(pathBasisShape, (float*)shapeBasisCPU, 4 * nVertices * NumberOfEigenvectors);
    loadVector(pathBasisExpression, (float*)expressionBasisCPU, 4 * nVertices * NumberOfExpressions);

    auto shapeDevCPU = new float[NumberOfEigenvectors];
    auto expressionDevCPU = new float[NumberOfExpressions];
    loadVector(pathStdDevShape, shapeDevCPU, NumberOfEigenvectors);
    loadVector(pathStdDevExpression, expressionDevCPU, NumberOfExpressions);

    convertBasis(shapeBasisEigen, nVertices, NumberOfEigenvectors, shapeBasisCPU);
    convertBasis(exprBasisEigen, nVertices, NumberOfExpressions, expressionBasisCPU);

    convertDeviations(shapeDevEigen, NumberOfEigenvectors, shapeDevCPU);
    convertDeviations(exprDevEigen, NumberOfExpressions, expressionDevCPU);

    loadAverageMesh();

    delete[] shapeBasisCPU;
    delete[] expressionBasisCPU;
    delete[] shapeDevCPU;
    delete[] expressionDevCPU;
}

void FaceModel::forwardPass(const Eigen::VectorXd& alpha,
                           const Eigen::VectorXd& delta,
                           Eigen::VectorXf& vertices_out)
{
    vertices_out = shapeBasisEigen * alpha.cast<float>() + exprBasisEigen * delta.cast<float>();
}

int FaceModel::loadAverageMesh() {
    OpenMesh::IO::Options ropt;

    // Set input options
    ropt += OpenMesh::IO::Options::VertexColor;

    m_avg_mesh.request_vertex_colors();

    // assure we have vertex normals
    if (!m_avg_mesh.has_vertex_colors())
    {
      std::cerr << "ERROR: Standard vertex property 'Colors' not available for average mesh!\n";
      return 1;
    }

    if (!OpenMesh::IO::read_mesh(m_avg_mesh, FILENAME_AVG_MESH, ropt)) {
        std::cerr << "ERROR: Could not load " << FILENAME_AVG_MESH << "\n";
        return 1;
    }

    // The average mesh provided by Justus Thies must be scaled by a factor
    // AVG_MESH_SCALE
    for (FaceMesh::VertexIter v_it = m_avg_mesh.vertices_begin();
         v_it != m_avg_mesh.vertices_end(); ++v_it)
    {
      m_avg_mesh.set_point( *v_it, m_avg_mesh.point(*v_it) * SCALE_AVG_MESH );
    }

    // Save synth mesh as average mesh in the beginning.
    m_synth_mesh = m_avg_mesh;

    return 0;
}

int FaceModel::synthesizeModel(const Eigen::VectorXf& diff_vertices) {
    FaceMesh synth_mesh = m_avg_mesh;

    size_t i = 0;

    for (FaceMesh::VertexIter v_it = synth_mesh.vertices_begin();
         v_it != synth_mesh.vertices_end(); ++v_it)
    {
      // Add diff_vertives (change in the vertex coordinates based
      // on the PCA model) to the existing vertices.
      synth_mesh.set_point( *v_it, synth_mesh.point(*v_it)
                          + FaceMesh::Point(diff_vertices(i), diff_vertices(i+1),
                                          diff_vertices(i+2)));

      i += 3;
    }

    m_synth_mesh = synth_mesh;

    return 0;
}

int FaceModel::writeSynthesizedModel() {
    OpenMesh::IO::Options wopt;
    wopt += OpenMesh::IO::Options::VertexColor;

    try
    {
      if (!OpenMesh::IO::write_mesh(m_synth_mesh, FILENAME_OUT_SYNTH_MESH, wopt))
      {
        std::cerr << "Cannot write mesh to file '" << FILENAME_OUT_SYNTH_MESH
                  << "'" << std::endl;
      }
    }
    catch(std::exception& x)
    {
      std::cerr << x.what() << std::endl;
      return 1;
    }

    return 0;
}

void FaceModel::loadVector(const std::string &filename, float *res, unsigned int length)
{
    std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
    if (!in)
    {
        std::cout << "ERROR:\tCan not open file: " << filename << std::endl;
        while (1);
    }
    unsigned int numberOfEntries;
    in.read((char*)&numberOfEntries, sizeof(unsigned int));
    if (length == 0) length = numberOfEntries;
    in.read((char*)(res), length * sizeof(float));

    in.close();
}

void FaceModel::progressBar(const char* str, int num)
{
    return;
}

float* FaceModel::loadEigenvectors(const std::string &filename, unsigned int components, unsigned int numberOfEigenvectors)
{
    float *res = new float[components*numberOfEigenvectors];

    for (unsigned int i = 0; i < numberOfEigenvectors; i++)
    {
        progressBar("Load Model Basis:", i / float(numberOfEigenvectors));
        std::stringstream ss;
        ss << filename << i << ".vec";

        loadVector(ss.str().c_str(), &(res[components*i]), 0);
    }
    progressBar("Load Model Basis:", 1.0f);

    return res;
}

int FaceModel::convertBasis(Eigen::MatrixXf& basis, int nVertices, int numCols, float4* basisCPU)
{
    for(int i = 0; i < numCols; i++) {
        for (int j = 0; j < nVertices; j++) {
            basis(j * 3 + 0, i) = basisCPU[i * nVertices + j].x;
            basis(j * 3 + 1, i) = basisCPU[i * nVertices + j].y;
            basis(j * 3 + 2, i) = basisCPU[i * nVertices + j].z;
        }
    }

    return 0;
}

int FaceModel::convertDeviations(Eigen::VectorXf& devs, int num_dims, float* devCPU)
{
    for(int i = 0; i < num_dims; i++) {
            devs(i, 0) = devCPU[i];
    }

    return 0;
}