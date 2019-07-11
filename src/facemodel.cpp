#include "facemodel.h"

#include <iostream>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/getopt.h>


constexpr const char* FILENAME_AVG_MESH = "../models/averageMesh.off";
constexpr const char* FILENAME_OUT_SYNTH_MESH = "synthesizedMesh.off";

FaceModel::FaceModel()
{
    loadAverageMesh();
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

    return 0;
}

int FaceModel::writeSynthesizedModel(const Eigen::VectorXf& diff_vertices) {
    FaceMesh synth_mesh = m_avg_mesh;

    OpenMesh::IO::Options wopt;
    wopt += OpenMesh::IO::Options::VertexColor;

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

    try
    {
      if (!OpenMesh::IO::write_mesh(synth_mesh, FILENAME_OUT_SYNTH_MESH, wopt))
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