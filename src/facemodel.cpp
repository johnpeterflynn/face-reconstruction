#include "facemodel.h"

#include <iostream>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/getopt.h>


constexpr const char* FILENAME_AVG_MESH = "../models/averageMesh.off";

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
    for (AvgMesh::VertexIter v_it = m_avg_mesh.vertices_begin();
         v_it != m_avg_mesh.vertices_end(); ++v_it)
    {
      m_avg_mesh.set_point( *v_it, m_avg_mesh.point(*v_it) * SCALE_AVG_MESH );
    }
}
