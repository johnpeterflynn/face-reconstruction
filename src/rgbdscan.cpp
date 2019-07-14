#include "rgbdscan.h"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/getopt.h>


RGBDScan::RGBDScan(const std::string& path) {
    loadMesh(path);
}

void RGBDScan::loadMatchIndices() {
    // Dummy data for now
    m_match_indices[2000] = 2000;
    m_match_indices[10000] = 10000;
    m_match_indices[25000] = 25000;
    m_match_indices[30000] = 30000;
    m_match_indices[50000] = 50000;
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
      return 1;
    }

    if (!OpenMesh::IO::read_mesh(m_scanned_mesh, path, ropt)) {
        std::cerr << "ERROR: Could not load " << path << "\n";
        return 1;
    }

    return 0;
}