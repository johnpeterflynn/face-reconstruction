#include "rgbdscan.h"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/getopt.h>


RGBDScan::RGBDScan(const std::string& path) {
    loadMesh(path);
}

void RGBDScan::loadMatchIndices() {
    // Dummy data for now
/*
    m_match_indices[12280] = 12280;
    m_match_indices[4280] = 4280;
    m_match_indices[8319] = 8319;
    m_match_indices[48246] = 48246;

    m_match_indices[482] = 482;
    m_match_indices[2881] = 2881;
    m_match_indices[50000] = 50000;
    m_match_indices[31000] = 31000;
*/

    m_match_indices[12280] = 9734;
    m_match_indices[4280] = 10045;
    m_match_indices[8319] = 12539;
    m_match_indices[48246] = 19679;

}

std::map<int, int> RGBDScan::getMatchIndices() {
    return m_match_indices;
}

int RGBDScan::loadMesh(const std::string& path) {
    OpenMesh::IO::Options ropt;

    // Set input options
    ropt += OpenMesh::IO::Options::VertexColor;

    m_scanned_mesh.request_vertex_colors();
    m_scanned_mesh.request_vertex_normals();

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

    if ( !ropt.check( OpenMesh::IO::Options::VertexNormal ) )
    {
      // we need face normals to update the vertex normals
      m_scanned_mesh.request_face_normals();
      // let the mesh update the normals
      m_scanned_mesh.update_normals();
      // dispose the face normals, as we don't need them anymore
      m_scanned_mesh.release_face_normals();
    }

    return 0;
}