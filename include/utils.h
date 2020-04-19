#ifndef UTILS_H
#define UTILS_H

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

typedef OpenMesh::TriMesh_ArrayKernelT<> FaceMesh;


static int loadMesh(const std::string& path, FaceMesh& mesh) {
    OpenMesh::IO::Options ropt;

    // Set input options
    ropt += OpenMesh::IO::Options::VertexColor;

    mesh.request_vertex_colors();

    // assure we have vertex normals
    if (!mesh.has_vertex_colors())
    {
      std::cerr << "ERROR: Standard vertex property 'Colors' not available for average mesh!\n";
      return 1;
    }

    if (!OpenMesh::IO::read_mesh(mesh, path, ropt)) {
        std::cerr << "ERROR: Could not load " << path << "\n";
        return 1;
    }

    return 0;
}

static double compareMeshes(FaceMesh &synth_mesh, FaceMesh &ideal_mesh) {
    // Compute the cumulative offset of all vertices between two meshes of the
    // same vertex number

    double sum = 0;

    FaceMesh::VertexIter vi_it = ideal_mesh.vertices_begin();
    for (FaceMesh::VertexIter vs_it = synth_mesh.vertices_begin();
         vs_it != synth_mesh.vertices_end(); ++vs_it)
    {
      FaceMesh::Point ps = synth_mesh.point(*vs_it);
      Eigen::Vector3d vs(ps[0], ps[1], ps[2]);
      FaceMesh::Point pi = ideal_mesh.point(*vi_it);
      Eigen::Vector3d vi(pi[0], pi[1], pi[2]);
      Eigen::Vector3d diff = vs - vi;

      sum += sqrt(diff.dot(diff));
      ++vi_it;
    }

    return sum;
}

static double compareToIdeal(FaceMesh &synth_mesh, const std::string file_ideal) {
    FaceMesh ideal_mesh;
    loadMesh(file_ideal, ideal_mesh);
    double sum_diff = compareMeshes(synth_mesh, ideal_mesh);

    return sum_diff;
}

#endif // UTILS_H
