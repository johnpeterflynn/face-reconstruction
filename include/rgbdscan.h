#ifndef RGBDSCAN_H
#define RGBDSCAN_H

#include <iostream>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>


typedef OpenMesh::TriMesh_ArrayKernelT<>  ScanMesh;

class RGBDScan {
public:
    RGBDScan(const std::string& path);

    void loadMatchIndices();
    std::map<int, int> getMatchIndices();

private:
    int loadMesh(const std::string& path);

public:
    ScanMesh m_scanned_mesh;
    std::map<int, int> m_match_indices;
};

#endif // RGBDSCAN_H
