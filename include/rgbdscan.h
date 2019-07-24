#ifndef RGBDSCAN_H
#define RGBDSCAN_H

#include <iostream>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>


typedef OpenMesh::TriMesh_ArrayKernelT<>  ScanMesh;

class RGBDScan {
public:
    RGBDScan(const std::string& face_path,
             const std::string& landmark_path);

    std::map<int, int> getMatchIndices();

private:
    void loadMatchIndices(const std::string& landmark_path);
    int loadMesh(const std::string& face_path);

public:
    ScanMesh m_scanned_mesh;
    ScanMesh m_landmarks;
    std::map<int, int> m_match_indices;
};

#endif // RGBDSCAN_H
