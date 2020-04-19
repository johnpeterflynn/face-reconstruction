#ifndef CONFIG_H
#define CONFIG_H

#include <string>

constexpr int NumberOfEigenvectors = 160;
constexpr int NumberOfExpressions = 76;
constexpr int nVertices = 53490;

// Number of parameters in alpha and delta to solve for
constexpr int NUM_PARAMS_ALPHA = 40;
constexpr int NUM_PARAMS_DELTA = 30;

const std::string PATH_SCANNED_MESH = "Face_Raw.ply";
const std::string PATH_SCANNED_LANDMARKS = "HD_Face.ply";
const std::string PATH_DEFAULT_SCANNED_MESG = "scan/kinectdata_example.off";

const std::string PATH_MESH = "mesh/";
const std::string PATH_BASIS = "basis/";

const std::string FILENAME_AVG_MESH = PATH_MESH + "averageMesh.off";
const std::string FILENAME_AVG_OPT_MESH = PATH_MESH + "averageMesh_blackface.off";
const std::string PATH_OUT_SYNTH_MESH = "synthesizedMesh.off";

const std::string FILENAME_BASIS_SHAPE = PATH_BASIS + "ShapeBasis.matrix";
const std::string FILENAME_BASIS_EXPRESSION = PATH_BASIS + "ExpressionBasis.matrix";
const std::string FILENAME_STDDEV_SHAPE = PATH_BASIS + "StandardDeviationShape.vec";
const std::string FILENAME_STDDEV_EXPRESSION = PATH_BASIS + "StandardDeviationExpression.vec";


#endif // CONFIG_H
