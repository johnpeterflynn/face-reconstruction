#ifndef CONFIG_H
#define CONFIG_H

constexpr int NumberOfEigenvectors = 160;
constexpr int NumberOfExpressions = 76;
constexpr int nVertices = 53490;

// Number of parameters in alpha and delta to solve for
constexpr int NUM_PARAMS_ALPHA = 10;
constexpr int NUM_PARAMS_DELTA = 10;

constexpr const char* MODEL_PATH = "../models";
constexpr const char* FILENAME_SCANNED_MESH = "../testData/fakekinectdata.off";

constexpr const char* FILENAME_AVG_MESH = "../models/averageMesh.off";
constexpr const char* FILENAME_OUT_SYNTH_MESH = "synthesizedMesh.off";

constexpr const char* filenameBasisShape = "ShapeBasis.matrix";
constexpr const char* filenameBasisExpression = "ExpressionBasis.matrix";
constexpr const char* filenameStdDevShape = "StandardDeviationShape.vec";
constexpr const char* filenameStdDevExpression = "StandardDeviationExpression.vec";


#endif // CONFIG_H
