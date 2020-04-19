#include <iostream>
#include <chrono>
#include <fstream>
#include <string>

#include <boost/program_options.hpp>
#include <Eigen/Dense>

#include "pch.h"
#include "facemodel.h"
#include "rgbdscan.h"
#include "facesolver.h"
#include "config.h"
#include "utils.h"

namespace po = boost::program_options;


void run(FaceSolver &face_solver, RGBDScan &face_scan, FaceModel &face_model,
         const std::string filename_out)
{
    Eigen::VectorXd alpha = Eigen::VectorXd::Zero(NumberOfEigenvectors);
    Eigen::VectorXd delta = Eigen::VectorXd::Zero(NumberOfExpressions);

    Sophus::SE3d T_xy;

    auto start = std::chrono::high_resolution_clock::now();
    face_solver.solve(face_model, face_scan, alpha, delta, T_xy);
    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << "Writing synthesized model to file\n";
    face_model.writeSynthesizedModel(filename_out, alpha, delta, T_xy);

    // Print alpha (geometry parameter we solved for)
    //std::cout << alpha << std::endl;

    std::cout << "Transform" << std::endl;
    std::cout << T_xy.matrix() << std::endl;

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "Face solver took: "
             << duration.count() / 1000000.0 << " seconds" << std::endl;

    // Compute sum of vertex displacements between this model and the ideal
    FaceMesh synth_mesh = face_model.synthesizeModel(alpha, delta, T_xy);
    double sum_diff = compareToIdeal(synth_mesh, "synth_ideal.off");
    std::cout << "Difference: " << sum_diff << "\n";
}

int main(int argc, char *argv[]) {
    std::ifstream config_File("config.ini");

    // Declare the supported options.
    po::options_description generic("Generic options");
    generic.add_options()
        ("help", "produce help message")
    ;

    // Declare the required options.
    po::options_description required("Required options");
    required.add_options()
        ("scan", po::value<std::string>(), "path to scan file (.off)")
        ("corr", po::value<std::string>(), "path to corr file (.corr)")
        ("out", po::value<std::string>(), "output name of synthesized mesh (.off)")
    ;

    // Declare the optional options.
    po::options_description optional("Optional");
    optional.add_options()
        ("model-path", po::value<std::string>(), "path to model contents (basis, mesh, landmarks)")
    ;

    // Declare the supported options.
    po::options_description face("FaceSolver");
    face.add_options()
        ("huber", po::value<double>(), "set Huber parameter")
        ("geo-reg", po::value<double>(), "set geometric regulatization constant")
        ("knn-dist-thresh", po::value<double>(), "set threshold to accept scan-model matches in meters")
        ("num-iters", po::value<int>(), "set number of knn -> loss minimization iterations")
        ("frac-used-vertices", po::value<double>(), "set fraction of vertices to randomly sample for optimization")
        ("ignore-borders", po::value<bool>(), "ignore borders of scan")
    ;

    po::options_description config("Config options");
    config.add(face).add(optional);
    po::options_description visible("Allowed options");
    visible.add(generic).add(required).add(config);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, visible), vm);
    po::store(po::parse_config_file(config_File, config), vm);
    po::notify(vm);

    // Display help and exit
    if (vm.count("help")) {
        std::cout << visible << "\n";
        return 1;
    }

    // Enforce inclusion of required parameters
    if (!vm.count("scan")) {
        std::cout << "Scan file not set.\n";
        return 1;
    }
    if (!vm.count("corr")) {
        std::cout << "Corr file not set.\n";
        return 1;
    }
    if (!vm.count("out")) {
        std::cout << "Output file not set.\n";
        return 1;
    }

    // Build model, scan and solver objects
    FaceModel face_model(vm["model-path"].as<std::string>());

    RGBDScan face_scan(vm["scan"].as<std::string>(),
                       vm["corr"].as<std::string>());

    FaceSolver face_solver(vm["geo-reg"].as<double>(),
                           vm["huber"].as<double>(),
                           vm["knn-dist-thresh"].as<double>(),
                           vm["num-iters"].as<int>(),
                           vm["frac-used-vertices"].as<double>(),
                           vm["ignore-borders"].as<bool>());

    // Execute solver with tests
    run(face_solver, face_scan, face_model, vm["out"].as<std::string>());

    return 0;
}
