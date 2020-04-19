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
#include"config.h"

namespace po = boost::program_options;


void run(FaceSolver &face_solver, RGBDScan &face_scan, FaceModel &face_model)
{
    Eigen::VectorXd alpha = Eigen::VectorXd::Zero(NumberOfEigenvectors);
    Eigen::VectorXd delta = Eigen::VectorXd::Zero(NumberOfExpressions);

    Sophus::SE3d T_xy;

    auto start = std::chrono::high_resolution_clock::now();
    face_solver.solve(face_model, face_scan, alpha, delta, T_xy);
    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << "Writing synthesized model to file\n";
    face_model.writeSynthesizedModel(alpha, delta, T_xy);

    // Print alpha (geometry parameter we solved for)
    //std::cout << alpha << std::endl;

    std::cout << "Transform" << std::endl;
    std::cout << T_xy.matrix() << std::endl;

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "Face solver took: "
             << duration.count() / 1000000.0 << " seconds" << std::endl;

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
        ("data-path", po::value<std::string>(), "path to data folder")
        ("out", po::value<std::string>(), "output name of synthesized mesh")
    ;

    // Declare the supported options.
    po::options_description config("FaceSolver");
    config.add_options()
        ("huber", po::value<double>(), "set Huber parameter")
        ("geo-reg", po::value<double>(), "set geometric regulatization constant")
        ("knn-dist-thresh", po::value<double>(), "set threshold to accept scan-model matches in meters")
        ("num-iters", po::value<int>(), "set number of knn -> loss minimization iterations")
        ("frac-used-vertices", po::value<double>(), "set fraction of vertices to randomly sample for optimization")
        ("ignore-borders", po::value<bool>(), "ignore borders of scan")
    ;

    po::options_description visible("Allowed options");
    visible.add(generic).add(required).add(config);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, visible), vm);
    po::store(po::parse_config_file(config_File, config), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << visible << "\n";
        return 1;
    }

    if (vm.count("data-path")) {
     //<< vm["data-path"].as<int>() << ".\n";
    } else {
        std::cout << "Data path not set.\n";
    }

    if (vm.count("out")) {
     //<< vm["data-path"].as<int>() << ".\n";
    } else {
        std::cout << "Output file not set.\n";
    }

    FaceModel face_model(vm["data-path"].as<std::string>() + MODEL_PATH);
    RGBDScan face_scan(vm["data-path"].as<std::string>() + FILENAME_SCANNED_MESH, vm["data-path"].as<std::string>() + FILENAME_SCANNED_LANDMARKS);
    FaceSolver face_solver(vm["geo-reg"].as<double>(), vm["huber"].as<double>(),
            vm["knn-dist-thresh"].as<double>(), vm["num-iters"].as<int>(),
            vm["frac-used-vertices"].as<double>(),
            vm["ignore-borders"].as<bool>());

    run(face_solver, face_scan, face_model);

    return 0;
}
