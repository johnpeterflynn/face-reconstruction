#include "facesolver.h"

#include <ceres/ceres.h>
#include "nabo/nabo.h"

#include "local_parameterization_se3.hpp"
#include "config.h"


struct ReconstructionCostFunctor {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  ReconstructionCostFunctor(const Eigen::Matrix<double, 3, 1>& v_face_avg_in,
                            const Eigen::Matrix<double, 3, 1>& v_scan_in,
                            const Eigen::Matrix<double, 3, NUM_PARAMS_ALPHA>& shapeBasisEigenRow_in,
                            const Eigen::Matrix<double, 3, NUM_PARAMS_DELTA>& exprBasisEigenRow_in,
                            double weight_in = 1.0)
      : v_face_avg(v_face_avg_in), v_scan(v_scan_in),
        shapeBasisEigenRow(shapeBasisEigenRow_in),
        exprBasisEigenRow(exprBasisEigenRow_in), weight(weight_in) {}

  template <class T>
  bool operator()(T const* const salpha,/* T const* const sdelta,*/
                  T const* const sT_xy,
                  T* sresiduals) const {
    Eigen::Map<Eigen::Matrix<T, NUM_PARAMS_ALPHA, 1> const> const alpha(salpha);
    //Eigen::Map<Eigen::Matrix<T, NUM_PARAMS_DELTA, 1> const> const delta(sdelta);
    Eigen::Map<Sophus::SE3<T> const> const T_xy(sT_xy);
    Eigen::Map<Eigen::Matrix<T, 3, 1>> residuals(sresiduals);

    Eigen::Matrix<T, 3, 1> v_model =
            (v_face_avg.cast<T>() + shapeBasisEigenRow.cast<T>() * alpha /*+ exprBasisEigenRow.cast<T>() * delta*/);

    Eigen::Matrix<T, 3, 1> v_scan_est = T_xy * v_model;

    residuals = sqrt(T(weight)) * (v_scan.cast<T>() - v_scan_est);

    return true;
  }

  const Eigen::Matrix<double, 3, 1> v_face_avg;
  const Eigen::Matrix<double, 3, 1> v_scan;

  // TODO: Again, these are copied so need another way.
  const Eigen::Matrix<double, 3, NUM_PARAMS_ALPHA> shapeBasisEigenRow;
  const Eigen::Matrix<double, 3, NUM_PARAMS_DELTA> exprBasisEigenRow;

  const double weight; // Weight of residual
};

struct GeometryRegularizationCostFunctor {
  GeometryRegularizationCostFunctor(double lambda_in, double stddev_in, int index_in)
      : lambda(lambda_in), index(index_in), stddev(stddev_in) {}

  template <class T>
  bool operator()(T const* const params_in, T* residuals) const {
    const T param = params_in[index];

    residuals[0] = sqrt(lambda) * param / stddev;

    return true;
  }

private:
    const double lambda;
    const int index;
    const double stddev;
};

FaceSolver::FaceSolver(double geo_regularization, int num_iterations) :
    m_geo_regularization(geo_regularization),
    m_num_iterations(num_iterations)
{
}

void FaceSolver::runCeres(const MyMesh& avg_face_mesh, const MyMesh& scanned_mesh,
              const Eigen::MatrixXi& indices,
              const std::set<int>& matches,
              const Eigen::MatrixXf& shapeBasisEigen,
              const Eigen::MatrixXf& exprBasisEigen,
              const Eigen::VectorXf& shapeDevEigen,
              const Eigen::VectorXf& exprDevEigen,
              Eigen::VectorXd& alpha, Eigen::VectorXd& delta,
              Sophus::SE3d& T_xy) {
    ceres::Problem problem;

    // Model transformation
    //Eigen::Matrix3f R;
    //Eigen::Vector3f t;

    problem.AddParameterBlock(T_xy.data(),
                              Sophus::SE3d::num_parameters,
                              new Sophus::test::LocalParameterizationSE3);

    for (MyMesh::VertexIter v_it = avg_face_mesh.vertices_begin();
         v_it != avg_face_mesh.vertices_end(); ++v_it)
    {
        // Get index of average face vertex
        int v_idx = v_it->idx();
/*
        if ((v_idx % 50) != 0) {
            continue;
        }
*/

        // Constants in optimization

        // Get vertex itself for average face mesh
        MyMesh::Point p3_face_avg = avg_face_mesh.point(*v_it);
        Eigen::Vector3d v3_face_avg(p3_face_avg[0], p3_face_avg[1], p3_face_avg[2]);

        // Get nearest neighbor vertex in scanned mesh
        MyMesh::VertexHandle v_scan_handle(indices(0, v_idx));
        MyMesh::Point p3_scan_nn = scanned_mesh.point(v_scan_handle);
        Eigen::Vector3d v3_scan_nn(p3_scan_nn[0], p3_scan_nn[1], p3_scan_nn[2]);

        // TODO:These are being copied but pointers should be used instead
        Eigen::Matrix<double, 3, NUM_PARAMS_ALPHA> shapeBasisEigenRows
                = shapeBasisEigen.block<3, NUM_PARAMS_ALPHA>(3 * v_idx, 0).cast<double>();
        Eigen::Matrix<double, 3, NUM_PARAMS_DELTA> exprBasisEigenRows
                = exprBasisEigen.block<3, NUM_PARAMS_DELTA>(3 * v_idx, 0).cast<double>();

        // Add more weight to sparse correspondences
        double weight = 1.0;
        if (matches.find(v_idx) != matches.end()) {
            weight = 100.0;
        }

        // For each entry in a vector v3
        problem.AddResidualBlock(
            // <dim of residual, dim of alpha, dim of delta>
            new ceres::AutoDiffCostFunction<ReconstructionCostFunctor, 3,
                    NUM_PARAMS_ALPHA/*, NUM_PARAMS_DELTA*/,
                    Sophus::SE3d::num_parameters>(
                new ReconstructionCostFunctor(v3_face_avg, v3_scan_nn,
                                              shapeBasisEigenRows,
                                              exprBasisEigenRows, weight)),
            nullptr, alpha.data()/*, delta.data()*/, T_xy.data());

        //std::cout << "Adding residual: " << v_idx << "/" << avg_face_mesh.n_vertices() << "\n";
    }

    // Q: TODO: Is it a waste here to pass the whole alpha and delta? Would it
    // be better to make one large residual for alpha and delta?
    for (int i = 0; i < NUM_PARAMS_ALPHA; i++) {
        problem.AddResidualBlock(
            // <dim of residual, dim of alpha, dim of delta>
            new ceres::AutoDiffCostFunction<GeometryRegularizationCostFunctor, 1,
                    NUM_PARAMS_ALPHA>(
                new GeometryRegularizationCostFunctor(m_geo_regularization, (double)shapeDevEigen(i), i)),
            nullptr, alpha.data());

    }

    /*
    for (int i = 0; i < NUM_PARAMS_DELTA; i++) {
        problem.AddResidualBlock(
            // <dim of residual, dim of alpha, dim of delta>
            new ceres::AutoDiffCostFunction<GeometryRegularizationCostFunctor, 1,
                    NUM_PARAMS_DELTA>(
                new GeometryRegularizationCostFunctor(m_geo_regularization, (double)exprDevEigen(i), i)),
            nullptr, delta.data());
        // Q: TODO: Is it a waste here to pass the whole delta?
    }
*/
    ceres::Solver::Options options;
    options.max_num_iterations = 1;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//DENSE_QR;
    options.minimizer_progress_to_stdout = true;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << std::endl;
    std::cout << summary.FullReport() << std::endl;
}

void FaceSolver::calculate_knn(const Eigen::MatrixXf& M, const Eigen::MatrixXf& q,
                  Eigen::MatrixXi& indices)
{
    const int K = indices.rows();

    Nabo::NNSearchF* nns = Nabo::NNSearchF::createKDTreeLinearHeap(M);
    //Eigen::VectorXi indices(K);
    Eigen::MatrixXf dists2(K, q.cols());

    // ALLOW_SELF_MATCH appears to be necessary to match vertices that are
    // exactly the same.
    nns->knn(q, indices, dists2, K,  0, Nabo::NNSearchF::ALLOW_SELF_MATCH);
    delete nns;

    return;
}

int FaceSolver::knn_model_to_scan(const FaceModel& face_model, const MyMesh& scanned_mesh, int K, Eigen::MatrixXi& indices) {
    Eigen::MatrixXf FM(3, face_model.m_synth_mesh.n_vertices());
    Eigen::MatrixXf SM(3, scanned_mesh.n_vertices());

    meshToMatrix(face_model.m_synth_mesh, FM);
    meshToMatrix(scanned_mesh, SM);

    std::cout << "Starting KNN\n";

    // Indices for  for each column vector in FN
    indices.resize(K, FM.cols());

    calculate_knn(SM, FM, indices);

    std::cout << "Finished KNN\n";
}

void FaceSolver::meshToMatrix(const MyMesh& mesh, Eigen::MatrixXf& M_out) {
    for (MyMesh::VertexIter v_it = mesh.vertices_begin();
         v_it != mesh.vertices_end(); ++v_it)
    {
      MyMesh::Point p3 = mesh.point(*v_it);
      Eigen::Vector3f v3(p3[0], p3[1], p3[2]);

      M_out.col(v_it->idx()) = v3;
    }
}

void FaceSolver::solve(FaceModel& face_model, RGBDScan face_scan,
                       Eigen::VectorXd& alpha, Eigen::VectorXd& delta,
                       Sophus::SE3d& T_xy) {
    const int K = 1;
    Eigen::MatrixXi indices;

    std::cout << "Optimization starting: " << std::endl;
    for (int i = 0; i < m_num_iterations; i++) {
        knn_model_to_scan(face_model, face_scan.m_scanned_mesh, K, indices);

        std::map<int, int> match_indices = face_scan.getMatchIndices();

        // Replace KNN mathes with manually selected matches
        std::set<int> matches;
        for (auto entry : match_indices) {
            matches.insert(entry.first);
            indices(0, entry.first) = entry.second;
        }

        std::cout << "Running ceres: " << std::endl;
        runCeres(face_model.m_synth_mesh, face_scan.m_scanned_mesh, indices, matches,
                 face_model.shapeBasisEigen, face_model.exprBasisEigen,
                 face_model.shapeDevEigen, face_model.exprDevEigen,
                 alpha, delta, T_xy);

        std::cout << "Ceres finished: " << std::endl;
    }

    std::cout << "Optimization finished: " << std::endl;
}