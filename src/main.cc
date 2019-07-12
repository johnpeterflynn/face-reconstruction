// ConsoleApplication1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <Eigen/Dense>
#include <vector>

//#include <opencv2/core/core.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
//typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyPolyMesh;
#include <OpenMesh/Tools/Utils/getopt.h>
#include <ceres/ceres.h>

#include "nabo/nabo.h"

#include "facemodel.h"

constexpr int NumberOfEigenvectors = 160;
constexpr int NumberOfExpressions = 76;
constexpr int nVertices = 53490;

int surfaceNormalsTest(void)
{

	MyMesh mesh;
	
    mesh.request_vertex_normals();

    //mesh.request_face_normals();

    if (!OpenMesh::IO::read_mesh(mesh, "../testData/kinectdata.off")) {
		std::cerr << "read error\n";
        exit(1);
	}
  
  OpenMesh::IO::Options opt;
  if (!opt.check(OpenMesh::IO::Options::VertexNormal )) {
    // we need face normals to update the vertex normals
    mesh.request_face_normals();
    
    // let the mesh update the normals
    mesh.update_normals();
    
    // dispose the face normals, as we don't need them anymore
    mesh.release_face_normals();
  }

  // move all vertices one unit length along it's normal direction
  for (MyMesh::VertexIter v_it = mesh.vertices_begin();
       v_it != mesh.vertices_end(); ++v_it)
  {
    std::cout << "Vertex #" << *v_it << ": " << mesh.point( *v_it );
    mesh.set_point( *v_it, mesh.point(*v_it) + 0.01 * mesh.normal(*v_it) );
    std::cout << " moved to " << mesh.point( *v_it ) << std::endl;
  }
  // don't need the normals anymore? Remove them!

  mesh.release_vertex_normals();

  // just check if it really works
  if (mesh.has_vertex_normals())
  {
    std::cerr << "Ouch! ERROR! Shouldn't have any vertex normals anymore!\n";
    return 1;
  }

  if (!OpenMesh::IO::write_mesh(mesh, "out.off")) 
  {
    std::cerr << "write error\n";
    exit(1);
  }

  return 0;
}

void calculate_knn(const Eigen::MatrixXf& M, const Eigen::MatrixXf& q,
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

constexpr const char* FILENAME_SCANNED_MESH = "../testData/fakekinectdata.off";
//"../testData/kinectdata.off";

int loadScannedMesh(MyMesh& scanned_mesh) {
    OpenMesh::IO::Options ropt;

    // Set input options
    ropt += OpenMesh::IO::Options::VertexColor;

    scanned_mesh.request_vertex_colors();

    // assure we have vertex normals
    if (!scanned_mesh.has_vertex_colors())
    {
      std::cerr << "ERROR: Standard vertex property 'Colors' not available for average mesh!\n";
      return 1;
    }

    if (!OpenMesh::IO::read_mesh(scanned_mesh, FILENAME_SCANNED_MESH, ropt)) {
        std::cerr << "ERROR: Could not load " << FILENAME_SCANNED_MESH << "\n";
        return 1;
    }

    return 0;
}

void projectionTest(const MyMesh& mesh, std::string s, int corr_index,
                    Eigen::MatrixXf& M_out) {
    const unsigned int IMG_WIDTH = 640;
    const unsigned int IMG_HEIGHT = 480;

    int index = 0;

    for (MyMesh::VertexIter v_it = mesh.vertices_begin();
         v_it != mesh.vertices_end(); ++v_it)
    {
      MyMesh::Point p3 = mesh.point(*v_it) * 2000.0;

      // 3D-2D projection

      float depth = 1;//p3[2];
      Eigen::Vector3f v3(p3[0], p3[1], depth);

      float f = 1; // focal length
      float ox = IMG_WIDTH / 2.0;
      float oy = IMG_HEIGHT / 2.0;
      Eigen::Matrix<float, 2, 3> Intrinsics;
      Intrinsics << f, 0, ox,
                    0, f, oy;
      Eigen::Matrix<float, 2, 3> Pi = Intrinsics / depth;

      Eigen::Vector2f v2 = Pi * v3;

      // Add 2d projection v2 to M_out
      M_out.col(index) = v2;

      // Let's allow off-screen projections for the moment, since they
      // presumably won't be selected by KNN anyway.
      /*
      if (0 <= v2(0) && v2(0) < IMG_WIDTH && 0 <= v2(1) && v2(1) < IMG_HEIGHT) {
        cv::Vec3b vcolor;
        if (index == corr_index) {
          vcolor = cv::Vec3b(255, 255, 0);
        }
        else {
          vcolor = cv::Vec3b(0, 0, 255);
        }
        M.at<cv::Vec3b>(cv::Point(round(v2(0)),round(v2(1)))) = vcolor;
      }
      */

      index++;
    }

    //cv::imwrite( "./images/" + s + ".jpg", M);
}

void writeMatrixToImg(const Eigen::MatrixXf& matrix, cv::Mat& img,
                 cv::Vec3b vcolor) {
    //if (0 <= v2(0) && v2(0) < IMG_WIDTH && 0 <= v2(1) && v2(1) < IMG_HEIGHT) {
    for (int c = 0; c < matrix.cols(); c++) {
        Eigen::Vector2f v2 = matrix.col(c);
      img.at<cv::Vec3b>(cv::Point(round(v2(0)),round(v2(1)))) = vcolor;
    }
    //}
}

void meshToMatrix(const MyMesh& mesh, Eigen::MatrixXf& M_out) {
    for (MyMesh::VertexIter v_it = mesh.vertices_begin();
         v_it != mesh.vertices_end(); ++v_it)
    {
      MyMesh::Point p3 = mesh.point(*v_it);
      Eigen::Vector3f v3(p3[0], p3[1], p3[2]);

      M_out.col(v_it->idx()) = v3;
    }
}

int knn_test(const FaceModel& face_model, const MyMesh& scanned_mesh, int K, Eigen::MatrixXi& indices) {
    Eigen::MatrixXf FM(3, face_model.m_synth_mesh.n_vertices());
    Eigen::MatrixXf SM(3, scanned_mesh.n_vertices());

    meshToMatrix(face_model.m_synth_mesh, FM);
    meshToMatrix(scanned_mesh, SM);
/*
    std::cout << "Num FM rows, cols: " << FM.rows() << ", " << FM.cols() << "\n";
    std::cout << "Num SM rows, cols: " << SM.rows() << ", " << SM.cols() << "\n";

    // Sample the matrices to see if they look okay.
    std::cout << "FM first 10 values\n";
    std::cout << FM.block<3,10>(0,0) << "\n";
    std::cout << "SM first 10 values\n";
    std::cout << SM.block<3,10>(0,0) << "\n";
*/
    std::cout << "Starting KNN\n";

    // Indices for  for each column vector in FN
    //Eigen::MatrixXi indices;
    indices.resize(K, FM.cols());

    calculate_knn(SM, FM, indices);

    std::cout << "Finished KNN\n";

    //std::cout << "incides:\n";
    //std::cout << indices << "\n";
/*
    for (int i = 0; i < indices.cols(); i++) {
        std::cout << "Face " << i << ": " << FM.col(i).transpose() << std::endl;
        std::cout << "Scan " << indices(0,i) << ": " << SM.col(indices(0,i)).transpose() << std::endl;
    }
    */

/*
    std::cout << "Writing results to image\n";

    //## Visualize the projections and nearest neighbors
    const unsigned int IMG_WIDTH = 640;
    const unsigned int IMG_HEIGHT = 480;

    cv::Mat img(IMG_HEIGHT, IMG_WIDTH, CV_8UC3, cv::Scalar(0, 0, 0));

    writeMatrixToImg(FM, img, cv::Vec3b(0, 0, 255));

    for (int index = 0; index < all_indices.size(); index++) {
      Eigen::Vector2f v2 = SM_proj.col(index);
      img.at<cv::Vec3b>(cv::Point(round(v2(0)),round(v2(1)))) = cv::Vec3b(0, 255, 0);
    }

    cv::imwrite( "./images/modeltest.jpg", img);
    //##
    */
}

//-- BEGIN CERES STUFF
struct ReconstructionCostFunctor {
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  ReconstructionCostFunctor(const Eigen::Matrix<double, 3, 1>& v_face_avg_in,
                            const Eigen::Matrix<double, 3, 1>& v_scan_in,
                            const Eigen::Matrix<double, 3, NumberOfEigenvectors>& shapeBasisEigenRow_in,
                            const Eigen::Matrix<double, 3, NumberOfExpressions>& exprBasisEigenRow_in,
                            double weight_in = 1.0)
      : v_face_avg(v_face_avg_in), v_scan(v_scan_in),
        shapeBasisEigenRow(shapeBasisEigenRow_in),
        exprBasisEigenRow(exprBasisEigenRow_in), weight(weight_in) {}

  template <class T>
  bool operator()(T const* const salpha,/* T const* const sdelta,*/
                  /*T const* const sangles, T const* const st,*/
                  T* sresiduals) const {
    Eigen::Map<Eigen::Matrix<T, NumberOfEigenvectors, 1> const> const alpha(salpha);
    //Eigen::Map<Eigen::Matrix<T, NumberOfExpressions, 1> const> const delta(sdelta);
    Eigen::Map<Eigen::Matrix<T, 3, 1>> residuals(sresiduals);

    residuals = sqrt(T(weight)) * (v_scan - (v_face_avg + shapeBasisEigenRow * alpha /*+ exprBasisEigenRow * delta*/));

    return true;
  }

  const Eigen::Matrix<double, 3, 1> v_face_avg;
  const Eigen::Matrix<double, 3, 1> v_scan;

  // TODO: Again, these are copied so need another way.
  const Eigen::Matrix<double, 3, NumberOfEigenvectors> shapeBasisEigenRow;
  const Eigen::Matrix<double, 3, NumberOfExpressions> exprBasisEigenRow;

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

void runCeres(const MyMesh& avg_face_mesh, const MyMesh& scanned_mesh,
              const Eigen::MatrixXi& indices,
              const std::set<int>& matches,
              const Eigen::MatrixXf& shapeBasisEigen,
              const Eigen::MatrixXf& exprBasisEigen,
              const Eigen::VectorXf& shapeDevEigen,
              const Eigen::VectorXf& exprDevEigen,
              Eigen::VectorXd& alpha, Eigen::VectorXd& delta) {
    ceres::Problem problem;

    // Model transformation
    //Eigen::Matrix3f R;
    //Eigen::Vector3f t;

    for (MyMesh::VertexIter v_it = avg_face_mesh.vertices_begin();
         v_it != avg_face_mesh.vertices_end(); ++v_it)
    {
        // Get index of average face vertex
        int v_idx = v_it->idx();

        if ((v_idx % 20) != 0) {
            continue;
        }


        // Constants in optimization

        // Get vertex itself for average face mesh
        MyMesh::Point p3_face_avg = avg_face_mesh.point(*v_it);
        Eigen::Vector3d v3_face_avg(p3_face_avg[0], p3_face_avg[1], p3_face_avg[2]);

        // Get nearest neighbor vertex in scanned mesh
        MyMesh::VertexHandle v_scan_handle(indices(0, v_idx));
        MyMesh::Point p3_scan_nn = scanned_mesh.point(v_scan_handle);
        Eigen::Vector3d v3_scan_nn(p3_scan_nn[0], p3_scan_nn[1], p3_scan_nn[2]);

        // TODO:These are being copied but pointers should be used instead
        Eigen::Matrix<double, 3, NumberOfEigenvectors> shapeBasisEigenRows
                = shapeBasisEigen.block<3, NumberOfEigenvectors>(3 * v_idx, 0).cast<double>();
        Eigen::Matrix<double, 3, NumberOfExpressions> exprBasisEigenRows
                = exprBasisEigen.block<3, NumberOfExpressions>(3 * v_idx, 0).cast<double>();

        // Add more weight to sparse correspondences
        double weight = 1.0;
        if (matches.find(v_idx) != matches.end()) {
            weight = 100.0;
        }

        // For each entry in a vector v3
        problem.AddResidualBlock(
            // <dim of residual, dim of alpha, dim of delta>
            new ceres::AutoDiffCostFunction<ReconstructionCostFunctor, 3,
                    NumberOfEigenvectors/*, NumberOfExpressions*/>(
                new ReconstructionCostFunctor(v3_face_avg, v3_scan_nn,
                                              shapeBasisEigenRows,
                                              exprBasisEigenRows, weight)),
            nullptr, alpha.data()/*, delta.data()*/);

        //std::cout << "Adding residual: " << v_idx << "/" << avg_face_mesh.n_vertices() << "\n";
    }

    // Add regularization for alpha
    double lambda = 0.000005; // regularization param
    // Q: TODO: Is it a waste here to pass the whole alpha and delta? Would it
    // be better to make one large residual for alpha and delta?
    for (int i = 0; i < NumberOfEigenvectors; i++) {
        problem.AddResidualBlock(
            // <dim of residual, dim of alpha, dim of delta>
            new ceres::AutoDiffCostFunction<GeometryRegularizationCostFunctor, 1,
                    NumberOfEigenvectors>(
                new GeometryRegularizationCostFunctor(lambda, (double)shapeDevEigen(i), i)),
            nullptr, alpha.data());

    }
    for (int i = 0; i < NumberOfExpressions; i++) {
        problem.AddResidualBlock(
            // <dim of residual, dim of alpha, dim of delta>
            new ceres::AutoDiffCostFunction<GeometryRegularizationCostFunctor, 1,
                    NumberOfExpressions>(
                new GeometryRegularizationCostFunctor(lambda, (double)exprDevEigen(i), i)),
            nullptr, delta.data());
        // Q: TODO: Is it a waste here to pass the whole delta?
    }

    ceres::Solver::Options options;
    options.max_num_iterations = 1;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//DENSE_QR;
    options.minimizer_progress_to_stdout = true;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << std::endl;
    std::cout << summary.FullReport() << std::endl;
}
//-- END CERES STUFF

void LoadVector(const std::string &filename, float *res, unsigned int length)
{
	std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
	if (!in)
	{
		std::cout << "ERROR:\tCan not open file: " << filename << std::endl;
		while (1);
	}
	unsigned int numberOfEntries;
	in.read((char*)&numberOfEntries, sizeof(unsigned int));
	if (length == 0) length = numberOfEntries;
	in.read((char*)(res), length * sizeof(float));

	in.close();
}

void ProgressBar(const char* str, int num)
{
	return;
}

float* LoadEigenvectors(const std::string &filename, unsigned int components, unsigned int numberOfEigenvectors)
{
	float *res = new float[components*numberOfEigenvectors];

	for (unsigned int i = 0; i < numberOfEigenvectors; i++)
	{
		ProgressBar("Load Model Basis:", i / float(numberOfEigenvectors));
		std::stringstream ss;
		ss << filename << i << ".vec";

		LoadVector(ss.str().c_str(), &(res[components*i]), 0);
	}
	ProgressBar("Load Model Basis:", 1.0f);

	return res;
}

typedef Eigen::DiagonalMatrix<double, NumberOfEigenvectors + NumberOfExpressions> JacobiPrecondMatrix;


typedef struct float4 {
	float x;
	float y;
	float z;
	float w;
} float4;
constexpr const char* filenameBasisShape = "../models/ShapeBasis.matrix";
constexpr const char* filenameBasisExpression = "../models/ExpressionBasis.matrix";
constexpr const char* filenameStdDevShape = "../models/StandardDeviationShape.vec";

constexpr const char* filenameStdDevExpression = "../models/StandardDeviationExpression.vec";


int convert_basis(Eigen::MatrixXf& basis, int nVertices, int NumberOfEigenvectors, float4* basisCPU)
{
    for(int i = 0; i < NumberOfEigenvectors; i++) {
        for (int j = 0; j < nVertices; j++) {
            basis(j * 3 + 0, i) = basisCPU[i * nVertices + j].x;
            basis(j * 3 + 1, i) = basisCPU[i * nVertices + j].y;
            basis(j * 3 + 2, i) = basisCPU[i * nVertices + j].z;
        }
    }

    return 0;
}

int convert_deviations(Eigen::VectorXf& devs, int num_dims, float* devCPU)
{
    for(int i = 0; i < num_dims; i++) {
            devs(i, 0) = devCPU[i];
    }

    return 0;
}

int forward_pass(Eigen::MatrixXf& shape_basis, Eigen::MatrixXf& expr_basis, Eigen::VectorXd& alpha, Eigen::VectorXd& delta, Eigen::VectorXf& vertices_out)
{
    vertices_out = shape_basis * alpha.cast<float>() + expr_basis * delta.cast<float>();
	return 0;
}


int get_residuals_current(Eigen::VectorXf& residuals, Eigen::VectorXf& model_data, Eigen::VectorXi& residual_vert_ids)
{
	residuals(0) = model_data(50*3+0) - 5;
	residuals(1) = model_data(50*3+1) - 5;
	residuals(2) = model_data(50*3+2) - 10;
	residual_vert_ids(0) = 50;
	return 0;
}

int get_jacobian(Eigen::MatrixXf& jacobian, Eigen::MatrixXf& shape_basis, Eigen::MatrixXf& expr_basis, Eigen::VectorXf& residuals, int residuals_per_vert, int nResidualVerts, Eigen::VectorXi& residual_vert_ids)
{
	
	for (int cc1 = 0; cc1 < 1; cc1++)
	{
		for (int cc2 = 0; cc2 < 160; cc2++)
		{
			jacobian(cc1*residuals_per_vert + 0, cc2) = shape_basis(residual_vert_ids(cc1)* 3+0,cc2) ;
			jacobian(cc1*residuals_per_vert + 1, cc2) = shape_basis(residual_vert_ids(cc1) * 3 + 1, cc2);
			jacobian(cc1*residuals_per_vert + 2, cc2) = shape_basis(residual_vert_ids(cc1) * 3 + 2, cc2);
		}
		for (int cc2 = 0; cc2 < 76; cc2++)
		{
			jacobian(cc1*residuals_per_vert + 0, cc2+160) = expr_basis(residual_vert_ids(cc1)*3+0, cc2);
			jacobian(cc1*residuals_per_vert + 1, cc2 + 160) = expr_basis(residual_vert_ids(cc1) * 3 + 1, cc2);
			jacobian(cc1*residuals_per_vert + 2, cc2 + 160) = expr_basis(residual_vert_ids(cc1) * 3 + 2, cc2);

		}

		//jacobian = 2 * residuals.asDiagonal() * jacobian;

	}
	return 0;
}

int compute_jacobi_precond(JacobiPrecondMatrix& Jacobi_precond_local, Eigen::MatrixXd& J_t, Eigen::MatrixXd& J)
{
	Jacobi_precond_local = ((J_t * J).diagonal()).asDiagonal();
	Jacobi_precond_local=Jacobi_precond_local.inverse();
	Eigen::VectorXd* test_vec = new Eigen::VectorXd(236);
	*test_vec = Eigen::VectorXd::Constant(236, 1);

	Jacobi_precond_local.applyThisOnTheLeft(*test_vec);

	double vec_norm=test_vec->norm();
	
	delete test_vec;

	return 0;
}

int conjugateGradientSolver(Eigen::MatrixXf& J_t, Eigen::MatrixXf& J,  Eigen::VectorXf& B, Eigen::VectorXf& X)
{
	double TOLERANCE = 1.0e-8;

	constexpr double NEARZERO = 1.0e-10;

	

	int n = J_t.rows();

	Eigen::VectorXd * X_local = new Eigen::VectorXd(X.rows());

	Eigen::MatrixXd * J_t_local = new Eigen::MatrixXd(J_t.rows(), J_t.cols());
	Eigen::MatrixXd * J_local = new Eigen::MatrixXd(J.rows(), J.cols());

	JacobiPrecondMatrix Jacobi_precond_local;

	

	*J_t_local = J_t.cast<double>();
	*J_local = J.cast<double>();

	compute_jacobi_precond(Jacobi_precond_local, *J_t_local, *J_local);

	*X_local = Eigen::VectorXd::Zero(X.rows());

	Eigen::VectorXd * X_best = new Eigen::VectorXd(X.rows());

	Eigen::VectorXd * R = new Eigen::VectorXd(B.rows());
	Eigen::VectorXd * Rold = new Eigen::VectorXd(B.rows());
	Eigen::VectorXd * P = new Eigen::VectorXd(B.rows());
	Eigen::VectorXd * AP = new Eigen::VectorXd(J_t.rows());

	*R = B.cast<double>();

	Jacobi_precond_local.applyThisOnTheLeft(*R);

	*P = *R;
	int k = 0;

	double bestnorm = -1;

	
	*X_best = *X_local;

	while (k < n)
	{
		*Rold = *R;                                         // Store previous residual
		*AP =  ( (*J_t_local * (*J_local * *P)));
		Jacobi_precond_local.applyThisOnTheLeft(*AP);

		double alpha = R->squaredNorm() / std::max((const double)(P->dot(*AP)), NEARZERO);
		*X_local = *X_local + (alpha * *P);            // Next estimate of solution
		*R = *R -(alpha * *AP);          // Residual 

		double cur_res_norm = R->norm();

		if (cur_res_norm < TOLERANCE)
		{
			*X_best= *X_local;

			break;             // Convergence test // Use norm?
		}

		if ((bestnorm == -1) || (cur_res_norm < bestnorm))
		{
			bestnorm = cur_res_norm;
			*X_best = *X_local;
		}

		

		double beta = R->squaredNorm() / std::max((double)(Rold->squaredNorm()), NEARZERO);
		*P = *R + (beta * *P);             // Next gradient
		k++;
	}
	X = X_best->cast<float>();
	delete R, Rold, P, AP, X_local, J_local, J_t_local;

	return 0;
}


int update_params(Eigen::VectorXf& alpha, Eigen::VectorXf& delta, Eigen::VectorXf& delta_P)
{
	for (int cc1 = 0; cc1 < alpha.rows(); cc1++)
	{
		alpha(cc1) += delta_P(cc1);
	}
	for (int cc1 = 0; cc1 < delta.rows(); cc1++)
	{
		delta(cc1) += delta_P(alpha.rows()+cc1);
	}
	return 0;
}

int cg_solver_helper(Eigen::MatrixXf& A, Eigen::VectorXf& B, Eigen::VectorXf& X);

int main()
{
    FaceModel face_model;
    MyMesh scanned_mesh;

    loadScannedMesh(scanned_mesh);

    const int K = 1;
    Eigen::MatrixXi indices;


	auto shapeBasisCPU = new float4[nVertices * NumberOfEigenvectors];
	auto expressionBasisCPU = new float4[nVertices * NumberOfExpressions];
	LoadVector(filenameBasisShape, (float*)shapeBasisCPU, 4 * nVertices * NumberOfEigenvectors);
	LoadVector(filenameBasisExpression, (float*)expressionBasisCPU, 4 * nVertices * NumberOfExpressions);

    auto shapeDevCPU = new float[NumberOfEigenvectors];
    auto expressionDevCPU = new float[NumberOfExpressions];
    LoadVector(filenameStdDevShape, shapeDevCPU, NumberOfEigenvectors);
    LoadVector(filenameStdDevExpression, expressionDevCPU, NumberOfExpressions);

    Eigen::MatrixXf * shapeBasisEigen = new Eigen::MatrixXf(3 * nVertices, NumberOfEigenvectors);
	Eigen::MatrixXf * exprBasisEigen = new Eigen::MatrixXf(3 * nVertices, NumberOfExpressions);
	
    Eigen::VectorXf * shapeDevEigen = new Eigen::VectorXf(NumberOfEigenvectors);
    Eigen::VectorXf * exprDevEigen = new Eigen::VectorXf(NumberOfExpressions);

	std::cout << "converting the basis: " << std::endl;

	convert_basis(*shapeBasisEigen, nVertices, NumberOfEigenvectors, shapeBasisCPU);
    convert_basis(*exprBasisEigen, nVertices, NumberOfExpressions, expressionBasisCPU);//shapeBasisCPU); // Q: Should be expressionBasisCPU?

    convert_deviations(*shapeDevEigen, NumberOfEigenvectors, shapeDevCPU);
    convert_deviations(*exprDevEigen, NumberOfExpressions, expressionDevCPU);

    //Eigen::VectorXf * alpha = new Eigen::VectorXf(NumberOfEigenvectors);
    //Eigen::VectorXf * delta = new Eigen::VectorXf(NumberOfExpressions);

    //*alpha = Eigen::VectorXf::Random(NumberOfEigenvectors) / 10;
    //*delta= Eigen::VectorXf::Random(NumberOfExpressions) / 10;
    Eigen::VectorXd alpha = Eigen::VectorXd::Zero(NumberOfEigenvectors);
    Eigen::VectorXd delta = Eigen::VectorXd::Zero(NumberOfExpressions);

	Eigen::VectorXf * vertices_out = new Eigen::VectorXf(3 * nVertices);
    //std::cout << "forward pass: " << std::endl;

    std::cout << "Optimization starting: " << std::endl;
    for (int i = 0; i < 3; i++) {
        knn_test(face_model, scanned_mesh, K, indices);

        std::map<int, int> match_indices;
        match_indices[2000] = 2000;
        match_indices[10000] = 10000;
        match_indices[25000] = 25000;
        match_indices[30000] = 30000;
        match_indices[50000] = 50000;

        std::set<int> matches;
        for (auto entry : match_indices) {
            matches.insert(entry.first);
            indices(0, entry.first) = entry.second;
        }

        std::cout << "Running ceres: " << std::endl;
        runCeres(face_model.m_synth_mesh, scanned_mesh, indices, matches,
                 *shapeBasisEigen, *exprBasisEigen, *shapeDevEigen, *exprDevEigen,
                 alpha, delta);

        std::cout << "Ceres finished: " << std::endl;

        forward_pass(*shapeBasisEigen, *exprBasisEigen, alpha, delta, *vertices_out);

        face_model.synthesizeModel(*vertices_out);
    }

    std::cout << "Optimization finished: " << std::endl;

    std::cout << alpha << std::endl;
/*

	constexpr int nResidualVerts = 1;
	constexpr int nResidualsPerVert = 3;

	Eigen::VectorXf * residuals = new Eigen::VectorXf(nResidualVerts*nResidualsPerVert);

	Eigen::VectorXi * resudial_vertex_ids = new Eigen::VectorXi(nResidualVerts);

	std::cout << "getting residuals and Jacobian: " << std::endl;

	Eigen::MatrixXf * jacobian = new Eigen::MatrixXf(nResidualVerts*nResidualsPerVert, NumberOfEigenvectors + NumberOfExpressions);

	Eigen::VectorXf * delta_P = new Eigen::VectorXf(NumberOfEigenvectors+ NumberOfExpressions);

	Eigen::VectorXf * resid = new Eigen::VectorXf(NumberOfEigenvectors + NumberOfExpressions);

	Eigen::VectorXf * B = new Eigen::VectorXf(NumberOfEigenvectors + NumberOfExpressions);

	get_residuals_current(*residuals, *vertices_out, *resudial_vertex_ids);

	get_jacobian(*jacobian, *shapeBasisEigen, *exprBasisEigen,*residuals, nResidualsPerVert, nResidualVerts, *resudial_vertex_ids);

	for(int cc1=0; cc1<20; cc1++)
	{


		// Newton method
		
		forward_pass(*shapeBasisEigen, *exprBasisEigen, *alpha, *delta, *vertices_out);

		get_residuals_current(*residuals, *vertices_out, *resudial_vertex_ids);

		std::cout << "loss: " << residuals->squaredNorm() << std::endl;

		get_jacobian(*jacobian, *shapeBasisEigen, *exprBasisEigen, *residuals, nResidualsPerVert, nResidualVerts, *resudial_vertex_ids);


		Eigen::MatrixXf * jacobianT = new Eigen::MatrixXf(NumberOfEigenvectors + NumberOfExpressions, nResidualVerts*nResidualsPerVert);

		Eigen::MatrixXf * jacobianTjacobian = new Eigen::MatrixXf(NumberOfEigenvectors + NumberOfExpressions, NumberOfEigenvectors + NumberOfExpressions);

		*jacobianT = jacobian->transpose();
		*jacobianTjacobian = *jacobianT * *jacobian; // inefficient

		*B = -1 * (*jacobianT * *residuals);

		conjugateGradientSolver(*jacobianT, *jacobian, *B, *delta_P);

		*resid = *jacobian * *delta_P + *residuals;

		update_params(*alpha, *delta, *delta_P);

	}
*/
    std::cout << "Writing synthesized model to file\n";
    face_model.writeSynthesizedModel();

	std::cout << "Hello World!\n"; 
}


// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
