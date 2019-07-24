// ConsoleApplication1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>

#include <string>
#include <Eigen/Dense>
#include <vector>

//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

#include "facemodel.h"
#include "rgbdscan.h"
#include "facesolver.h"
#include"config.h"


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

/*
void writeMatrixToImg(const Eigen::MatrixXf& matrix, cv::Mat& img,
                 cv::Vec3b vcolor) {
    //if (0 <= v2(0) && v2(0) < IMG_WIDTH && 0 <= v2(1) && v2(1) < IMG_HEIGHT) {
    for (int c = 0; c < matrix.cols(); c++) {
        Eigen::Vector2f v2 = matrix.col(c);
      img.at<cv::Vec3b>(cv::Point(round(v2(0)),round(v2(1)))) = vcolor;
    }
    //}
}
*/
typedef Eigen::DiagonalMatrix<double, NumberOfEigenvectors + NumberOfExpressions> JacobiPrecondMatrix;


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

int optimize() {
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
}

Sophus::SE3d createTransformTest() {
    Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();

    Eigen::Matrix3d m;
    m = Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitX())
        * Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(-M_PI/3, Eigen::Vector3d::UnitZ());

    transform.block<3, 3>(0, 0) = m;
    transform.block<3, 1>(0, 3) = Eigen::Vector3d(0.5, 0.1, 0.1);

    Sophus::SE3d T_xy(transform);

    return T_xy;
}


int main()
{
    Eigen::VectorXd alpha = Eigen::VectorXd::Zero(NumberOfEigenvectors);
    Eigen::VectorXd delta = Eigen::VectorXd::Zero(NumberOfExpressions);

    Sophus::SE3d T_xy;

    FaceModel face_model(MODEL_PATH);
    RGBDScan face_scan(FILENAME_SCANNED_MESH, FILENAME_SCANNED_LANDMARKS);
    FaceSolver face_solver(0.0005, 0.004, 0.004, 6, 1, true);// percent = 0.1 works pretty okay too

    face_solver.solve(face_model, face_scan, alpha, delta, T_xy);

    std::cout << "Writing synthesized model to file\n";
    face_model.writeSynthesizedModel(alpha, delta, T_xy);

    // Print alpha (geometry parameter we solved for)
    std::cout << alpha << std::endl;

    std::cout << T_xy.matrix() << std::endl;
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
