// ConsoleApplication1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <Eigen/Dense>
#include <vector>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

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
    mesh.set_point( *v_it, mesh.point(*v_it) + 0.05 * mesh.normal(*v_it) );
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


constexpr int NumberOfEigenvectors = 160;
constexpr int NumberOfExpressions = 76;

typedef Eigen::DiagonalMatrix<double, NumberOfEigenvectors + NumberOfExpressions> JacobiPrecondMatrix;

constexpr int nVertices = 53490;
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
	for (int cc1 = 0; cc1 < nVertices; cc1++)
	{
			for (int cc3 = 0; cc3 < NumberOfEigenvectors; cc3++)
			{
				basis(cc1 * 3 + 0,cc3) = basisCPU[cc1*NumberOfEigenvectors + cc3].x;
				basis(cc1 * 3 + 1,cc3) = basisCPU[cc1*NumberOfEigenvectors + cc3].y;
				basis(cc1 * 3 + 2,cc3) = basisCPU[cc1*NumberOfEigenvectors + cc3].z;
			}
	}
	return 0;
}
int forward_pass(Eigen::MatrixXf& shape_basis, Eigen::MatrixXf& expr_basis, Eigen::VectorXf& alpha, Eigen::VectorXf& delta, Eigen::VectorXf& vertices_out)
{
	vertices_out = shape_basis * alpha + expr_basis * delta;
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
	auto shapeBasisCPU = new float4[nVertices * NumberOfEigenvectors];
	auto expressionBasisCPU = new float4[nVertices * NumberOfExpressions];
	LoadVector(filenameBasisShape, (float*)shapeBasisCPU, 4 * nVertices * NumberOfEigenvectors);
	LoadVector(filenameBasisExpression, (float*)expressionBasisCPU, 4 * nVertices * NumberOfExpressions);

	auto shapeDevCPU = new float[NumberOfEigenvectors];
	auto expressionDevCPU = new float[NumberOfExpressions];

	Eigen::MatrixXf * shapeBasisEigen = new Eigen::MatrixXf(3 * nVertices, NumberOfEigenvectors);
	Eigen::MatrixXf * exprBasisEigen = new Eigen::MatrixXf(3 * nVertices, NumberOfExpressions);
	
	std::cout << "converting the basis: " << std::endl;

	convert_basis(*shapeBasisEigen, nVertices, NumberOfEigenvectors, shapeBasisCPU);
	convert_basis(*exprBasisEigen, nVertices, NumberOfExpressions, shapeBasisCPU);

	Eigen::VectorXf * alpha = new Eigen::VectorXf(NumberOfEigenvectors);
	Eigen::VectorXf * delta = new Eigen::VectorXf(NumberOfExpressions);

	*alpha = Eigen::VectorXf::Random(NumberOfEigenvectors)/10;
	*delta= Eigen::VectorXf::Random(NumberOfExpressions)/10;

	Eigen::VectorXf * vertices_out = new Eigen::VectorXf(3 * nVertices);
	std::cout << "forward pass: " << std::endl;

	forward_pass(*shapeBasisEigen, *exprBasisEigen, *alpha, *delta, *vertices_out);
	
 	LoadVector(filenameStdDevShape, shapeDevCPU, NumberOfEigenvectors);
	LoadVector(filenameStdDevExpression, expressionDevCPU, NumberOfExpressions);

	

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
