#include <pcl/PCLPointCloud2.h>
#include <pcl/io/ply_io.h>
#include <pcl/console/print.h>
#include <pcl/console/parse.h>
#include <pcl/console/time.h>
#include <pcl/surface/poisson.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/integral_image_normal.h>

using namespace pcl;
using namespace pcl::io;
using namespace pcl::console;

int default_depth = 8;
int default_solver_divide = 8;
int default_iso_divide = 8;
float default_point_weight = 4.0f;

bool
loadCloud (const std::string &filename, pcl::PCLPointCloud2 &cloud)
{
  TicToc tt;
  print_highlight ("Loading "); print_value ("%s ", filename.c_str ());

  tt.tic ();
  if (pcl::io::loadPLYFile (filename, cloud) < 0)
    return (false);
  print_info ("[done, "); print_value ("%g", tt.toc ()); print_info (" ms : "); print_value ("%d", cloud.width * cloud.height); print_info (" points]\n");
  print_info ("Available dimensions: "); print_value ("%s\n", pcl::getFieldsList (cloud).c_str ());

  return (true);
}

void
compute_normals (const pcl::PCLPointCloud2::ConstPtr &input, pcl::PCLPointCloud2 &output,
         int k, double radius)
{
  // Convert data to PointCloud<T>
  PointCloud<PointXYZ>::Ptr xyz (new PointCloud<PointXYZ>);
  fromPCLPointCloud2 (*input, *xyz);

  TicToc tt;
  tt.tic ();
 
  PointCloud<Normal> normals;

  // Try our luck with organized integral image based normal estimation
  if (xyz->isOrganized ())
  {
    IntegralImageNormalEstimation<PointXYZ, Normal> ne;
    ne.setInputCloud (xyz);
    ne.setNormalEstimationMethod (IntegralImageNormalEstimation<PointXYZ, Normal>::COVARIANCE_MATRIX);
    ne.setNormalSmoothingSize (float (radius));
    ne.setDepthDependentSmoothing (true);
    ne.compute (normals);
  }
  else
  {
    NormalEstimation<PointXYZ, Normal> ne;
    ne.setInputCloud (xyz);
    ne.setSearchMethod (search::KdTree<PointXYZ>::Ptr (new search::KdTree<PointXYZ>));
    ne.setKSearch (k);
    ne.setRadiusSearch (radius);
    ne.compute (normals);
  }

  print_highlight ("Computed normals in "); print_value ("%g", tt.toc ()); print_info (" ms for "); print_value ("%d", normals.width * normals.height); print_info (" points.\n");

  // Convert data back
  pcl::PCLPointCloud2 output_normals;
  toPCLPointCloud2 (normals, output_normals);
  concatenateFields (*input, output_normals, output);
}

void
compute (const pcl::PCLPointCloud2::ConstPtr &input, PolygonMesh &output,
         int depth, int solver_divide, int iso_divide, float point_weight)
{
  PointCloud<PointNormal>::Ptr xyz_cloud (new pcl::PointCloud<PointNormal> ());
  fromPCLPointCloud2 (*input, *xyz_cloud);

  print_info ("Using parameters: depth %d, solverDivide %d, isoDivide %d\n", depth, solver_divide, iso_divide);

	Poisson<PointNormal> poisson;
	poisson.setDepth (depth);
	poisson.setSolverDivide (solver_divide);
	poisson.setIsoDivide (iso_divide);
  poisson.setPointWeight (point_weight);
  poisson.setInputCloud (xyz_cloud);

  TicToc tt;
  tt.tic ();
  print_highlight ("Computing ...");
  poisson.reconstruct (output);

  print_info ("[Done, "); print_value ("%g", tt.toc ()); print_info (" ms]\n");
}

void
saveCloud (const std::string &filename, const PolygonMesh &output)
{
  TicToc tt;
  tt.tic ();

  print_highlight ("Saving "); print_value ("%s ", filename.c_str ());
  pcl::io::savePLYFile (filename, output);

  print_info ("[done, "); print_value ("%g", tt.toc ()); print_info (" ms]\n");
}

/* ---[ */
int
main ()
{
  // Command line parsing
  int depth = default_depth;

  int solver_divide = default_solver_divide;

  int iso_divide = default_iso_divide;

  float point_weight = default_point_weight;

  // Load the first file
  pcl::PCLPointCloud2::Ptr cloud (new pcl::PCLPointCloud2);
  if (!loadCloud ("/Users/ardakeskiner/Downloads/2019_07_19_16_37_07/FaceRaw/Face_Raw_000000.ply", *cloud))
    return (-1);

  // Apply the Poisson surface reconstruction algorithm
  PolygonMesh output;
  pcl::PCLPointCloud2 intermediate;
  compute_normals(cloud, intermediate, 5, 0.0);

  pcl::PCLPointCloud2::Ptr intermediate_2 (new pcl::PCLPointCloud2(intermediate));
  compute (intermediate_2, output, depth, solver_divide, iso_divide, point_weight);

  // Save into the second file
  saveCloud ("asd.ply", output);
}