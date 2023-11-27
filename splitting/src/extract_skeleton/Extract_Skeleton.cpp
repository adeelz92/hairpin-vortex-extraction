#define CGAL_EIGEN3_ENABLED
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                        Kernel;
typedef Kernel::Point_3                                       Point;
typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;
typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;
//only needed for the display of the skeleton as maximal polylines
struct Display_polylines {
	const Skeleton& skeleton;
	std::ofstream& out;
	int polyline_size;
	std::stringstream sstr;
	Display_polylines(const Skeleton& skeleton, std::ofstream& out)
		: skeleton(skeleton), out(out)
	{}
	void start_new_polyline() {
		polyline_size = 0;
		sstr.str("");
		sstr.clear();
	}
	void add_node(Skeleton_vertex v) {
		++polyline_size;
		sstr << skeleton[v].point << "\n";
	}
	void end_polyline()
	{
		out << polyline_size << "\n" << sstr.str();
	}
};
// This example extracts a medially centered skeleton from a given mesh.
int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		std::cerr << "Please specify the output filename." << std::endl;
		return EXIT_FAILURE;
	}

	const std::string filename = argv[1];
	const std::string output_filename = argv[2];

	Polyhedron tmesh;
	if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, tmesh))
	{
		std::cerr << "Invalid input." << std::endl;
		return EXIT_FAILURE;
	}

	Skeleton skeleton;
	Skeletonization mcs(tmesh);
	// CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);
	mcs.set_quality_speed_tradeoff(0.5);
	mcs.set_medially_centered_speed_tradeoff(0.08);
	mcs.contract_until_convergence();
	mcs.convert_to_skeleton(skeleton);
	// std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
	// std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";
	// Output all the edges of the skeleton.
	std::ofstream output(output_filename);
	Display_polylines display(skeleton, output);
	CGAL::split_graph_into_polylines(skeleton, display);
	output.close();

	return EXIT_SUCCESS;
}
