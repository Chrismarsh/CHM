// 	Copyright (C) 2011  Chris Marsh
// 
// 	This program is free software: you can redistribute it and/or modify
// 	it under the terms of the GNU General Public License as published by
// 	the Free Software Foundation, either version 3 of the License, or
// 	(at your option) any later version.
// 
// 	This program is distributed in the hope that it will be useful,
// 	but WITHOUT ANY WARRANTY; without even the implied warranty of
// 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// 	GNU General Public License for more details.
// 
// 	You should have received a copy of the GNU General Public License
// 	along with this program.  If not, see <http://www.gnu.org/licenses/>.


#pragma once

//for valgrind, remove
#define CGAL_DISABLE_ROUNDING_MATH_CHECK

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <set>
#include <stack>
#include <fstream>
#include <utility>
#define ARMA_DONT_USE_CXX11 //intel on linux breaks otherwise
#include <armadillo>

#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/ptr_container/ptr_map.hpp>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Triangulation_2.h>


#include <tbb/concurrent_vector.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;

//required for the spatial searching
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <CGAL/Euclidean_distance.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>

//http://www.paraview.org/Bug/print_bug_page.php?bug_id=14164
//http://review.source.kitware.com/#/c/11956/
//until 6.0.1 comes out
//#define VTK_HAS_STD_ISNAN
//#define VTK_HAS_STD_ISINF

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h> //pointdata
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkProbeFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPointData.h>
#include <vtkPlaneSource.h>
#include <vtkGeometryFilter.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkElevationFilter.h>
#include <vtkCurvatures.h>
#include <vtkXMLPolyDataWriter.h>
#ifdef NOMATLAB
#include "libmaw.h"
#endif

#include "vertex.hpp"
//#include "face.hpp"
#include "timeseries.hpp"
#include "math/distance.hpp"

/**
* \struct face_info
* A way of embedding arbirtrary data into the face. This is how modules should store their data.
*/
struct face_info
{

    virtual ~face_info()
    {
    };
};

//fwd decl
class segmented_AABB;
class triangulation;


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Triangle_3 Triangle_3;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef CGAL::Projection_traits_xy_3<K> Gt; //allows for using 2D algorithms on the 3D points

typedef ex_vertex<Gt> Vb; //custom vertex class
/**
* \class face
* \brief Defines the triangle face
*/

template < class Gt, class Fb = CGAL::Constrained_triangulation_face_base_2<Gt> >
class face : public Fb
{
    friend class triangulation;
public:

    /**
     * Multiuse flag for algorithms to use.
     * Does not make any gaurantee about state.
     */
    bool coloured;
//    face_info* info;

    /**
    * \typedef Vertex_handle Handle to a vertex
    */
    typedef typename Fb::Vertex_handle Vertex_handle;

    /**
    * \typedef Face_handle Handle to a face
    */
    typedef typename Fb::Face_handle Face_handle;


    template < typename TDS2 >
    struct Rebind_TDS
    {
        typedef typename Fb::template Rebind_TDS<TDS2>::Other Fb2;
        typedef face<Gt, Fb2> Other;
    };

    face();
    face(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2);
    face(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Face_handle n0, Face_handle n1, Face_handle n2);
    face(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Face_handle n0, Face_handle n1, Face_handle n2, bool c0, bool c1, bool c2);

    ~face();

    /**
    * Aspect of the face. North = 0, CW . Calculated on first use, subsequent usages will not recalculate
    * \return Face aspect [rad]
    */
    double aspect();

    /**
    * Slope of the face. Calculated on first use, subsequent usages will not recalculate
    * \return slope [rad]
    */
    double slope();

    /**
    * Normalized face normal. Calculated on first use, subsequent usages will not recalculate
    */
    Vector_3 normal();

    /**
    * Center of the face as defined by a centroid. Calculated on first use, subsequent usages will not recalculate
    */
    Point_3 center();

    /**
     * Exactly the same as find_closest_face in triangulation but uses the current face's center
     * @param azmimuth
     * @param distance
     * @return
     */
    const Face_handle find_closest_face(double azmimuth, double distance);

    /**
    * Checks if a point x,y is within the face
    * \return true if this face contains the point x,y
    */
    bool contains(double x, double y);

    /**
    * Checks if a point p is within the face
    * \return true if this face contains the point p
    */
    bool contains(Point_3 p);

    /**
    * Checks if the given face intersects this face in 2D. The algorithm checks if the specified face contains any of the 3 vertexs or center of this face.
    * \param fh Face handle to the other face
    * \return true if this face and fh intersects
    */
    bool intersects(Face_handle fh);

    /**
    * Sets the face data, for the given variable, at the current timestep
    * \param variable Variable to assign the data to
    * \param data Data
    */
    void set_face_data(const std::string& variable, double data);

    /**
     * Sets the vector for the given variable.
     * Does not support timeseries output.
     * There are no checks on inter-module vector guarantees like normal variables
     * Proceed with caution
     * @param variable
     * @param v
     */
    void set_face_vector(const std::string& variable, Vector_3 v);

    /**
    * Returns the face data
    * \param variable  variable being queried
    * \return  this timesteps value for the given variable
    */
    double face_data(const std::string& variable);

    /**
     * Returns the face vector for a specified variable
     * @param variable
     * @return
     */
    Vector_3 face_vector(const std::string& variable);


    /**
     * Removes the face data associated with a given ID. Throws if ID not found.
     * Assumes memory was allocated in a module. This just makes it easier to cleanup.
     * \param ID Module ID
     * \return none
     */
    void remove_face_data(std::string ID);

    /**
    * Initializes  this faces timeseires with the given variables, for the given datetime series with the given size
    * \param variables Names of the variables to add
    * \param datetime Vector of boost::ptimes for the entire duration of the timesries
    */
    void init_time_series(std::set<std::string> variables, timeseries::date_vec datetime);

    /**
    * Obtains the timeseries associated with the given variable
    * \param ID variable
    */
    timeseries::variable_vec face_time_series(std::string ID);

    /**
     * Returns true if the variable is present
     */
    bool has(const std::string& variable);

    /**
     * Returns a list of variables in this face's timeseries
     * \return Vector of variable names
     */
    std::vector<std::string> variables();

    /**
    * Returns the underlying timeseries object
    * \return Pointer to the underlying timeseries
    */
    boost::shared_ptr<timeseries> get_underlying_timeseries();

    /**
    * Returns the iterator of the current timestep.
    */
    timeseries::iterator now();

    /**
     * Increments the face to the next timestep
     */
    void next();

    /**
     * Resets the face to the start of the timeseries
     */
    void reset_to_begining();

    /**
     * Returns the xcoordinates of the center of the face
     * @return x UTM coordinate
     */
    double get_x();

    /**
     * Returns the ycoordinates of the center of the face
     * @return y UTM coordinate
     */
    double get_y();

    /**
     * Returns the zcoordinates of the center of the face
     * @return  elevation
     */
    double get_z();

    /**
     * Get triangle area (m)
     */
    double get_area();
    /**
     * Saves this face's timeseries to a file
     * @param fname specified file name
     */
    void to_file(std::string fname);

    template<typename T>
    T*get_module_data(std::string module);

    void set_module_data(std::string module, face_info *fi);

    template<typename T>
    T*make_module_data(std::string module);

    std::string _debug_name; //for debugging to find the elem that we want
    int _debug_ID; //also for debugging. ID == the position in the output order, starting at 0
    size_t cell_id;
    void set_parameter(std::string key,double value);
    double get_parameter(std::string key);
    std::vector<std::string>  parameters();
    bool has_parameter(std::string key);

    void set_initial_condition(std::string key,double value);
    double get_initial_condition(std::string key);

    /**
     * Returns a vector of names of all intial conditions
     * @return
     */
    std::vector<std::string> initial_conditions();

    /**
     * Returns a vector of all the names of all the vector data
     * @return
     */
    std::vector<std::string> vectors() ;

    bool has_initial_condition(std::string key);

    bool _is_geographic;
private:


    double _slope;
    double _azimuth;
    double _area;
    double _x;
    double _y;
    double _z;

    //hold a pointer *back* to the triangulation. This let's use query triangles at distance X, etc
    //that allows for using data::parallel modules w/o having to use domain parallel.
    //const so we can't modify the domain via this as thar be dragons
    triangulation* _domain;

    boost::shared_ptr<Point_3> _center;
    boost::shared_ptr<Vector_3> _normal;

    //boost::ptr_map<std::string,face_info> _module_face_data;
    std::map<std::string,face_info* > _module_face_data;
    std::map<std::string,double> _parameters;
    std::map<std::string,double> _initial_conditions;
    std::map<std::string,Vector_3> _module_face_vectors; //holds vector components, currently no checks on anything. Proceed with caution.

    boost::shared_ptr<timeseries> _data;
    timeseries::iterator _itr;

};

typedef face<Gt> Fb; //custom face class
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Delaunay;

typedef Delaunay::Face_handle mesh_elem;
typedef boost::shared_ptr<tbb::concurrent_vector<double>  > vector;

//search tree typedefs
typedef K::Point_2 Point_2;
typedef boost::tuple<Point_2, Delaunay::Face_handle > Point_and_face;
typedef CGAL::Search_traits_2<K>                       Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_face,
		CGAL::Nth_of_tuple_property_map<0, Point_and_face>,
		Traits_base>                                              Traits;


typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
typedef K_neighbor_search::Tree                             Tree;
typedef K_neighbor_search::Distance                         Distance;




/**
*
*/
class triangulation
: public Delaunay
{
public:
    triangulation();
#ifdef MATLAB
/**
* If the Matlab engine integration is enabled, a pointer to the engine instance must be passed in.
* \param engine Pointer to a valid matlab engine instant
*/
    triangulation(boost::shared_ptr<maw::matlab_engine> engine);
#endif
    ~triangulation();

    /**
    * Creates an axis alligned bounding box that is segmented into rows x cols
    * \param rows number of rows in the AABB
    * \param cols number of cols in the AABB
    */
    boost::shared_ptr<segmented_AABB> AABB(size_t rows, size_t cols);

    /**
    * Loads a mesh from file. Should by x y z values with no header, space delimited.
    * \param file Fully qualified path to a file.
    */
	std::set<std::string> from_json(pt::ptree& mesh);

	/**
	 * Serializes a mesh attribute to file so it can be read into the model.
	 * \param output_path Full qualified path to a file to output to.
	 * \param variable The variable from a module to be output.
	 */
	void serialize_parameter(std::string output_path, std::string parameter);
    /**
    * Initializes a triangulation from the given x,y,z vectors
    * \param x X values
    * \param y Y values
    * \param z Z values
    */
//    void init(vector x, vector y, vector z);

    /**
    * Return the number of faces in the triangluation
    * \return Number of triangle faces
    */
    size_t size_faces();

    /**
    * Return the number of verticies in the triangulation
    * \return Number of vertices
    */
    size_t size_vertex();

    /**
     * Locates the closest triangle (based on centers) that the point x y lies on in 2 dimensions.
     * This will always return a triangle, even if the query point is outside of the domain.
     * The triangle found does not nessessairly contain the point specified if another triangle's center is closer.
     * @param x x coordinate of query point
     * @param y y coordinate of query point
     * @return
     */
    mesh_elem find_closest_face(double x, double y) const;

    /**
     * Locates the closest triangle (based on centers) that the point x y lies on in 2 dimensions.
     * This will always return a triangle, even if the query point is outside of the domain
     * The triangle found does not nessessairly contain the point specified if another triangle's center is closer.
     * @param query Query point
     * @return
     */
	mesh_elem find_closest_face(Point_2 query) const;

    /**
     * Lcoates the triangle that contains the query point. Guaranteed that if a triangle is found, the point lies inside the triangle.
     *
     * @param x
     * @param y
     * @return triangle if found, otherwise nullptr
     */
    mesh_elem locate_face(double x, double y);
    mesh_elem locate_face(Point_2 query);

    /**
    * Returns the finite face at index i. A given index will always return the same face.
    * \param i Index
    * \return A face handle to the ith face
    */
    Delaunay::Face_handle face(size_t i);

    /**
    * Returns the finite vertex at index i. A given index will always return the same vertex.
    * \param i Index
    * \return A vertex handle to the ith vertex
    */
    Delaunay::Vertex_handle vertex(size_t i);
#ifdef NOMATLAB
    /**
    * If Matlab integration is enabled, plots the given variable at the current timestep.
    * \param ID variable
    */
    void plot(std::string ID);

    /**
    * If Matlab integration is enabled, plots the given variable over time for the given triangle that contains the point x,y
    * \param x x coord
    * \param y y coord
    * \param ID variable to plot
    */

    void plot_time_series(double x, double y, std::string ID);
#endif

    /**
    * Saves the timeseries as a csv for the triangle that contains the point x,y to file
    * \param x x coord
    * \param y y coord
    * \param fname File same
    */
    void timeseries_to_file(double x, double y, std::string fname);

    /**
    * Saves the timeseries as a csv for the given triangle to file
    * \param m mesh element
    * \param fname File same
    */
    void timeseries_to_file(mesh_elem m, std::string fname);

	/**
	 * If output to the mesh vtk/vtu format is required, this will be allocate the vtk data structure.
	 */
	void init_vtkUnstructured_Grid(std::vector<std::string> output_variables);


	/**
	 * Updates the internal vtk structure with this timesteps data.
	 * Must be called prior to calling the write_vt* functions.
	 * The write_vt* functions could call this, however it makes them not threadsafe.
	 * If output_variables is empty, it will write all variables out
	 * @param output_variables Selected variables to write out.
	 */
    void update_vtk_data(std::vector<std::string> output_variables);

    /**
    * Saves the mesh with this timesteps values to a vtu file for visualization in Paraview
    */
	void write_vtu(std::string fname);

    void write_vtp(std::string file_name);

	/**
	 * Returns true if this is a geogrphic mesh
	 * @return
	 */
	bool is_geographic();

    /**
     * Returns the proj4 description of the projection used
     * @return
     */
	std::string proj4();

	//holds the spatial search tree
	//http://doc.cgal.org/latest/Spatial_searching/index.html
	boost::shared_ptr<Tree> dD_tree;

    bool _terrain_deformed;

    /**
     * Minimium z elevation in mesh
     * @return
     */
    double min_z();

    /**
     * Maximum elevation in mesh
     * @return
     */
    double max_z();
private:
    size_t _num_faces; //number of faces
    size_t _num_vertex; //number of rows in the original data matrix. useful for exporting to matlab, etc
    K::Iso_rectangle_2 _bbox;
	bool _is_geographic;
	int _UTM_zone;

	std::string _srs_wkt;
	//holds the vtk ugrid if we are outputing to vtk formats
	vtkSmartPointer<vtkUnstructuredGrid> _vtk_unstructuredGrid;

	//holds the vectors we use to create the vtu file
	std::map<std::string, vtkSmartPointer<vtkFloatArray> > data;
	std::map<std::string, vtkSmartPointer<vtkFloatArray> > vectors;

    // min and max elevations
    double _min_z;
    double _max_z;

    //If the triangulation is traversed using the finite_faces_begin/end iterators, the determinism of the order of traversal is not guaranteed
    //as well, it seems to prevent openmp for applying parallism to the for-loops. Therefore, we will just store a predefined list of faces and vertex handles
    //that allows us to traverse the triangulation in a deterministic order, as well as play nice with openmp
    std::vector< Delaunay::Face_handle  > _faces;
    std::vector< Delaunay::Vertex_handle > _vertexes;
    
#ifdef NOMATLAB
    //ptr to the matlab engine
    boost::shared_ptr<maw::matlab_engine> _engine;
    boost::shared_ptr<maw::graphics> _gfx;
#endif

};

/**
* \typedef mesh
* Provides a convenience typedef for passing around mesh pointers.
*/
typedef boost::shared_ptr<triangulation> mesh;

/**
* \class rect
* Sub-bounding box rectangle.
*/
class rect
{
public:
    /**
    * Creates a bounding box rectangle with the given coordinates
    * \param coord x,y coordinate pairs
    */
	rect( arma::mat* coord )
	{
		this->coord = coord;
	}
	~rect()
	{
		delete coord;
	}
	arma::mat* coord;
    tbb::concurrent_vector<mesh_elem> triangles;
};

/**
* \class segmented_AABB
* \brief segmented axis aligned bounding box
*
* A segmented axis aligned bounding box for quickly spatially binning triangles. A poor man's k-tree.
* 	         4       top      3
*	         +--------+--------+
*	         |        |        |
*	         |        |        |
*	   left  +--------+--------+    right
*	         |        |        |
*	         |        |        |
*	         +-----------------+
*	         1     bottom      2
*	 +y
*	 ^
*	 |
*	 +-> +x
*
*/
class segmented_AABB
{
public:
	segmented_AABB();
	~segmented_AABB();

    /**
    * Creates a segmented AABB for the given triangulation of size rows x cols
    * \param domain A triangulation mesh
    * \param rows number of rows
    * \param cols number of cols
    */
	void make( triangulation* domain,  size_t rows, size_t cols);

    /**
    * Returns the sub rectangle
    * \param row
    * \param col
    * \return The sub rectangle
    */
	rect* get_rect( size_t row, size_t col);

    /**
    * \var n_rows
    * Number of rows in the AABB
    */
	size_t n_rows;

    /**
    * \var n_cols
    * Number of cols in the AABB
    */
	size_t n_cols;

    /**
    * Determins if a point exists within a given rectangle
    * \param v Vertex handle toa  point (vertex) in the triangulation
    * \param r A subrect
    * \returns true if point lies whitin r
    */
    bool pt_in_rect(Delaunay::Vertex_handle v, rect* r);

    /**
    * Determines if a point given by x y coords exists in a given rectangle
    * \param x coord
    * \param y coord
    * \param r sub rectangle
    * \returns true if point is within r
    */
	bool pt_in_rect(double x, double y, rect* r);

private:
	tbb::concurrent_vector<tbb::concurrent_vector<rect*> > m_grid;

};


template < class Gt, class Fb >
bool face<Gt, Fb>::has_parameter(std::string key)
{
    return _parameters.find( key ) != _parameters.end();
}

template < class Gt, class Fb >
void face<Gt, Fb>::set_parameter(std::string key, double value)
{
    _parameters[key] = value;
}

template < class Gt, class Fb >
double face<Gt, Fb>::get_parameter(std::string key)
{
    return _parameters[key];
};

template < class Gt, class Fb >
bool face<Gt, Fb>::has_initial_condition(std::string key)
{
    return _initial_conditions.find( key ) != _initial_conditions.end();
}

template < class Gt, class Fb >
void face<Gt, Fb>::set_initial_condition(std::string key, double value)
{
    _initial_conditions[key] = value;
}

template < class Gt, class Fb >
double face<Gt, Fb>::get_initial_condition(std::string key)
{
    return _initial_conditions[key];
};

template < class Gt, class Fb >
std::vector<std::string>  face<Gt, Fb>::parameters()
{
    std::vector<std::string> params;
    for(auto& itr : _parameters)
    {
        params.push_back(itr.first);
    }
    return params;
};

template < class Gt, class Fb >
std::vector<std::string>  face<Gt, Fb>::initial_conditions()
{
    std::vector<std::string> ics;
    for(auto& itr : _initial_conditions)
    {
        ics.push_back(itr.first);
    }
    return ics;
};

template < class Gt, class Fb >
std::vector<std::string>  face<Gt, Fb>::vectors()
{
    std::vector<std::string> vecs;
    for(auto& itr : _module_face_vectors)
    {
        vecs.push_back(itr.first);
    }
    return vecs;
};

template < class Gt, class Fb >
face<Gt, Fb>::~face()
{

    for(auto itr : _module_face_data)
    {
        delete itr.second;
    }
};
template < class Gt, class Fb >
face<Gt, Fb>::face()
{
    _slope = -1;
    _azimuth = -1;
    _data = boost::make_shared<timeseries>();

    _center = NULL;
    _normal = NULL;

    _area = -1.;

    _is_geographic = false;


}

template < class Gt, class Fb>
face<Gt, Fb>::face(Vertex_handle v0,
                   Vertex_handle v1,
                   Vertex_handle v2)
        : Fb(v0, v1, v2)
{
    _slope = -1;
    _azimuth = -1;
    _data = boost::make_shared<timeseries>();
    _center = NULL;
    _normal = NULL;
    _area = -1.;
    _is_geographic = false;
}

template < class Gt, class Fb >
face<Gt, Fb>::face(Vertex_handle v0,
                   Vertex_handle v1,
                   Vertex_handle v2,
                   Face_handle n0,
                   Face_handle n1,
                   Face_handle n2)
        : Fb(v0, v1, v2, n0, n1, n2)
{
    _slope = -1;
    _azimuth = -1;
    _data = boost::make_shared<timeseries>();
    _center = NULL;
    _normal = NULL;
    _area = -1.;
    _is_geographic = false;

}

template < class Gt, class Fb >
face<Gt, Fb>::face(Vertex_handle v0,
                   Vertex_handle v1,
                   Vertex_handle v2,
                   Face_handle n0,
                   Face_handle n1,
                   Face_handle n2,
                   bool c0,
                   bool c1,
                   bool c2)
        : Fb(v0, v1, v2, n0, n1, n2)
{
    _slope = -1;
    _azimuth = -1;
    _data = boost::make_shared<timeseries>();
    _center = NULL;
    _normal = NULL;
    _area = -1.;
    _is_geographic = false;
}

template < class Gt, class Fb>
double face<Gt, Fb>::aspect()
{
    if (_azimuth == -1)
    {
        if(!_normal)
            this->normal();
        //convert normal to spherical
//        double r = sqrt(_normal[0] * _normal[0] + _normal[1] * _normal[1] + _normal[2] * _normal[2]);
//        double theta = acos(_normal[2] / r);

        //y=north
        double phi = atan2((*_normal)[1], (*_normal)[0]); // + M_PI /*-3.14159/2*/; //south == 0
        _azimuth = phi - M_PI / 2.0; //set north = 0

        if (_azimuth < 0.0)
            _azimuth += 2.0 * M_PI;
    }

    return _azimuth;
}

template < class Gt, class Fb>
double face<Gt, Fb>::slope()
{
    if (_slope == -1)
    {
        if(!_normal)
            this->normal();

        //z surface normal
        arma::vec n(3);
        n(0) = 0.0; //x
        n(1) = 0.0; //y
        n(2) = 1.0;

        arma::vec normal(3);
        normal(0) = (*_normal)[0];
        normal(1) = (*_normal)[1];
        normal(2) = (*_normal)[2];

        _slope = acos(arma::norm_dot(normal, n));
    }

    return _slope;
}

template < class Gt, class Fb>
bool face<Gt, Fb>::intersects(face<Gt, Fb>::Face_handle fh)
{
    //does t contain my points?
    bool intersect = fh->contains(this->center());
    if (!intersect)
    {
        if (fh->contains(this->vertex(0)->point()) ||
            fh->contains(this->vertex(1)->point()) ||
            fh->contains(this->vertex(2)->point()))
        {
            intersect = true;
        }
    }
    return intersect;
}

template < class Gt, class Fb>
const typename face<Gt, Fb>::Face_handle face<Gt, Fb>::find_closest_face(double azmimuth,double distance)
{
    return _domain->find_closest_face(math::gis::point_from_bearing(center(), azmimuth, distance));
};

template < class Gt, class Fb>
Vector_3 face<Gt, Fb>::normal()
{
    if(!_normal)
    {
        if(_is_geographic)
        {

//            OGRSpatialReference monUtm;
//
//            OGRSpatialReference monGeo;
//            monGeo.SetWellKnownGeogCS("WGS84");


            CGAL::Point_3<K> v0(this->vertex(0)->point()[0]*100000., this->vertex(0)->point()[1]*100000.,this->vertex(0)->point()[2]);
            CGAL::Point_3<K> v1(this->vertex(1)->point()[0]*100000., this->vertex(1)->point()[1]*100000.,this->vertex(1)->point()[2]);
            CGAL::Point_3<K> v2(this->vertex(2)->point()[0]*100000., this->vertex(2)->point()[1]*100000.,this->vertex(2)->point()[2]);

            _normal = boost::make_shared<Vector_3>(CGAL::unit_normal(v0, v1, v2));

        }
        else
            _normal = boost::make_shared<Vector_3>(CGAL::unit_normal(this->vertex(0)->point(), this->vertex(1)->point(), this->vertex(2)->point()));
    }

//    CGAL::Point_3<K> v0(this->vertex(0)->point()[0]*100000., this->vertex(0)->point()[1]*100000.,this->vertex(0)->point()[2]);
//    CGAL::Point_3<K> v1(this->vertex(1)->point()[0]*100000., this->vertex(1)->point()[1]*100000.,this->vertex(1)->point()[2]);
//    CGAL::Point_3<K> v2(this->vertex(2)->point()[0]*100000., this->vertex(2)->point()[1]*100000.,this->vertex(2)->point()[2]);

//    CGAL::Point_3<CGAL::Exact_prediDocuments/PhD/code/CHM/cates_exact_constructions_kernel > v0_noscale(this->vertex(0)->point()[0], this->vertex(0)->point()[1],this->vertex(0)->point()[2]);
//    CGAL::Point_3<CGAL::Exact_predicates_exact_constructions_kernel > v1_noscale(this->vertex(1)->point()[0], this->vertex(1)->point()[1],this->vertex(1)->point()[2]);
//    CGAL::Point_3<CGAL::Exact_predicates_exact_constructions_kernel > v2_noscale(this->vertex(2)->point()[0], this->vertex(2)->point()[1],this->vertex(2)->point()[2]);
//
//    CGAL::Exact_predicates_exact_constructions_kernel::Vector_3 un1 = CGAL::unit_normal(v0_noscale, v1_noscale, v2_noscale);
//    Vector_3 un2 = CGAL::unit_normal(v0, v1, v2);
    return *_normal;
}

template < class Gt, class Fb>
Point_3 face<Gt, Fb>::center()
{
    if (!_center)
    {
        _center = boost::make_shared<Point_3>(CGAL::centroid(this->vertex(0)->point(), this->vertex(1)->point(), this->vertex(2)->point()));
        _x=_center->x();
        _y=_center->y();
        _z=_center->z();
    }
    //return CGAL::centroid(this->vertex(0)->point(), this->vertex(1)->point(), this->vertex(2)->point());
    return *_center;

}
template < class Gt, class Fb>
bool face<Gt, Fb>::contains(Point_3 p)
{
    return this->contains(p.x(), p.y());
}

template < class Gt, class Fb>
bool face<Gt, Fb>::contains(double x, double y)
{

    double x1 = this->vertex(1)->point().x();
    double y1 = this->vertex(1)->point().y();

    double x2 = this->vertex(2)->point().x();
    double y2 = this->vertex(2)->point().y();

    double x3 = this->vertex(0)->point().x();
    double y3 = this->vertex(0)->point().y();

    double lambda1 = ((y2 - y3)*(x - x3)+(x3 - x2)*(y - y3)) / ((y2 - y3)*(x1 - x3)+(x3 - x2)*(y1 - y3));

    if( !(lambda1 > 0.0 && lambda1 < 1.0) )
        return false; //bail early if possible

    double lambda2 = ((y3 - y1)*(x - x3)+(x1 - x3)*(y - y3)) / ((y3 - y1)*(x2 - x3)+(x1 - x3)*(y2 - y3));
    double lambda3 = 1.0 - lambda1 - lambda2;

    return lambda1 > 0.0 && lambda1 < 1.0
           && lambda2 > 0.0 && lambda2 < 1.0
           && lambda3 > 0.0 && lambda3 < 1.0;

}
template < class Gt, class Fb>
std::vector<std::string> face<Gt, Fb>::variables()
{
    return _data->list_variables();
}


template < class Gt, class Fb>
bool face<Gt, Fb>::has(const std::string& variable)
{
    return _itr->has(variable);
};

template < class Gt, class Fb>
void face<Gt, Fb>::set_face_data(const std::string& variable, double data)
{
    _itr->set(variable, data);
}

template < class Gt, class Fb>
double face<Gt, Fb>::face_data(const std::string& variable)
{
    return _itr->get(variable);
}

template < class Gt, class Fb>
Vector_3 face<Gt, Fb>::face_vector(const std::string& variable)
{
    auto res = _module_face_vectors.find(variable);
    if(res == _module_face_vectors.end())
    {
        return Vector_3(nan(""),nan(""),nan(""));
    }
//        BOOST_THROW_EXCEPTION( forcing_lookup_error() << errstr_info("Variable " + variable + " does not exist."));

    return _module_face_vectors[variable];
};

template < class Gt, class Fb>
void face<Gt, Fb>::init_time_series(std::set<std::string> variables, timeseries::date_vec datetime)
{
    _data->init(variables, datetime);

    _itr = _data->begin();
}

template < class Gt, class Fb>
timeseries::variable_vec face<Gt, Fb>::face_time_series(std::string ID)
{
    return _data->get_time_series(ID);
}

template < class Gt, class Fb>
void face<Gt, Fb>::next()
{
    ++_itr;
}

template < class Gt, class Fb>
void face<Gt, Fb>::reset_to_begining()
{
    _itr = _data->begin();
}

template < class Gt, class Fb>
void face<Gt, Fb>::to_file(std::string fname)
{
    _data->to_file(fname);
}

template < class Gt, class Fb>
double face<Gt, Fb>::get_x()
{
    if(!_center)
        this->center();

    return _x;
}

template < class Gt, class Fb>
double face<Gt, Fb>::get_y()
{
    if(!_center)
        this->center();


    return _y;
}

template < class Gt, class Fb>
double face<Gt, Fb>::get_z()
{
    if(!_center)
        this->center();

    return _z;
}
template < class Gt, class Fb>
boost::shared_ptr<timeseries> face<Gt, Fb>::get_underlying_timeseries()
{
    return _data;
}

template < class Gt, class Fb>
timeseries::iterator face<Gt, Fb>::now()
{
    return _itr;
}


template < class Gt, class Vb>
template<typename T>
T* face<Gt, Vb>::make_module_data(std::string module)
{


    auto it = _module_face_data.find(module);

    //we don't already have this, make a new one.
    if(it == _module_face_data.end())
    {
        T* data = new T;
        _module_face_data[module] = data;
    }

    return get_module_data<T>(module);
}


template < class Gt, class Fb>
template < typename T>
T* face<Gt, Fb>::get_module_data(std::string module)
{
    auto it = _module_face_data.find(module);

    //we don't already have this
    if(it == _module_face_data.end())
    {
        BOOST_THROW_EXCEPTION(module_data_error() << errstr_info ("No data for module " + module));
    }

    return dynamic_cast<T*>(it->second);


}
template < class Gt, class Fb>
void face<Gt, Fb>::set_face_vector(const std::string& variable, Vector_3 v)
{
    _module_face_vectors[variable] = v;
};

template < class Gt, class Fb>
void face<Gt, Fb>::remove_face_data(std::string module)
{
    try
    {
        auto ptr =  _module_face_data[module];
        //split it up to catch the exception if ID not found
        delete ptr;
        _module_face_data.erase(module);
    }
    catch(...)
    {
        BOOST_THROW_EXCEPTION(module_data_error() << errstr_info ("No data for module " + module));
    }
};
template < class Gt, class Fb>
void face<Gt, Fb>::set_module_data(std::string module, face_info *fi)
{
    _module_face_data[module] = fi;
}
template < class Gt, class Fb>
double face<Gt, Fb>::get_area()
{
    if(_area == -1)
    {
        // supports geographic
        if(has_parameter("area"))
        {
            _area = get_parameter("area");
        }
        else
        {
            auto x1 = this->vertex(0)->point().x();
            auto y1 = this->vertex(0)->point().y();

            auto x2 = this->vertex(1)->point().x();
            auto y2 = this->vertex(1)->point().y();

            auto x3 = this->vertex(2)->point().x();
            auto y3 = this->vertex(2)->point().y();

//        auto lol= CGAL::area(this->vertex(0)->point(),this->vertex(1)->point(),this->vertex(2)->point());
            _area = 0.5 * fabs( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );
        }
    }


    return _area;
};