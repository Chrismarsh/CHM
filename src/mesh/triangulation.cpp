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


#include "triangulation.h"

void triangulation::create_delaunay(arma::vec* x, arma::vec* y, arma::vec* z)
{
    if (m_engine)
    {

        int num_nodes = x->n_rows;

        mxArray* xy = mxCreateDoubleMatrix(num_nodes, 3, mxREAL);
        double* ptr = mxGetPr(xy);

        //create temp arrays to send to matlab
        for (int row = 0; row < num_nodes; row++)
        {
            //matlab is col major storage	
            ptr[row + 0 * num_nodes] = (*x)(row);
            ptr[row + 1 * num_nodes] = (*y)(row);
        }
        m_engine->put("xy", xy);
        m_engine->evaluate("tri=DelaunayTri(xy(:,1),xy(:,2))");
        //clean up our temp array
        mxDestroyArray(xy);
        xy = NULL;
        ptr = NULL;

        //change this later to the struct lookup
        m_engine->evaluate("t=tri.Triangulation");
        //get our triangulation structure from matlab

        maw::d_mat tri = m_engine->get_double_matrix("t");
        //	mxArray* tri = m_engine->get("t");

        //will be n * 3
        // 		const mwSize* size = mxGetDimensions(tri);
        //  		m_size = size[0]; //first element is the # of rows
        m_size = tri->n_rows;
        m_data_size = num_nodes;

        //double* t = mxGetPr(tri);

        for (size_t i = 0; i < m_size; i++)
        {
            //col major lookup!!!
            //get the row index data from the matlab triangulation structure
            //note: all these indexes will be off by 1 because matlab indexing starts at 1 and not 0
            int v1 = int((*tri)(i + 0 * m_size));
            v1 -= 1;
            int v2 = int((*tri)(i + 1 * m_size));
            v2 -= 1;
            int v3 = int((*tri)(i + 2 * m_size));
            v3 -= 1;

            point vertex1, vertex2, vertex3;
            vertex1.x = (*x)(v1);
            vertex1.y = (*y)(v1);
            vertex1.z = (*z)(v1);

            vertex2.x = (*x)(v2);
            vertex2.y = (*y)(v2);
            vertex2.z = (*z)(v2);

            vertex3.x = (*x)(v3);
            vertex3.y = (*y)(v3);
            vertex3.z = (*z)(v3);

            m_triangles.push_back(new triangle(vertex1, vertex2, vertex3, 0));
            m_triangles[i]->global_id[0] = v1 + 1; //+1 for matlab
            m_triangles[i]->global_id[1] = v2 + 1;
            m_triangles[i]->global_id[2] = v3 + 1;

            m_triangles[i]->triangle_id = i;
        }

        //clean up in matlab
        m_engine->evaluate("clear t xy tri;");
    }
    LOG_DEBUG << "Created a mesh with " + boost::lexical_cast<std::string>(m_size) +" triangles";
}

size_t triangulation::size()
{
    return m_size;
}

triangulation::triangulation(maw::matlab_engine* engine)
{
    m_engine = engine;
    m_size = 0;

}

triangulation::~triangulation()
{
    for (auto it = m_triangles.begin(); it != m_triangles.end(); it++)
    {
        delete *it;
    }

}

triangle& triangulation::operator()(size_t t)
{
    return *(m_triangles[t]);
}

void triangulation::set_vertex_data(arma::mat& data)
{
    if (data.n_cols != 3)
        throw std::runtime_error("Wrong number of columns");

    //loop over all triangles

    for (auto it = m_triangles.begin(); it != m_triangles.end(); ++it)
    {
        size_t v1 = (*it)->global_id[0];
        size_t v2 = (*it)->global_id[1];
        size_t v3 = (*it)->global_id[2];

        point vertex1, vertex2, vertex3;
        vertex1.x = data(v1 - 1, 0);
        vertex1.y = data(v1 - 1, 1);
        vertex1.z = data(v1 - 1, 2);

        vertex2.x = data(v2 - 1, 0);
        vertex2.y = data(v2 - 1, 1);
        vertex2.z = data(v2 - 1, 2);

        vertex3.x = data(v3 - 1, 0);
        vertex3.y = data(v3 - 1, 1);
        vertex3.z = data(v3 - 1, 2);

        (*it)->set_vertex_values(vertex1, vertex2, vertex3);

    }

}

void triangulation::compute_face_normals()
{
    LOG_DEBUG << "Sending elevation data to matlab";
    m_engine->put_double_matrix("mxDomain",this->mf_elevation_data());
    LOG_DEBUG << "Sending triangulation data to matlab";
    m_engine->put_double_matrix("tri",this->mf_tri_matrix());

    
    m_engine->evaluate("[NormalVx NormalVy NormalVz PosVx PosVy PosVz]=computeNormalVectorTriangulation(mxDomain,tri,'center-cells');");
    m_engine->evaluate("normals=[NormalVx NormalVy NormalVz];");
    m_engine->evaluate("center=[PosVx PosVy PosVz]");
    m_engine->evaluate("clear PosVx PosVy PosVz NormalVx NormalVy NormalVz");

    maw::d_mat normals = m_engine->get_double_matrix("normals");
    maw::d_mat centers = m_engine->get_double_matrix("center");

    size_t counter = 0;
    for (auto it = m_triangles.begin(); it != m_triangles.end(); it++)
    {
        arma::vec n(3);
        n(0) = normals->row(counter)(0);
        n(1) = normals->row(counter)(1);
        n(2) = normals->row(counter)(2);

        arma::vec c(3);
        c(0) = centers->row(counter)(0);
        c(1) = centers->row(counter)(1);
        c(2) = centers->row(counter)(2);

        (*it)->set_facenormal(n);
        point p(c(0), c(1), c(2));
        (*it)->center = p;
        ++counter;
    }
    m_engine->evaluate("clear normals center mxDomain tri");
    LOG_DEBUG << "Done calculating face normals";
}

triangle* triangulation::find_containing_triangle(double x, double y)
{
    for (std::vector<triangle*>::iterator it = m_triangles.begin(); it != m_triangles.end(); it++)
    {
        if ((*it)->contains(x, y))
            return *it;
    }

    return NULL;
}

maw::d_mat triangulation::mf_tri_matrix()
{
    maw::d_mat tri(new arma::mat(m_size, 3));
    size_t i = 0;
    for (auto it = m_triangles.begin(); it != m_triangles.end(); it++)
    {
        (*tri)(i, 0) = (*it)->global_id[0];
        (*tri)(i, 1) = (*it)->global_id[1];
        (*tri)(i, 2) = (*it)->global_id[2];
        i++;
    }

    return tri;
}

maw::d_mat triangulation::mf_elevation_data()
{
    maw::d_mat data(new arma::mat(m_data_size, 3)); //special case as we actually have *less*

    for (auto& it : m_triangles)
    {

        (*data)(it->global_id[0] - 1, 0) = it->get_vertex(0).x;
        (*data)(it->global_id[0] - 1, 1) = it->get_vertex(0).y;
        (*data)(it->global_id[0] - 1, 2) = it->get_vertex(0).z;

        (*data)(it->global_id[1] - 1, 0) = it->get_vertex(1).x;
        (*data)(it->global_id[1] - 1, 1) = it->get_vertex(1).y;
        (*data)(it->global_id[1] - 1, 2) = it->get_vertex(1).z;

        (*data)(it->global_id[2] - 1, 0) = it->get_vertex(2).x;
        (*data)(it->global_id[2] - 1, 1) = it->get_vertex(2).y;
        (*data)(it->global_id[2] - 1, 2) = it->get_vertex(2).z;
    }

    return data;
}

maw::d_vec triangulation::mf_face_data(std::string ID)
{
    maw::d_vec data(new arma::vec(m_size));

    int i = 0;
    for (auto& it : m_triangles)
    {
        double d = it->get_face_data(ID);
        (*data)(i) = d;
        ++i;

    }

    return data;
}