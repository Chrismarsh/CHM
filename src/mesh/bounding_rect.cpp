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

#include "bounding_rect.h"



void bounding_rect::make(const maw::d_vec x, const maw::d_vec y, int n_rows, int n_cols)
{
	int n_segments = n_cols;
	int n_v_segments = n_rows;

	arma::vec midpoint_bottom_x(n_segments);
	arma::vec midpoint_top_x(n_segments);

	arma::vec midpoint_bottom_y(n_segments);
	arma::vec midpoint_top_y(n_segments);

	m_engine->put_double_vector("BBR_x",x);
	m_engine->put_double_vector("BBR_y",y);

	
	m_engine->evaluate("[bbx,bby,~,~]=minboundrect(BBR_x(:),BBR_y(:));");


	maw::d_vec bbx = (m_engine->get_double_vector("bbx"));
	maw::d_vec bby = (m_engine->get_double_vector("bby"));

	m_engine->evaluate("clear bbx bby BBR_x BBR_y");

	

	this->n_rows = n_rows;
	this->n_cols = n_cols;

// 	bbx = new arma::vec(bbx);
// 	bby = new arma::vec(bby);

	//need to construct new axis aligned BBR
	arma::u32 index;
	
	//left most pt
	double x_left = bbx->min(index);
	double y_left = bby->operator()(index);

	//right most pt
	double x_right = bbx->max(index);
	double y_right = bby->operator()(index);

	//bottom most pt
	double y_bottom = bby->min(index);
	double x_bottom = bbx->operator()(index);

	//top most pt
	double y_top = bby->max(index);
	double x_top = bbx->operator()(index);


	//horizontal step size
	double h_dx = (x_right - x_left) / n_segments;
	double v_dy = (y_top - y_bottom) / n_v_segments;


	//resize the row dimension
	m_grid.resize(n_v_segments);

	//start top left
	//loop over rows
	for(int i = 0; i < n_v_segments; i++)
	{
		//set the y coordinate to the current row top left coordinate
		double h_y = y_top - v_dy*i;
		double h_x = x_left;

		//set number of cols
		m_grid[i].resize(n_segments);


		//loop over columns
		for(int j = 0; j < n_segments; j++)
		{
			
			arma::mat* t = new arma::mat(5,2);
				

			*t << h_x << h_y - v_dy << arma::endr // bottom left
			   << h_x + h_dx << h_y - v_dy << arma::endr //bottom right
			   << h_x + h_dx << h_y << arma::endr // top right
			    << h_x << h_y << arma::endr //top left
				<< h_x << h_y - v_dy << arma::endr; // bottom left


			m_grid[i][j] = new rect(t);
			h_x = h_x + h_dx;
		}
	}

// 	delete bbx;
// 	delete bby;
	
	
/*
	double m = (bby(1)-bby(0))/(bbx(1)-bbx(0));

	double step=(bbx(2)-bbx(3))/n_segments;
	for (int i = 0; i<n_segments;i++)
	{
		double xpos_bottom = bbx(0)+step*(i+1);
		double xpos_top     =bbx(3)+step*(i+1);
		         
		midpoint_bottom_x[i]=xpos_bottom;
		midpoint_bottom_y[i]=m*(xpos_bottom-bbx(0))+bby(0); //2pt line eqn
			       
		midpoint_top_x[i]=xpos_top;
		midpoint_top_y[i]=m*(xpos_top-bbx(3))+bby(3); //2pt line eqn

	}

	m_rectangles.reserve(n_segments);

	for (int i = 0; i<n_segments;i++)
	{
		if (i == 0)
		{
			arma::mat* c = new arma::mat(5,2);
			arma::mat& coord = *c;
			coord	<< bbx(0)				 << bby(0)					<< arma::endr  //bottom left
					<< midpoint_bottom_x[0]	 << midpoint_bottom_y[0]	<< arma::endr  //bottom right mid point
					<< midpoint_top_x[0]	 << midpoint_top_y[0]		<< arma::endr  //top right mid point
					<< bbx(3)				 << bby(3)					<< arma::endr  //top right
					<< bbx(0)				 << bby(0)					<< arma::endr;
			
			m_rectangles.push_back(new rect(c));
		}
		//       last rect
		else if( i==n_segments-1)
		{
			arma::mat* c = new arma::mat(5,2);
			arma::mat& coord = *c;

			coord	<< 	midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr  //bottom left
					<< 	bbx(1)					<<  bby(1)						<< arma::endr  //bottom right mid point
					<<  bbx(2)					<< 	bby(2)						<< arma::endr  //top right mid point
					<< 	midpoint_top_x[i-1]		<< 	midpoint_top_y[i-1]			<< arma::endr  //top right
					<<  midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr;

			m_rectangles.push_back(new rect(c));

		}
		else
		{

			arma::mat* c = new arma::mat(5,2);
			arma::mat& coord = *c;

			coord	<< 	midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr  //bottom left
					<< 	midpoint_bottom_x[i]	<<  midpoint_bottom_y[i]		<< arma::endr  //bottom right mid point
					<<  midpoint_top_x[i]		<< 	midpoint_top_y[i]			<< arma::endr  //top right mid point
					<< 	midpoint_top_x[i-1]		<< 	midpoint_top_y[i-1]			<< arma::endr  //top right
					<<  midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr;

			m_rectangles.push_back(new rect(c));
		}
	
	}
	}*/
}

rect* bounding_rect::get_rect( int i, int j )
{
		return m_grid[i][j];
//	 return m_rectangles.at(i);
}

bounding_rect::bounding_rect( maw::matlab_engine* m_engine )
{
	this->m_engine = m_engine;

}

bool bounding_rect::pt_in_rect( double x, double y, rect* r )
{
	for (int i = 0; i<4; i++)
	{
		double x0 = r->coord->operator()(i,0);
		double y0 = r->coord->operator()(i,1);
		double x1 = r->coord->operator()(i+1,0);
		double y1 = r->coord->operator()(i+1,1);

		double pt = ((y - y0)*(x1 - x0) - (x - x0)*(y1 - y0));
		if (pt < 0)
			return false;
		else if (pt ==0)
			return false;
	}
	return true;
}

bounding_rect::~bounding_rect()
{
	for(auto it=m_grid.begin();it!=m_grid.end();it++)
	{
		for(auto jt=it->begin(); jt!= it->end(); jt++)
		{
			delete *jt;
		}
	}

	m_engine = NULL;
}