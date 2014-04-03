// % Copyright (c) 2011, Chris Marsh
// 	% All rights reserved.
// 	% 
// 	% Redistribution and use in source and binary forms, with or without 
// 	% modification, are permitted provided that the following conditions are 
// 	% met:
// % 
// 	%     * Redistributions of source code must retain the above copyright 
// 	%       notice, this list of conditions and the following disclaimer.
// 	%     * Redistributions in binary form must reproduce the above copyright 
// 	%       notice, this list of conditions and the following disclaimer in 
// 	%       the documentation and/or other materials provided with the distribution
// 	%     * Neither the name of the University of Saskatchewan nor the names 
// 	%       of its contributors may be used to endorse or promote products derived 
// 	%       from this software without specific prior written permission.
// 	%       
// 	% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// 	% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// 	% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// 	% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
// 	% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
// 	% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
// 	% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
// 	% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
// 	% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
// 	% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// 	% POSSIBILITY OF SUCH DAMAGE.


#include "graphics.h"

namespace maw
{


//returns the handle from matlab
//vertices should be in the form [x(:) y(:) z(:)]
//faces should be in the form  tri.Triangulation
double graphics::plot_patch( std::string vertices, std::string faces, std::string face_data )
{
	std::string command = std::string("patch_handle = patch('Vertices',") + 	vertices + 
		std::string(",'Faces',") + 	faces +
		std::string(",'facevertexcdata',") + face_data +
		std::string(",'facecolor','flat', 'edgecolor','none');");//[201/256 201/256 201/256]
	m_engine->evaluate(command.c_str());

	mxArray* handle =  m_engine->get("patch_handle");

	if(handle)
		return (mxGetScalar(handle));
	else
		return -1;
}

graphics::graphics( matlab_engine *engine )
{
	m_engine = engine;
}

double graphics::update_patch( double handle, std::string vertices, std::string face_data )
{
	std::string command = std::string("set(")+ boost::lexical_cast<std::string,double>(handle) + std::string(",'Vertices',") + vertices  + std::string(",'facevertexcdata',") + face_data + std::string(")");
	m_engine->evaluate(command.c_str());
	//m_engine->evaluate("refreshdata");
	mxArray* h =  m_engine->get("patch_handle");

	return (mxGetScalar(h));

}

double graphics::add_title( std::string title, int fontsize /*= 14*/,std::string color /*='black'*/ )
{
	m_engine->evaluate("if exist('ht','var')==1 delete(ht.th); end");
	std::string command = std::string("ht = mtit('") + title + std::string("','fontsize',") + boost::lexical_cast<std::string,int>(fontsize) + std::string(")");
	m_engine->evaluate(command);
	mxArray* ht =  m_engine->get("ht");
	m_engine->evaluate(std::string("set(ht.th,'color','") + color + std::string("')"));
	return (mxGetScalar(ht));
}

void graphics::hold_on()
{
	m_engine->evaluate("hold on");
}

void graphics::hold_off()
{
	m_engine->evaluate("hold off");
}

double graphics::plot_line( std::string x, std::string y, std::string options/*=""*/ )
{
	double ret = -1;
	if (m_engine)
	{
		std::string command = std::string("plot_handle=plot(") + x +std::string(",") + y + std::string(",'-'");

		if (options == "")
		{
			command += std::string(")");
		}
		else
		{
			command += std::string(",") + options + std::string(")");
		}


		m_engine->evaluate(command);
		ret = m_engine->get_scaler("plot_handle");
		m_engine->evaluate("clear plot_handle");
		
	}
	return ret;
}

double graphics::plot_line(const d_vec x, const d_vec y,std::string options/*=""*/ )
{
	m_engine->put_double_vector("x",x);
	m_engine->put_double_vector("y",y);

	std::string x_name = "x";
	std::string y_name = "y";

	double handle = plot_line(x_name,y_name,options);

	m_engine->evaluate("clear x y plot_handle");

	return handle;

}

graphics::~graphics()
{

}

void graphics::save_to_file( std::string fname )
{
	m_engine->evaluate(std::string("export_figure('") + fname +std::string("','-png');"));
}

void graphics::colorbar( std::string toggle/*="on"*/ )
{
	if(toggle == "on")
		m_engine->evaluate("colorbar");
	else
		m_engine->evaluate("colorbar('off')");
}

}