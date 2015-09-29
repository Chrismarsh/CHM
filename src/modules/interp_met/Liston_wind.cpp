#include "Liston_wind.hpp"

Liston_wind::Liston_wind( std::string ID)
{

   _depends_from_met->push_back("u");
    //_depends_from_met->push_back("u");


   _provides->push_back("u");
    _provides->push_back("Liston Curvature");

    this->ID = ID;
    _parallel_type = parallel::data;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

void Liston_wind::init(mesh domain)
{
    for (Delaunay::Finite_faces_iterator fit = domain->finite_faces_begin();
         fit != domain->finite_faces_end(); ++fit)
    {
        Delaunay::Face_handle face = fit;
        face->coloured = false;
    }


    double curmax = -99999.0;

    for (Delaunay::Finite_faces_iterator fit = domain->finite_faces_begin();
         fit != domain->finite_faces_end(); ++fit)
    {
        Delaunay::Face_handle face = fit;



        std::set<Point_3> myPoints;
        std::stack<Delaunay::Face_handle> neighbours;


        std::vector<Delaunay::Face_handle> neighbours_to_uncolor;
        int cur_depth = 0;
        int max_depth = 10;

        neighbours.push(face);

        while (!neighbours.empty() && (cur_depth <= max_depth))
        {
            Delaunay::Face_handle nface = neighbours.top();
            neighbours.pop();
            //this triangle
            //origin of fitting coord system = first input data point
            myPoints.insert(nface->center());
            for (int i = 0; i < 3; i++)
            {
                myPoints.insert(nface->vertex(i)->point());
            }

            //add all of f's neighbours
            for (int k = 0; k < 3; k++)
            {
                auto n = nface->neighbor(k);
                if (!domain->is_infinite(n) && !n->coloured)
                {
                    n->coloured = true;
                    neighbours_to_uncolor.push_back(n);
                    neighbours.push(n);
                }

            }
            cur_depth++;
        }

        for (auto it:neighbours_to_uncolor)
        {
            it->coloured = false;
        }
        neighbours_to_uncolor.clear();

        Point_3 me = face->center();
        Delaunay::Face_handle north;
        Delaunay::Face_handle south;
        Delaunay::Face_handle east;
        Delaunay::Face_handle west;

        Delaunay::Face_handle northeast;
        Delaunay::Face_handle northwest;
        Delaunay::Face_handle southeast;
        Delaunay::Face_handle southwest;


        double distance = 300.0;


        north = domain->locate_face(me.x(), me.y() + distance);
        if (north == NULL || domain->is_infinite(north))
            north = face;

        south = domain->locate_face(me.x(), me.y() - distance);
        if (south == NULL || domain->is_infinite(south))
            south = face;

        west = domain->locate_face(me.x() - distance, me.y());
        if (west == NULL || domain->is_infinite(west))
            west = face;

        east = domain->locate_face(me.x() + distance, me.y());
        if (east == NULL || domain->is_infinite(east))
            east = face;

        double z = face->get_z();
        double zw = west->get_z();
        double ze = east->get_z();
        double zs = south->get_z();
        double zn = north->get_z();

        double znw = 0.;
        double zne = 0.;
        double zse = 0.;
        double zsw = 0.;

        northeast = domain->locate_face(me.x() + distance, me.y() + distance);
        if (northeast == NULL || domain->is_infinite(northeast))
            zne = (ze + zn) / 2.;
        else
            zne = northeast->get_z();

        northwest = domain->locate_face(me.x() - distance, me.y() + distance);
        if (northwest == NULL || domain->is_infinite(south))
            znw = (zw + zn) / 2.;
        else
            znw = northwest->get_z();

        southeast = domain->locate_face(me.x() + distance, me.y() - distance);
        if (southeast == NULL || domain->is_infinite(west))
            zse = (ze + zs) / 2.;
        else
            zse = southeast->get_z();

        southwest = domain->locate_face(me.x() - distance, me.y() - distance);
        if (southwest == NULL || domain->is_infinite(east))
            zsw = (zw + zs) / 2.;
        else
            zsw = southwest->get_z();

        double curve = .25 * ((z - .5 * (zw + ze)) / (2 * distance) + (z - .5 * (zs + zn)) / (2 * distance) +
                              (z - .5 * (zsw + zne)) / (2 * sqrt(2 * distance)) +
                              (z - .5 * (znw + zse)) / (2 * sqrt(2 * distance)));

        face->set_face_data("Liston Curvature",curve);

        if ( fabs(curve) > curmax)
            curmax = fabs(curve);

    }


    for (Delaunay::Finite_faces_iterator fit = domain->finite_faces_begin();
         fit != domain->finite_faces_end(); ++fit)
    {
        double value = fit->face_data("Liston Curvature");
        value =  value / curmax / 2.0;//rescale to [-0.5,+0.5];
        fit->set_face_data("Liston Curvature",value);
    }

}
void Liston_wind::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{

    std::vector< boost::tuple<double, double, double> > values;
    for (auto& s : global_param->stations)
    {
        double u = s->get("u");
        values.push_back( boost::make_tuple(s->x(), s->y(), u ) );
    }

    interp_base* interp=nullptr;
    std::string interp_method = "spline";
    if(interp_method == "spline")
        interp = new thin_plate_spline();

    auto query = boost::make_tuple(elem->get_x(), elem->get_y(), elem->get_z());
    double value = (*interp)(values, query);

    elem->set_face_data("u", value);

    delete interp;
}

Liston_wind::~Liston_wind()
{

}