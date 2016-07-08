#pragma once

#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/Delaunay_mesh_criteria_2.h>
#include <utility>
#include <ostream>
#include <algorithm>
#include <map>
#include <boost/shared_ptr.hpp>
#include <vector>
#include "triangle.h"

    template<class CDT>
    class mesh_2_criterion_area
    {
    protected:
        typedef typename CDT::Geom_traits Geom_traits;
        double max_area;
        double min_area;
        double B;
        Geom_traits traits;
        const std::vector< std::pair< boost::shared_ptr<raster>,double> >& r;
        const std::vector< std::pair< boost::shared_ptr<raster>,double> >& category_rasters;
    public:

        mesh_2_criterion_area(const double aspect_bound = 0.125,
                              const double max_area = 0,
                              const double min_area = 1,
                              const std::vector< std::pair< boost::shared_ptr<raster>,double> >& rasters = std::vector< std::pair< boost::shared_ptr<raster>,double>> (),
                              const std::vector< std::pair< boost::shared_ptr<raster>,double> >& category_rasters = std::vector< std::pair< boost::shared_ptr<raster>,double>>(),
                              const Geom_traits &traits = Geom_traits())
        : r(rasters),category_rasters(category_rasters)

        {
            this->max_area = max_area;
            this->min_area = min_area;
            B = aspect_bound;
            this->traits=traits;

        }

        inline
        double maxarea() const
        { return max_area; }

        inline
        void set_area_bound(const double ab)
        { max_area = ab; }

        inline
        double bound() const
        { return B; }

        inline
        void set_bound(const double bound)
        { B = bound; }


        class Quality
        {
        public:

            Quality()
            {
                _sine = 0;
                _area = 0;
            };
//
//            Quality(double sine, double area, double tolerance)
//            {
//                _sine = sine;
//
//                _area = area;
//
//            }

            const double &area() const
            { return _area; }

            const double &sine() const
            { return _sine; }

            const double &tol() const
            { return _tolerance; }

//        // q1<q2 means q1 is prioritised over q2
//        // ( q1 == *this, q2 == q )
            bool operator<(const Quality &q) const
            {
                if ( this->_area < q._area)
                    return true;
                else
                {
                    if( this->_sine < q._sine)
                        return true;
                    else
                    {
                        if(_tolerance.size() != 0 && q._tolerance.size() == 0)
                            return false; //they're empty, we aren't, we should check their tolerance
                        if(_tolerance.size() == 0 && q._tolerance.size() != 0)
                            return true; //we're empty, they aren't, we should check our tolerance

                        if(_tolerance.size() != 0 && q._tolerance.size() != 0)
                        {
                            //compare who has the the greatest current error
                            auto my_max_tol = std::max_element(_tolerance.begin(), _tolerance.end(),
                                                               [] (std::pair< double,double> const& lhs, std::pair< double,double> const& rhs)
                                                               {return lhs.first < rhs.first;});

                            auto their_max_tol = std::max_element(q._tolerance.begin(), q._tolerance.end(),
                                                                  [] (std::pair< double,double> const& lhs, std::pair< double,double> const& rhs)
                                                                  {return lhs.first < rhs.first;});

                            if(my_max_tol->first > their_max_tol->first )
                                return true; //remember we want to be preferential if our max tolerance is worse
                        }

                        return (sine() < q.sine()); //check angles as tie breaker
                    }
                }
            }

            double _area;
            double _sine;
            std::vector< std::pair<double,double> > _tolerance; // first = tol, seond = max tol

            std::vector<std::pair<double,double>> _category_tol; //first = current percent, second = threshold percent

        };

        class Is_bad
        {

        protected:
            const double B;
            const double max_area;
            const double min_area;
            const std::vector< std::pair< boost::shared_ptr<raster>,double> >& r;
            const std::vector< std::pair< boost::shared_ptr<raster>,double> >& category_rasters;
            const Geom_traits &traits;
        public:

            typedef typename CDT::Point Point_2;

            Is_bad(const double aspect_bound,
                   const double area_bound,
                   const double min_area,
                   const std::vector< std::pair< boost::shared_ptr<raster>,double> >& r,
                   const std::vector< std::pair< boost::shared_ptr<raster>,double> >& category_rasters,
                   const Geom_traits &traits)
                    : B(aspect_bound), max_area(area_bound), min_area(min_area),
                      r(r), category_rasters(category_rasters),traits(traits)
            {

            }

            double elevation_tolerance(const typename CDT::Face_handle &fh, const raster& r) const
            {
                triangle t;
                vertex _v0;
                vertex _v1;
                vertex _v2;

                _v0[0] = fh->vertex(0)->point().x();
                _v0[1] = fh->vertex(0)->point().y();

                _v1[0] = fh->vertex(1)->point().x();
                _v1[1] = fh->vertex(1)->point().y();

                _v2[0] = fh->vertex(2)->point().x();
                _v2[1] = fh->vertex(2)->point().y();

                t.make_rasterized(_v0, _v1, _v2, r);
                if(t.is_nan)
                    return 0; //bail

                auto pxpy = t.rasterized_triangle->xy_to_pxpy(t.v0[0],t.v0[1]);
                t.v0[0] = pxpy.first;
                t.v0[1] = pxpy.second;

                pxpy = t.rasterized_triangle->xy_to_pxpy(t.v1[0],t.v1[1]);
                t.v1[0] = pxpy.first;
                t.v1[1] = pxpy.second;


                pxpy = t.rasterized_triangle->xy_to_pxpy(t.v2[0],t.v2[1]);
                t.v2[0] = pxpy.first;
                t.v2[1] = pxpy.second;


                //Xs or Y s are collinear. This seems to happen if we've set a min triangle area less than the pixel resolution of our dem
                //when this happens, all the points end up collinear in 1 axis, and anything that depends on this
                if(  (t.v0[0] == 0 && t.v1[0] == 0 && t.v2[0] == 0) || (t.v0[1] == 0 && t.v1[1] == 0 && t.v2[1] == 0) )
                {
                    return 0; //bail out nothing we can do
                }

                //create the vectors veco0 and
                double u1,u2,u3;
                double v1,v2,v3;
                double o1,o2,o3; //origin of tri

                o1 = t.v0[0];
                o2 = t.v0[1];
                o3 = t.v0[2];

                //if we have only nan z values, bail.
                if(isnan(t.v0[2]) || isnan(t.v1[2]) || isnan(t.v2[2]))
                {
                    t.is_nan=true;
                    return 0;
                }

                //following http://www.had2know.com/academics/equation-plane-through-3-points.html
                //create the two vectors
                u1 = t.v1[0] - o1;
                u2 = t.v1[1] - o2;
                u3 = t.v1[2] - o3;

                v1 = t.v2[0] - o1;
                v2 = t.v2[1] - o2;
                v3 = t.v2[2] - o3;

                //calculate the normal vector via cross product
                double a, b, c;
                a = u2 * v3 - v2 * u3;
                b = v1 * u3 - u1 * v3;
                c = u1 * v2 - v1 * u2;


                //solve for d
                double d =a*o1 + b*o2 + c*o3;

                double rmse = 0;
                double n = 0;

                std::map<int,int> lc;
                for (int y = 0; y < t.rasterized_triangle->getDs()->GetRasterYSize(); y++)
                {
                    for (int x = 0; x < t.rasterized_triangle->getDs()->GetRasterXSize(); x++)
                    {
                        double value = t.rasterized_triangle->getpXpY(x, y);
                        if (!isnan(value))
                        {
                            double z = -(a*x+b*y-d)/c; //plane eqn solved for z. allows us to predict z values via x,y coords
                            double diff = z - value;

                            rmse += diff * diff;
                            n++;
                        }
                    }
                }

                //bail, somehow we have no raster cells under our triangle.
                if (n == 0.)
                    return 0;

                rmse /= n;

                rmse = sqrt(rmse);

                return rmse;

            }
            double categoryraster_isok(const typename CDT::Face_handle &fh,const raster& r) const
            {
                triangle t;
                vertex _v0;
                vertex _v1;
                vertex _v2;

                _v0[0] = fh->vertex(0)->point().x();
                _v0[1] = fh->vertex(0)->point().y();

                _v1[0] = fh->vertex(1)->point().x();
                _v1[1] = fh->vertex(1)->point().y();

                _v2[0] = fh->vertex(2)->point().x();
                _v2[1] = fh->vertex(2)->point().y();

                t.make_rasterized(_v0, _v1, _v2, r);


                if(t.is_nan)
                    return true; //bail

                if( t.v0[2] != t.v1[2] ||
                    t.v0[2] != t.v2[2] ||
                    t.v1[2] != t.v2[2])
                    return false;

                double n=0;//total tri
                std::map<int,int> lc; //holds the count of each category (e.g., landcover type)
                for (int y = 0; y < t.rasterized_triangle->getDs()->GetRasterYSize(); y++)
                {
                    for (int x = 0; x < t.rasterized_triangle->getDs()->GetRasterXSize(); x++)
                    {
                        double value = t.rasterized_triangle->getpXpY(x, y);
                        if (!isnan(value))
                        {
                            int key = (int)value;
                            if(lc.find(key) != lc.end())
                                lc[key]++;
                            else
                                lc[key] = 1;
                            n++;
                        }
                    }
                }

                if (n==0)
                    return true; //bail

                bool isok=false;

                std::vector<double> percents;

                for(auto& itr : lc)
                {
                    double pCover =  itr.second / n;
                    percents.push_back(pCover);

                }
                return *std::max_element(percents.begin(),percents.end()); // the most dominate landcover

            }
            double mean_elevation_diff(const typename CDT::Face_handle &fh,const raster & r) const
            {
                triangle t;
                vertex v0;
                vertex v1;
                vertex v2;

                v0[0] = fh->vertex(0)->point().x();
                v0[1] = fh->vertex(0)->point().y();

                v1[0] = fh->vertex(1)->point().x();
                v1[1] = fh->vertex(1)->point().y();

                v2[0] = fh->vertex(2)->point().x();
                v2[1] = fh->vertex(2)->point().y();

                t.make_rasterized(v0, v1, v2, r);
                if(t.is_nan)
                    return 0; //bail

                // Initialize triangle mean elevation (m)
                double triangle_z_mean = 0;
                // Initialize count of triangle vertices found not-nan
                double tri_count = 0;

                // Check if elevations of vertex is not-nan
                if (!isnan(t.v0[2]))
                {
                    triangle_z_mean += t.v0[2];
                    ++tri_count;
                }

                if (!isnan(t.v1[2]))
                {
                    triangle_z_mean += t.v1[2];
                    ++tri_count;
                }

                if (!isnan(t.v2[2]))
                {
                    triangle_z_mean += t.v2[2];
                    ++tri_count;
                }


                //planar mean elevation
                triangle_z_mean /= tri_count;

                // If no non-nan vertices were found, return zero
                if (tri_count == 0)
                    return 0;

                // Initialize sum and count of grid cell elevations
                double sum = 0;
                double count = 0;

                for (int i = 0; i < t.rasterized_triangle->getDs()->GetRasterYSize(); i++)
                {
                    for (int j = 0; j < t.rasterized_triangle->getDs()->GetRasterXSize(); j++)
                    {
                        double value = t.rasterized_triangle->getpXpY(j, i);
                        if (!isnan(value))
                        {
                            sum += value;
                            ++count;
                        }
                    }
                }

                // If none were found return zero
                if (count == 0)
                    return 0;

                // Take mean of all grid cell elevations found
                double mean = sum / count;

                // Take difference of means
                double diff = fabs(triangle_z_mean - mean);

                return diff;

            }

            CGAL::Mesh_2::Face_badness operator()(const Quality q) const
            {
                if (q.area() > max_area)
                    return CGAL::Mesh_2::IMPERATIVELY_BAD; //IMPERATIVELY_BAD

                if (q.sine() < this->B)
                    return CGAL::Mesh_2::BAD;

                if (q.area() <= min_area )
                    return CGAL::Mesh_2::NOT_BAD;

                for(auto& itr : q._tolerance)
                {
                    //first == current tol
                    //second == max tol
                    if(itr.first > itr.second) //do we violate the max tol
                        return CGAL::Mesh_2::BAD;
                }

                for(auto& itr : q._category_tol)
                {
                    //first == current tol
                    //second == max tol
                    if(itr.first < itr.second)  // is our max land cover below the required threshold to keep this triangle?
                        return CGAL::Mesh_2::BAD;
                }

                return CGAL::Mesh_2::NOT_BAD;
            }

            CGAL::Mesh_2::Face_badness operator()(const typename CDT::Face_handle &fh,
                                                  Quality &q) const
            {
                typedef typename CDT::Geom_traits Geom_traits;
                typedef typename Geom_traits::Compute_area_2 Compute_area_2;
                typedef typename Geom_traits::Compute_squared_distance_2
                        Compute_squared_distance_2;

                Geom_traits traits;

                Compute_area_2 area_2 = traits.compute_area_2_object();
                Compute_squared_distance_2 squared_distance = traits.compute_squared_distance_2_object();

                const Point_2 &pa = fh->vertex(0)->point();
                const Point_2 &pb = fh->vertex(1)->point();
                const Point_2 &pc = fh->vertex(2)->point();

                double a = CGAL::to_double(squared_distance(pb, pc));
                double b = CGAL::to_double(squared_distance(pc, pa));
                double c = CGAL::to_double(squared_distance(pa, pb));

                double area = CGAL::to_double(area_2(pa, pb, pc));

                // area = 4 * area^2(triangle)
                double area2 = 2 * CGAL::to_double(area_2(pa, pb, pc));
                area2 = area2 * area2;

                if (a < b) if (a < c)
                    q._sine = area2 / (b * c);
                else
                    q._sine = area2 / (a * b);
                else if (b < c)
                    q._sine = area2 / (a * c);
                else
                    q._sine = area2 / (a * b);

                if (max_area != 0)
                {
                    q._area = area;
                }

                //if our triangle already fails either: area or angle, bail now.
                auto current_badness = operator()(q);
                if (current_badness != CGAL::Mesh_2::NOT_BAD)
                    return current_badness;

                //do numeric rasters
                for(auto& itr : r)
                {
                    auto t = elevation_tolerance(fh,*(itr.first));
                    q._tolerance.push_back(std::make_pair(t,itr.second));

                    current_badness = operator()(q);
                    if (current_badness != CGAL::Mesh_2::NOT_BAD)
                        return current_badness;
                }

                for(auto& itr : category_rasters)
                {
                    auto t = categoryraster_isok(fh,*(itr.first));
                    q._category_tol.push_back(std::make_pair(t,itr.second));

                    current_badness = operator()(q);
                    if (current_badness != CGAL::Mesh_2::NOT_BAD)
                        return current_badness;
                }
                return operator()(q);
            }
        };

        Is_bad is_bad_object() const
        {
            return Is_bad(this->bound(), max_area, min_area, r, category_rasters, this->traits);
        }
    };
