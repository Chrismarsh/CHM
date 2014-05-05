#include "spline.hpp"

double spline::operator()(station_list const&  stations, mesh_elem& elem, std::string variable, boost::shared_ptr<interp_visitor> visitor, boost::shared_ptr<global> global_param)
{

//    unsigned int size = 0; //number of points
//    size++; // need to make room for the constants
//
//    double *A = new double[ size * size ];
//    double pi = 3.14159;
//    double c = 0.577215; //euler constant
//    double *B = new double[size]; // known values - constant value of 0 goes in b[size-1]
//
//
////    double dij = sqrt( pow((int)(points->at(i).col) - (int)(points->at(j).col),2.0) + pow((int)(points->at(i).row) - (int)(points->at(j).row),2.0) );
//    double dij = 1; // TODO Fix this hardcoded distance value
//    double Rd = 0.0;
//
//    if ( (dij - 0.0) < 0.000001)
//        Rd = 0.0;
//    else
//        Rd = -1.0/(2.0*pi*0.01)*( log(dij*0.1/2.0) + c + gsl_sf_bessel_K0(dij*0.1));
//
//    //A(i,j+1) = Rd;
//    //A[i * size + (j+1)] = Rd;
//
//    //A(j,i+1) = Rd;
//    //A[j * size + (i+1)] = Rd;			
//
//
//   //set constants and build b values
//   for(unsigned int i = 0; i< size;i++)
//   {
//           A[i*size+0] = 1;
//           A[(size-1)*size+i] = 1;
//   }
//
//   for(unsigned int i = 0; i<size -1; i++)
//   {
//           B[i] = 0;//points->at(i).z;
//   }
//   B[size-1] = 0.0; //constant
//
//
//   A[(size-1)*(size)] = 0.0;
//
//   //solve equation via gsl via lu decomp
//   gsl_matrix_view m 	= gsl_matrix_view_array (A, size, size);
//   gsl_vector_view b	= gsl_vector_view_array (B, size);
//   gsl_vector *x		= gsl_vector_alloc (size);
//   int s				= 0;
//   gsl_permutation * p = gsl_permutation_alloc (size);
//   gsl_linalg_LU_decomp (&m.matrix, p, &s);
//   gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
//
//   //location of the unknown
//   unsigned int ux = 0;//col;
//   unsigned int uy = 0;//row;
//
//   double z0 = x->data[0];//little a
//
//   //-1 to skip the x[0] we already pulled off above
//   for (unsigned int i = 0;i < (x->size - 1) ;i++)
//   {
//           //TODO fix this
//           double dxy = 0;//sqrt( pow((int)(ux) - (int)(points->at(i).col),2.0) + pow((int)(uy) - (int)(points->at(i).row),2.0) );
//           double Rd = 0.0;
//
//           if ((dxy - 0.0) < 0.000001)
//                   Rd = 0.0;
//           else
//                   Rd = -1.0/(2.0*pi*0.01)*( log(dxy*0.1/2.0) + c + gsl_sf_bessel_K0(dxy*0.1));
//
//           z0 = z0 + x->data[i+1]*Rd; 
//   }
//
//
//   gsl_permutation_free (p);
//   gsl_vector_free (x);
//   delete B;
//   delete A;
//
//   //loss of precision
//   float fz0 = (float)z0;

//   if(visitor)
//   {
//           visitor->operator ()(&fz0,row,col);
//   }

//   raster->set_raster_at(row,col,fz0);

   return 0;
}

spline::spline()
{
    
}

spline::~spline()
{
    
}