//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#include "LinearAlgebra.hpp"

namespace math
{
  namespace LinearAlgebra
  {

    NearestNeighborProblem::NearestNeighborProblem(mesh& domain, int nLayer) :
      m_domain(domain), m_nLayer(nLayer)
    {

      // TODO Accept a communicator on construction
      m_comm = Tpetra::getDefaultComm();

      // Sizes of the domain
      size_t ntri = domain->size_faces();
      size_t n_global_tri = domain->size_global_faces();

      /*
	Initialize the local-global index map for the extruded mesh system.
      */
      auto global_IDs = domain->get_global_IDs();
      std::vector<int> extruded_global_IDs(ntri*nLayer);
      // Create the global IDs for the extruded system
      // Ordering:
      // - mesh elements and then layers successively
      auto extruded_ID_iterator = extruded_global_IDs.begin();
      for (int i=0; i<nLayer; ++i) {
	std::transform(global_IDs.begin(),global_IDs.end(),extruded_ID_iterator,
		       [=](int id) -> int { return i*n_global_tri + id; } );
	extruded_ID_iterator += ntri;
      }



      const int numGlobalElements = n_global_tri*nLayer;
      const int* data_extruded_IDs = extruded_global_IDs.data();
      const int indexListSize = ntri*nLayer;
      int indexBase = 0;

      m_map = rcp(new map_type(numGlobalElements, data_extruded_IDs, indexListSize, indexBase, m_comm));

      // loop over locally owned rows, figure out number of neighbors (owned or
      // otherwise!), and what their global indices are.
      std::vector<size_t> num_entries(ntri*nLayer,1);
      std::vector<std::array<int,6>> neighbor_global_idx(ntri*nLayer);
#pragma omp parallel for
      for (size_t i = 0; i < ntri; ++i) {
	auto face = domain->face(i);
	int face_bottom_idx = face->cell_global_id;
	int face_bottom_local_idx = face->cell_local_id;
	// Lateral neighbors and self
	for (int layer=0; layer<nLayer; ++layer ) {
	  int element_idx = n_global_tri*layer + face_bottom_idx;
	  int local_array_idx = ntri*layer + face_bottom_local_idx;
	  neighbor_global_idx[local_array_idx][0] = element_idx;
	  for (int f = 0; f < 3; f++){
	    auto neighbor = face->neighbor(f);
	    if (neighbor != nullptr)  {
	      int neigh_bottom_idx = neighbor->cell_global_id;
	      int neigh_global_idx = n_global_tri*layer + neigh_bottom_idx;
	      neighbor_global_idx[local_array_idx][num_entries[local_array_idx]] = neigh_global_idx;
	      ++num_entries[local_array_idx];
	    }
	  }
	}
	/*
	  Above and below neighbor loops are null when single layer
	*/
	// Neighbor below
	for (int layer=1; layer<nLayer; ++layer ) {
	  int element_idx = n_global_tri*layer + face_bottom_idx;
	  int local_array_idx = ntri*layer + face_bottom_local_idx;
	  int below_idx = n_global_tri*(layer-1) + face_bottom_idx;
	  neighbor_global_idx[local_array_idx][num_entries[local_array_idx]] = below_idx;
	  ++num_entries[local_array_idx];
	}
	// Neighbor above
	for (int layer=0; layer<nLayer-1; ++layer ) {
	  int element_idx = n_global_tri*layer + face_bottom_idx;
	  int local_array_idx = ntri*layer + face_bottom_local_idx;
	  int above_idx = n_global_tri*(layer+1) + face_bottom_idx;
	  neighbor_global_idx[local_array_idx][num_entries[local_array_idx]] = above_idx;
	  ++num_entries[local_array_idx];
	}
      }

      /*
	Set up the CrsGraph structure for creating the distributed CrsMatrix and
	Vectors for the linear nearest neighbor system

	Graph for nearest neighbor connectivity.
	6 entries (max) per row: self, three neighbors, above and below
	(fewer entries for boundary elements)

	Preferred construction of CrsGraph uses a Teuchos::ArrayView<T> for num entries/row
      */
      Teuchos::ArrayView<size_t> num_entries_view(num_entries.data(),ntri*nLayer);
      m_graph = rcp (new graph_type (m_map, num_entries_view, Tpetra::StaticProfile));

      // Set all of the desired nonzero columns in the graph
      // due to insertGlobalIndices args
      // DO NOT DO THIS THREAD PARALLEL
      for (size_t i = 0; i < ntri; ++i) {
	auto face = domain->face(i);
	int face_bottom_idx = face->cell_global_id;
	int face_bottom_local_idx = face->cell_local_id;
	for (int layer = 0; layer<nLayer; ++layer) {
	  int element_idx = n_global_tri*layer + face_bottom_idx;
	  int local_array_idx = ntri*layer + face_bottom_local_idx;
	  // std::cout << "insertGlobal : " << element_idx << " : " << num_suspension_entries[element_idx] << "\n";
	  m_graph->insertGlobalIndices(element_idx, num_entries[local_array_idx], &(neighbor_global_idx.data()[local_array_idx][0]));
	}
      }
      m_graph->fillComplete();

      // Create a Tpetra::Matrix using the Map, with a static allocation
      // dictated by NumNz.  (We know exactly how many elements there will
      // be in each row, so we use static profile for efficiency.)
      m_matrix = rcp (new crs_matrix_type (m_graph));
      m_matrix->fillComplete();
      m_rhs = rcp(new MV(m_map, 1));
      m_solution = rcp(new MV(m_map, m_rhs->getNumVectors()));

      /*
	Belos solver and Ifpack2 preconditioner setup
	- TODO add optional input to set these params
      */
      // Create Belos iterative linear solver.
      RCP<ParameterList> solverParams(new ParameterList()); // solve parameters go in here
      solverParams->set( "Block Size", 1 );
      solverParams->set( "Num Blocks", 30 );
      solverParams->set( "Maximum Iterations", 1000 );
      solverParams->set( "Convergence Tolerance", 1e-8 );
      {
        Belos::SolverFactory<scalar_type, MV, OP> belosFactory;
        m_solver = belosFactory.create("GMRES", solverParams);
      }
      if (m_solver.is_null())
	{
	  BOOST_THROW_EXCEPTION(module_error() << errstr_info("PBSM3D failed to create solver"));
	}

      m_preconditioner = Ifpack2::Factory::create<row_matrix_type>("ILUT", m_matrix);
      if (m_preconditioner.is_null())
	{
	  BOOST_THROW_EXCEPTION(module_error() << errstr_info("PBSM3D failed to create preconditioner"));
	}
      ParameterList precondOptions;
      precondOptions.set("fact: drop tolerance", 1e-4);
      precondOptions.set("fact: ilut level-of-fill", 3.0); // Note this is different from num_entries_per_row: https://docs.trilinos.org/dev/packages/ifpack2/doc/html/classIfpack2_1_1ILUT.html#aee2011b313e3070ee43b2cfc2d183634
      m_preconditioner->setParameters(precondOptions);
      m_preconditioner->initialize();

      // Specify the deposition problem
      m_problem = rcp (new problem_type (m_matrix, m_solution, m_rhs));
      if (! m_preconditioner.is_null ()) {
	m_problem->setRightPrec (m_preconditioner);
      }
      m_problem->setProblem ();
      m_solver->setProblem (m_problem);

    } // end constructor

    void NearestNeighborProblem::zeroSystem()
    {
      m_matrix->resumeFill();
      // Zero out suspension system
      m_matrix->setAllToScalar(0.0);
      m_rhs->putScalar(0.0);
      m_solution->putScalar(0.0);
    }

    void NearestNeighborProblem::matrixReplaceGlobalValues(int global_row_idx, int global_col_idx, double val)
    {
      m_matrix->replaceGlobalValues(global_row_idx, tuple(global_col_idx), tuple(val));
    }

    void NearestNeighborProblem::matrixSumIntoGlobalValues(int global_row_idx, int global_col_idx, double val)
    {
      m_matrix->sumIntoGlobalValues(global_row_idx, tuple(global_col_idx), tuple(val));
    }

    void NearestNeighborProblem::rhsSumIntoGlobalValue(int global_idx, double val)
    {
      m_rhs->sumIntoGlobalValue(global_idx, 0, val);
    }

    SolveConverge NearestNeighborProblem::Solve()
    {
      m_matrix->fillComplete();
      m_preconditioner->compute();

      // Solve the linear system.
      m_solver->reset(Belos::Problem);
      {        Belos::ReturnType solveResult = m_solver->solve();
        if (solveResult != Belos::Converged)
	  {
            if (m_comm->getRank() == 0)
	      {
		BOOST_THROW_EXCEPTION(module_error() << errstr_info("Belos solver failed to converge"));
	      }
            // return EXIT_FAILURE;
	  }
      }

      // Get (and return) convergence info
      SolveConverge tmp;
      tmp.numIters = m_solver->getNumIters();
      tmp.residual = m_solver->achievedTol();
      return tmp;

    }

    // Solution's maximum value can be computed by InfNorm
    double NearestNeighborProblem::getSolutionMax()
    {
      double tempStorage;
      Teuchos::ArrayView<double> max_value(&tempStorage,1);
      m_solution->normInf(max_value);
      return max_value[0];
    }

    // RHS's maximum value can be computed by InfNorm
    double NearestNeighborProblem::getRhsMax()
    {
      double tempStorage;
      Teuchos::ArrayView<double> max_value(&tempStorage,1);
      m_rhs->normInf(max_value);
      return max_value[0];
    }

    ArrayRCP<const double> NearestNeighborProblem::getSolutionView()
    {
      return m_solution->get1dView();
    }

    void NearestNeighborProblem::writeSystemMatrixMarket(std::string file_prefix)
    {
      std::string matrix_file = file_prefix + "_matrix.mm";
      std::string matrix_name = file_prefix + " system matrix";
    Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparseFile( matrix_file, m_matrix, matrix_name, matrix_name );

    std::string rhs_file = file_prefix + "_rhs.mm";;
    std::string rhs_name= file_prefix + " system rhs";
    Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeDenseFile( rhs_file, m_rhs, rhs_name, rhs_name );

    }

    void NearestNeighborProblem::writeSolutionMatrixMarket(std::string file_prefix)
    {
      std::string solution_file = file_prefix + "_solution.mm";
      std::string solution_name = file_prefix + " system solution";
      Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeDenseFile( solution_file, m_solution, solution_name, solution_name );
    }

    NearestNeighborProblem::~NearestNeighborProblem(){}

  }
}
