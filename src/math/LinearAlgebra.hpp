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

#pragma once

#define HAVE_TEUCHOS_DEBUG
#define TEUCHOS_DEBUG

#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "triangulation.hpp"

namespace math
{
    namespace LinearAlgebra
    {

      // Typedefs/aliases for ease of Trilinos use
      using Teuchos::arcp;
      using Teuchos::ArrayRCP;
      using Teuchos::Comm;
      using Teuchos::ParameterList;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::Time;
      using Teuchos::tuple;
      typedef Tpetra::CrsGraph<> graph_type;
      typedef Tpetra::CrsMatrix<> crs_matrix_type;
      typedef Tpetra::Map<> map_type;
      typedef Tpetra::MultiVector<> MV;
      typedef Tpetra::Operator<> OP;
      typedef Tpetra::RowMatrix<> row_matrix_type;
      typedef MV::scalar_type scalar_type;
      typedef Ifpack2::Preconditioner<> prec_type;
      typedef Belos::LinearProblem<scalar_type, MV, OP> problem_type;
      typedef Belos::SolverManager<scalar_type, MV, OP> solver_type;
      typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;

      typedef Tpetra::Details::DefaultTypes::global_ordinal_type global_ordinal_type;
      typedef Tpetra::Details::DefaultTypes::local_ordinal_type local_ordinal_type;

      struct SolveConverge
      {
	int numIters;
	double residual;
      };

      class NearestNeighborProblem
      {

	// Unstructured triangular mesh with nLayer extrusion
	mesh& m_domain;
	int m_nLayer;    // defaults to 1 in constructor

	// Trilinos structures for linear system variables
	RCP<const Comm<int>> m_comm;
	RCP<crs_matrix_type> m_matrix;
	RCP<MV> m_solution, m_rhs;

	// Trilinos structures for linear system sparsity
	RCP<const map_type> m_map;
	RCP<graph_type> m_graph;

	// Trilinos structures for linear system solve
	// - needed for factories
	RCP<solver_type> m_solver;
	RCP<prec_type> m_preconditioner;
	RCP<problem_type> m_problem;

      public:
	NearestNeighborProblem(mesh& domain, int nLayer=1);
	~NearestNeighborProblem();

	void zeroSystem();

	void matrixReplaceGlobalValues(global_ordinal_type global_row_idx, global_ordinal_type global_col_idx, double val);
	void matrixSumIntoGlobalValues(global_ordinal_type global_row_idx, global_ordinal_type global_col_idx, double val);

	void matrixResumeFill();
	void matrixFillComplete();

	void rhsSumIntoGlobalValue(global_ordinal_type global_idx, double val);

	SolveConverge Solve();

	double getSolutionMax();
	double getRhsMax();

	ArrayRCP<const double> getSolutionView();

	// Dumping the problem and solution to MatrixMarket format for inspection
	void writeSystemMatrixMarket(std::string file_prefix);
	void writeSolutionMatrixMarket(std::string file_prefix);
      };

    }
}
