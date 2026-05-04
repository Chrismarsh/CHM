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
#include <Tpetra_Core.hpp>

namespace math
{
    namespace LinearAlgebra
    {
        size_t Sizes::total_elements() const { return global * vert_layers; }
        size_t Sizes::local_elements() const { return local * vert_layers; }

        void NearestNeighborProblem::set_comm_for_parallel()
        {
                m_comm = Tpetra::getDefaultComm();
        };

        void NearestNeighborProblem::build_Belos_solver()
        {
            /*
          Belos solver and Ifpack2 preconditioner setup
          - TODO add optional input to set these params
            */
            // Create Belos iterative linear solver.
            RCP<ParameterList> solverParams(new ParameterList()); // solve parameters go in here
            solverParams->set("Block Size", 1);
            solverParams->set("Num Blocks", 30);
            solverParams->set("Maximum Iterations", 1000);
            solverParams->set("Convergence Tolerance", 1e-8);
            {
                Belos::SolverFactory<scalar_type, MV, OP> belosFactory;
                m_solver = belosFactory.create("GMRES", solverParams);
            }
            if (m_solver.is_null())
            {
                CHM_THROW_EXCEPTION(module_error, "PBSM3D failed to create solver");
            }
        }

        void NearestNeighborProblem::zeroSystem()
        {
            m_matrix->resumeFill();
            // Zero out suspension system
            m_matrix->setAllToScalar(0.0);
            m_rhs->putScalar(0.0);
            m_solution->putScalar(0.0);
        }

        void NearestNeighborProblem::matrixReplaceGlobalValues(global_ordinal_type global_row_idx,
                                                               global_ordinal_type global_col_idx, double val)
        {
            m_matrix->replaceGlobalValues(global_row_idx, tuple(global_col_idx), tuple(val));
        }

        void NearestNeighborProblem::matrixSumIntoGlobalValues(global_ordinal_type global_row_idx,
                                                               global_ordinal_type global_col_idx, double val)
        {
            m_matrix->sumIntoGlobalValues(global_row_idx, tuple(global_col_idx), tuple(val));
        }

        void NearestNeighborProblem::rhsSumIntoGlobalValue(global_ordinal_type global_idx, double val)
        {
            // Critical section needed because sumIntoGlobalValues is not respecting
            // the 4th arg (force atomic update)
            // - likely a Trilinos/Kokkos bug
            // trilinos/Trilinos/issues/9519
#pragma omp critical
            m_rhs->sumIntoGlobalValue(global_idx, 0, val, true);
        }

        SolveConverge NearestNeighborProblem::Solve()
        {
            m_matrix->fillComplete();
            m_preconditioner->compute();

            // Solve the linear system.
            m_solver->reset(Belos::Problem);
            {
                Belos::ReturnType solveResult = m_solver->solve();
                if (solveResult != Belos::Converged)
                {
                    if (m_comm->getRank() == 0)
                    {
                        CHM_THROW_EXCEPTION(module_error, "Belos solver failed to converge");
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
            Teuchos::ArrayView<double> max_value(&tempStorage, 1);
            m_solution->normInf(max_value);
            return max_value[0];
        }

        // RHS's maximum value can be computed by InfNorm
        double NearestNeighborProblem::getRhsMax()
        {
            double tempStorage;
            Teuchos::ArrayView<double> max_value(&tempStorage, 1);
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
            Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparseFile(
                matrix_file, m_matrix, matrix_name, matrix_name);

            std::string rhs_file = file_prefix + "_rhs.mm";;
            std::string rhs_name = file_prefix + " system rhs";
            Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeDenseFile(rhs_file, m_rhs, rhs_name, rhs_name);
        }

        void NearestNeighborProblem::writeSolutionMatrixMarket(std::string file_prefix)
        {
            std::string solution_file = file_prefix + "_solution.mm";
            std::string solution_name = file_prefix + " system solution";
            Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeDenseFile(
                solution_file, m_solution, solution_name, solution_name);
        }

        NearestNeighborProblem::~NearestNeighborProblem()
        {
        }

    } // namespace LinearAlgebra
}
