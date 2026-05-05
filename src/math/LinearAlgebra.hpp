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

#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <concepts>
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
        using graph_type = Tpetra::CrsGraph<>;
        using crs_matrix_type = Tpetra::CrsMatrix<>;
        using map_type = Tpetra::Map<>;
        using MV = Tpetra::MultiVector<>;
        using OP = Tpetra::Operator<>;
        using row_matrix_type = Tpetra::RowMatrix<>;
        using scalar_type = MV::scalar_type;
        using prec_type = Ifpack2::Preconditioner<>;
        using problem_type = Belos::LinearProblem<scalar_type, MV, OP>;
        using solver_type = Belos::SolverManager<scalar_type, MV, OP>;
        using reader_type = Tpetra::MatrixMarket::Reader<crs_matrix_type>;
        // The type used for indexes on the local rank
        using global_index_type = Tpetra::Map<>::global_ordinal_type;
        using local_index_type = Tpetra::Map<>::local_ordinal_type;

        struct IndexTracker;

        struct SolveConverge
        {
            int numIters;
            double residual;
        };
		
        template<typename T>
        concept MeshObject = requires(T& t,size_t i)
        {
                {t->size_local_faces()} -> std::convertible_to<size_t>;
                {t->size_global_faces()} -> std::convertible_to<size_t>;
                {t->get_global_IDs()} -> std::ranges::range;
            requires std::convertible_to<
                        std::ranges::range_value_t<decltype(t->get_global_IDs())>,
                        int
                        >;
                {t->face(i)} -> std::same_as<mesh_elem>;
        };

        struct Sizes
        {
            size_t local;
            size_t global;
            size_t vert_layers;
            size_t total_elements() const;
            size_t local_elements() const;
        };

        class NearestNeighborProblem
        {
            void set_comm_for_parallel();
            template<MeshObject M> static Sizes set_domain_sizes(M&,size_t N);

            template<MeshObject M>
            RCP<const map_type> init_local_global_index_map(M&,const Sizes&);
            template <MeshObject M> RCP<graph_type> init_m_graph(M& domain,IndexTracker&, Sizes& sizes);
            template <MeshObject M> void set_nonzero_elements(M& domain, const IndexTracker&, const Sizes& sizes);
            template <MeshObject M> void fill_matrix_rhs_solution();
            template <MeshObject M> void set_problem();

            void build_Belos_solver();

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
			template<MeshObject M>
            NearestNeighborProblem(M& domain, int nLayer = 1);
            ~NearestNeighborProblem();

            void zeroSystem();

            void matrixReplaceGlobalValues(global_index_type global_row_idx, global_index_type global_col_idx,
                                           double val);
            void matrixSumIntoGlobalValues(global_index_type global_row_idx, global_index_type global_col_idx,
                                           double val);

            void matrixResumeFill();
            void matrixFillComplete();

            void rhsSumIntoGlobalValue(global_index_type global_idx, double val);

            SolveConverge Solve();

            double getSolutionMax();
            double getRhsMax();

            ArrayRCP<const double> getSolutionView();

            // Dumping the problem and solution to MatrixMarket format for inspection
            void writeSystemMatrixMarket(std::string file_prefix);
            void writeSolutionMatrixMarket(std::string file_prefix);
        };
		
        template<MeshObject M>
        Sizes NearestNeighborProblem::set_domain_sizes(M& domain,size_t numLayer)
        {
            Sizes s;
            s.local = domain->size_local_faces();
            s.global = domain->size_global_faces();
            s.vert_layers = numLayer;
            return s;
        }

        template <MeshObject M>
        RCP<const map_type> NearestNeighborProblem::init_local_global_index_map(M& domain, const Sizes& sizes)
        {
	    auto global_IDs = domain->get_global_IDs();
            std::vector<global_index_type> extruded_global_IDs(sizes.local_elements());
            // Create the global IDs for the extruded system
            // Ordering:
            // - mesh elements and then layers successively
            auto extruded_ID_iterator = extruded_global_IDs.begin();
            for (int i = 0; i < sizes.vert_layers; ++i)
            {
                std::transform(global_IDs.begin(), global_IDs.end(), extruded_ID_iterator,
                               [=](int id) -> global_index_type { return i * sizes.global + id; });
                extruded_ID_iterator += sizes.local;
            }


            const size_t numGlobalElements = sizes.global * sizes.vert_layers;
            const size_t indexListSize = sizes.local * sizes.vert_layers;
            int indexBase = 0;

            return rcp(new map_type(numGlobalElements, extruded_global_IDs.data(), indexListSize, indexBase, m_comm));
	};

        struct IndexTracker
        {
            std::vector<size_t> entries_per_element;
            std::vector<std::array<global_index_type,6>> neighbor_global_idx;

            template <class M> IndexTracker(M& domain, Sizes& sizes);
        };

        template <class M> IndexTracker::IndexTracker(M& domain, Sizes& sizes)
        {
#pragma omp parallel for
            for (size_t i = 0; i < sizes.local; ++i)
            {
                auto face = domain->face(i);
                int face_bottom_idx = face->cell_global_id;
                int face_bottom_local_idx = face->cell_local_id;

                // Lateral neighbors and self
                for (int layer = 0; layer < sizes.vert_layers; ++layer)
                {
                    int element_idx = sizes.global * layer + face_bottom_idx;
                    int local_array_idx = sizes.local * layer + face_bottom_local_idx;
                    neighbor_global_idx.at(local_array_idx).at(0) = element_idx;

                    for (int f = 0; f < 3; f++)
                    {
                        auto neighbor = face->neighbor(f);

                        if (neighbor != nullptr)
                        {
                            int neigh_bottom_idx = neighbor->cell_global_id;
                            int neigh_global_idx = sizes.global * layer + neigh_bottom_idx;
                            neighbor_global_idx.at(local_array_idx).at(entries_per_element.at(local_array_idx)) =
                                neigh_global_idx;
                            ++entries_per_element.at(local_array_idx);
                        }
                    }
                }
                /*
                  Above and below neighbor loops are null when single layer
                */
                // Neighbor below
                for (int layer = 1; layer < sizes.vert_layers; ++layer)
                {
                    int element_idx = sizes.global * layer + face_bottom_idx;
                    int local_array_idx = sizes.local * layer + face_bottom_local_idx;
                    int below_idx = sizes.global * (layer - 1) + face_bottom_idx;
                    neighbor_global_idx.at(local_array_idx).at(entries_per_element.at(local_array_idx)) = below_idx;
                    ++entries_per_element.at(local_array_idx);
                }
                // Neighbor above
                for (size_t layer = 0; layer < sizes.vert_layers - 1; ++layer)
                {
                    size_t element_idx = sizes.global * layer + face_bottom_idx;
                    size_t local_array_idx = sizes.local * layer + face_bottom_local_idx;
                    size_t above_idx = sizes.global * (layer + 1) + face_bottom_idx;
                    neighbor_global_idx.at(local_array_idx).at(entries_per_element.at(local_array_idx)) = above_idx;
                    ++entries_per_element.at(local_array_idx);
                }
            }
        }
        template <MeshObject M>
	RCP<graph_type> NearestNeighborProblem::init_m_graph(M& domain,IndexTracker& idx_tracker, Sizes& sizes)
	{

            /*
          Set up the CrsGraph structure for creating the distributed CrsMatrix and
          Vectors for the linear nearest neighbor system

          Graph for nearest neighbor connectivity.
          6 entries (max) per row: self, three neighbors, above and below
          (fewer entries for boundary elements)

          Preferred construction of CrsGraph uses a Teuchos::ArrayView<T> for num entries/row
            */
            Teuchos::ArrayView<size_t> num_entries_view(idx_tracker.entries_per_element.data(), sizes.local * sizes.vert_layers);
	    return rcp(new graph_type(m_map,num_entries_view));
	}

        template <MeshObject M> void NearestNeighborProblem::set_nonzero_elements(M& domain, const IndexTracker& idx_tracker, const Sizes& sizes)
        {
            // DO NOT DO THIS THREAD PARALLEL
            for (size_t i = 0; i < sizes.local; ++i)
            {
                auto face = domain->face(i);
                int face_bottom_idx = face->cell_global_id;
                int face_bottom_local_idx = face->cell_local_id;

                for (int layer = 0; layer < sizes.vert_layers; ++layer)
                {
                    int element_idx = sizes.global * layer + face_bottom_idx;
                    int local_array_idx = sizes.local * layer + face_bottom_local_idx;
                    // std::cout << "insertGlobal : " << element_idx << " : " << num_suspension_entries[element_idx] <<
                    // "\n";
                    m_graph->insertGlobalIndices(element_idx, idx_tracker.entries_per_element.at(local_array_idx),
                                                 &(idx_tracker.neighbor_global_idx.at(local_array_idx)[0]));
                }
            }
            m_graph->fillComplete();
        }

        template <MeshObject M> void NearestNeighborProblem::fill_matrix_rhs_solution()
        {
            // Create a Tpetra::Matrix using the Map, with a static allocation
            // dictated by NumNz.  (We know exactly how many elements there will
            // be in each row, so we use static profile for efficiency.)
            m_matrix = rcp(new crs_matrix_type(m_graph));
            m_matrix->fillComplete();
            m_rhs = rcp(new MV(m_map, 1));
            m_solution = rcp(new MV(m_map, m_rhs->getNumVectors()));
        }

        template <MeshObject M> void NearestNeighborProblem::set_problem()
        {
            m_preconditioner = Ifpack2::Factory::create<row_matrix_type>("ILUT", m_matrix);
            if (m_preconditioner.is_null())
            {
                CHM_THROW_EXCEPTION(module_error, "PBSM3D failed to create preconditioner");
            }
            ParameterList precondOptions;
            precondOptions.set("fact: drop tolerance", 1e-4);
            precondOptions.set("fact: ilut level-of-fill", 3.0);
            // Note this is different from num_entries_per_row:
            // https://docs.trilinos.org/dev/packages/ifpack2/doc/html/classIfpack2_1_1ILUT.html#aee2011b313e3070ee43b2cfc2d183634
            m_preconditioner->setParameters(precondOptions);
            m_preconditioner->initialize();

            // Specify the deposition problem
            m_problem = rcp(new problem_type(m_matrix, m_solution, m_rhs));
            if (!m_preconditioner.is_null())
            {
                m_problem->setRightPrec(m_preconditioner);
            }
            m_problem->setProblem();
            m_solver->setProblem(m_problem);
        }
        template <MeshObject M>
	NearestNeighborProblem::NearestNeighborProblem(M& domain, int nLayer)
        {
            set_comm_for_parallel();

            auto sizes = set_domain_sizes(domain,nLayer);

            m_map = init_local_global_index_map(domain,sizes);

            IndexTracker idx_tracker{domain,sizes};

            m_graph = init_m_graph(domain,idx_tracker,sizes);

            set_nonzero_elements<M>(domain, idx_tracker, sizes);

            fill_matrix_rhs_solution<M>();

            build_Belos_solver();

            set_problem<M>();
        }
    }
}
