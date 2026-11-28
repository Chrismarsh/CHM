#include "gtest/gtest.h"
#include <Kokkos_Core.hpp>

#ifdef USE_MPI
    #include <boost/mpi.hpp>
#endif

int main(int argc, char* argv[])
{
    int result = 0;


    ::testing::InitGoogleTest(&argc, argv);

#ifdef USE_MPI
    boost::mpi::environment _mpi_env;
    boost::mpi::communicator _comm_world;

    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
    if (_comm_world.rank() != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }
#endif
//    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
    result = RUN_ALL_TESTS();
    Kokkos::finalize();

//    MPI_Finalize();

    return result;
}
