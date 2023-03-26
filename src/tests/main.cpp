#include "gtest/gtest.h"

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
    result = RUN_ALL_TESTS();

//    MPI_Finalize();

    return result;
}
