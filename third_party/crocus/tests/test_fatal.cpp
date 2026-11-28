#include <cstdlib>

extern "C" void raise_crocus_error();

int main()
{
    // CTest expects a failure. Reaching an exit handler would leave a CHM
    // rank inside collective cleanup, so make that path fail the test.
    if (std::atexit([] { std::_Exit(EXIT_SUCCESS); }) == 0)
        raise_crocus_error();
    return EXIT_SUCCESS;
}
