# Repository Guidelines

## Project Structure & Module Organization
- `src/` contains the C++20 engine, process modules, filters, and GoogleTest fixtures in `src/tests/`.
- `docs/` runs the Sphinx toolchain; supporting figures and manuals live in `docs/images/`, `manual/`, and `wiki/`.
- `functional_tests/` and `test_data/` store meshes, forcing inputs, and regression baselines—keep generated outputs outside git.
- `tools/`, `build-test/`, `spack-env.sh`, `third_party/`, and `CMake/` provide build helpers and vendored modules; touch them only when upgrading dependencies.

## Build, Test, and Development Commands
Activate the spack toolchain before configuring.
- `spack env activate chm` (create it once from `spack.yaml` if needed) loads the compiler/MPI stack.
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON -DUSE_MPI=ON` configures an out-of-tree build.
- `cmake --build build -j$(nproc)` produces `build/bin/CHM`; add `--target install` only when packaging.
- `ctest --test-dir build --output-on-failure` runs the GoogleTest suite; heavier `functional_tests` scenarios stay manual.
- `READTHEDOCS=True make -C docs html` regenerates developer docs.
- `mpirun -np 4 build/bin/CHM path/to/CHM.config` launches a run; tune ranks to the mesh size.

## Coding Style & Naming Conventions
- Use the repo `.clang-format`: Allman braces, four spaces, 120-column limit, left-aligned pointers; run `clang-format -i` before pushing.
- Follow existing naming: lower_snake_case for files/functions, CamelCase for types, and keep headers mirrored by `.cpp` implementations.
- Gate MPI/OpenMP code with `USE_MPI`/`USE_OPENMP` and prefer configuration through `config_file` parameters rather than hard-coded constants.

## Testing Guidelines
- Add unit coverage in `src/tests/test_*.cpp` using `TEST(Module, Scenario)`; favor small deterministic fixtures.
- Configure with `-DBUILD_TESTS=ON` and call `ctest`; keep integration workflows (`functional_tests/mesh_versioning`) in dedicated scripts or CI jobs.
- Update `test_data/` when introducing new physics inputs so regressions stay reproducible offline.

## Commit & Pull Request Guidelines
- Match the history style: concise, imperative summaries (`fix old path in pvd creation`), grouped by logical change.
- Reference issues (`Fixes #123`) and note any new runtime flags, configs, or mesh requirements.
- Pull requests should outline physics impact, confirm `spack env activate` + `cmake --build` + `ctest` results, and attach logs or screenshots when behavior shifts.
