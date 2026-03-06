# Repository Guidelines

## Project Structure & Module Organization
`src/` contains the C++ modules (`Generator`, `MeshingTree`, `VoronoiMesh`, `Methods`) that sample, mesh, and export. `include/` hosts their headers so other tools can reuse the API. `data/` holds OBJ/CSV fixtures, and the top-level OBJ files (e.g., `surface_mesh.obj`, `voronoi.obj`) document reference cases. Experimental helpers live in folders such as `block/` and `face_method/`. `build/` is the out-of-source CMake output—regenerate it via the configure command below rather than checking in its contents.

## Build, Test, and Development Commands
- `cmake -S . -B build` configures CMake, including the vcpkg toolchain, and writes build files into `build/`.
- `cmake --build build` compiles `vorocrust` (honoring the `USE_CGAL` flag).
- `./build/vorocrust` runs the binary; feed it meshes or seeds from `data/` or `block/`.
- `cmake --build build --target clean` removes artifacts so you can reconfigure from scratch.

## Coding Style & Naming Conventions
Follow the existing C++17 surface: 4-space indentation, braces at the opening line of functions, and headers grouped with project includes before third-party. Class names remain `PascalCase` (e.g., `MeshingTree`), and helpers match the casing already present in the file. Use descriptive variable names for geometry data (`Eigen::MatrixXd V`) and keep helper methods short to limit nesting. There is no enforced formatter yet; if you run `clang-format`, keep the diff readable and consistent with surrounding files.

## Testing Guidelines
No automated tests exist yet. When adding coverage, place executables under `tests/`, register them with `add_test`, and name the source `test_<feature>.cpp` or `FeatureTest.cpp`. Store any required fixtures in `data/` or `tests/data/`. Until tests are available, validate changes by running `./build/vorocrust` on representative OBJ/CSV inputs and inspecting the generated meshes or seed CSV.

## Commit & Pull Request Guidelines
Commits stay short and present tense (history shows single-topic summaries like `color error` or `holes`). Reference related issue numbers in the body if available. Pull requests should describe what changed, why, and how to reproduce the results (data path, command, expected output). Document mesh differences in words or with reference OBJ files; add screenshots only when a visual comparison clarifies the impact.
