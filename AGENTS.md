# Repository Guidelines

## Project Structure & Module Organization
- `src/` contains implementation code, organized by domain: `coarsening/`, `partitioning/`, `refinement/`, `datastructures/`, and `utils/`.
- `include/` mirrors `src/` with public headers. Keep declarations in `include/...` and definitions in `src/...`.
- `CMakeLists.txt` defines a single executable target: `ConsMLP_lw`.
- `build/` is an out-of-source CMake build directory and should be treated as generated artifacts.

## Build, Test, and Development Commands
- Configure (Release by default):
  ```bash
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
  ```
- Build:
  ```bash
  cmake --build build -j
  ```
- Run locally (example):
  ```bash
  ./build/ConsMLP_lw <input.hgr> -k 4 -mode direct -init rand
  ```
- Show CLI options:
  ```bash
  ./build/ConsMLP_lw
  ```
  (prints usage when required args are missing)

## Coding Style & Naming Conventions
- Language standard is **C++11** (`CMAKE_CXX_STANDARD 11`).
- Use 4-space indentation and braces on the same line for functions/conditionals.
- Types/classes use `PascalCase` (e.g., `HypergraphHierarchy`); functions use `camelCase` (e.g., `parseArguments`); private members use trailing underscore (e.g., `config_`).
- Prefer matching header/source paths, e.g., `include/partitioning/Partitioner.h` with `src/partitioning/Partitioner.cpp`.

## Testing Guidelines
- No automated test target is currently defined (`CTest`/`add_test` not present).
- For changes, run a smoke test with a representative `.hgr` input and verify output metrics/logs are sane.
- If you add tests, place them under a dedicated `tests/` directory and wire them into CMake/CTest in the same PR.

## Commit & Pull Request Guidelines
- Git history is not available in this workspace snapshot, so no repository-specific commit pattern can be inferred.
- Use clear, imperative commit messages, preferably Conventional Commits style (e.g., `feat: add recursive split guard`).
- PRs should include: purpose, key design choices, commands run for validation, and sample output diffs when behavior changes.
