# Credit — agent instructions

## Project structure

- `src/pyvol/` — Python-facing package (import as `pyvol`). Contains `market_data.py` and the compiled extension `cpy_credit.cpython-*-darwin.so` (build output, gitignored — regenerate with CMake).
- `src/credit/` — C++20 source for the pybind11 extension module `cpy_credit`. Not a Python package despite the `__init__.py`. Headers: `params.hpp` (parameter structs), `stvol.hpp` (Heston model), `risk_neutral.hpp` (Merton/FPT/Black-Scholes), `utils.hpp` (bracketing root solver via Boost TOMS 748), `thread_pool.hpp` (fixed-size pool used by `fitHeston`).
- `tests/` — `cpp_tests/` (Boost.Test suites: `HestonTest.cpp`, `RootSolverTest.cpp`), `py_tests/` (scripts, no framework; hit live market data).

## Build (C++ extension)

Requires NLopt and Boost via Homebrew:

```sh
cmake -B build && cmake --build build
```

Post-build copies the `.so` into `src/pyvol/`. Rebuild when changing C++ sources or after fresh clone. CMake picks whatever Python it finds first; to build against the project venv pass `-DPython3_EXECUTABLE=$PWD/.venv/bin/python`.

## Python environment

```sh
pdm install
```

Python pinned to `==3.12` in `pyproject.toml`. Virtualenv is at `.venv/`.

## Tests

- C++: `ctest --test-dir build --output-on-failure` (or run `./build/HestonTest`, `./build/RootSolverTest` directly)
- Python: `python tests/py_tests/heston_fit.py` (needs the built `.so` and live Yahoo Finance / US Treasury access)

No lint, no formatter, no typecheck config.

## CI

GitHub Actions at `.github/workflows/ci.yml`. Runs on `macos-latest`. Steps: install Homebrew deps (nlopt, boost, pybind11) → Python 3.12 → PDM → build C++ extension → CTest → Python script (with `continue-on-error: true` — it depends on live Yahoo Finance / US Treasury APIs and may be flaky).

## Conventions

- C++20 (`CMakeLists.txt` sets the standard; `std::numbers`, concepts in use).
- Prefer value semantics for parameter structs (`HestonUnderlying` is held by value in `HestonCallMdl`; no `unique_ptr` at the pybind boundary).
- Root finding goes through `bracket_root` in `utils.hpp` — a bracketing solver chosen because Black-Scholes objectives have flat regions that break secant-style iterations; keep new solvers bracketing unless there is a measured reason not to.
