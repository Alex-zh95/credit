# Credit — agent instructions

## Project structure

- `src/pyvol/` — Python-facing package (import as `pyvol`). Contains `cpy_credit.cpython-312-darwin.so` (pre-built), `garch11.py`, `market_data.py`.
- `src/credit/` — C++17 source for pybind11 extension module `cpy_credit`. Not a Python package despite the `__init__.py`.
- `tests/` — `cpp_tests/` (raw C++ executables), `py_tests/` (scripts, no framework).

## Build (C++ extension)

Requires NLopt and Boost.Math via Homebrew (Apple Silicon paths hardcoded in `CMakeLists.txt`):

```sh
cmake -B build && cmake --build build
```

Post-build copies the `.so` into `src/pyvol/`. Rebuild when changing C++ sources or after fresh clone.

## Python environment

```sh
pdm install
```

Python pinned to `==3.12` in `pyproject.toml`. Virtualenv is at `.venv/`.

## Tests

- C++: `cd build && ctest` (or `./build/HestonTest`, `./build/SecantTest`)
- Python: `python tests/py_tests/garch_test.py` or `python tests/py_tests/heston_fit.py`

No lint, no formatter, no typecheck config.

## CI

GitHub Actions at `.github/workflows/ci.yml`. Runs on `macos-latest` (Apple Silicon) because `CMakeLists.txt` hardcodes `/opt/homebrew/` paths. Steps: install Homebrew deps (nlopt, boost) → Python 3.12 → PDM → build C++ extension → CTest. Python tests run with `continue-on-error: true` — they depend on live Yahoo Finance / US Treasury APIs and may be flaky.

## Build artifacts

`build/` is gitignored. The compiled `.so` at `src/pyvol/cpy_credit.cpython-312-darwin.so` is tracked — regenerate it after C++ changes.
