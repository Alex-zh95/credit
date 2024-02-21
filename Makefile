# Variables
PYREQS := $(shell python3 -m pybind11 --includes)

# Commands to conduct builds of various tests and output libraries as needed.
pymodule:
	g++ -O3 -Wall -shared -std=c++20 -fPIC $(PYREQS) source/cpy_credit.cpp source/risk_neutral.cpp -o out/cpy_credit.so

pytests:
	g++ -Wall -shared -std=c++20 -fPIC $(PYREQS) source/cpy_credit.cpp source/risk_neutral.cpp -o out/cpy_credit.so

pyinit:
	mkdir -p .venv
	pipenv install

ctests: 
	g++ -Wall -std=c++20 test/option_price_unit.cpp source/risk_neutral.cpp -o out/option_price_unit.out
	g++ -Wall -std=c++20 test/wang_transform_unit.cpp source/risk_neutral.cpp -o out/wang_transform_unit.out
	g++ -Wall -std=c++20 test/ROL_unit.cpp source/risk_neutral.cpp -o out/ROL_unit.out
	g++ -Wall -std=c++20 test/implied_asset_unit.cpp source/risk_neutral.cpp -o out/implied_asset_unit.out

