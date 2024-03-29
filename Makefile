# Variables
PYREQS := $(shell python3 -m pybind11 --includes)

# Compiler flags dependent on underlying OS
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Linux)
	PYLINK_FLAGS := -fPIC
else ifeq ($(UNAME_S), Darwin) # macOS
	PYLINK_FLAGS := -undefined dynamic_lookup
else
	$(error Unspported OS: $(UNAME_S))
endif

# Commands to conduct builds of various tests and output libraries as needed.
pymodule:
	g++ -O3 -Wall -shared -std=c++20 $(PYLINK_FLAGS) $(PYREQS) src/credit/cpy_credit.cpp src/credit/risk_neutral.cpp -o out/cpy_credit.so

pytests:
	g++ -Wall -shared -std=c++20 $(PYLINK_FLAGS) $(PYREQS) src/credit/cpy_credit.cpp src/credit/risk_neutral.cpp -o out/cpy_credit.so -g

pyinit:
	mkdir -p .venv
	pipenv install

ctests: 
	g++ -Wall -std=c++20 test/option_price_unit.cpp src/credit/risk_neutral.cpp -o out/option_price_unit.out -g
	g++ -Wall -std=c++20 test/wang_transform_unit.cpp src/credit/risk_neutral.cpp -o out/wang_transform_unit.out -g
	g++ -Wall -std=c++20 test/ROL_unit.cpp src/credit/risk_neutral.cpp -o out/ROL_unit.out -g
	g++ -Wall -std=c++20 test/implied_asset_unit.cpp src/credit/risk_neutral.cpp -o out/implied_asset_unit.out -g
	g++ -Wall -std=c++20 test/secant_test.cpp -o out/secant_test.out -g

