# Commands to conduct builds of various tests and output libraries as needed.

tests: 
	g++ -Wall -std=c++20 test/option_price_unit.cpp source/risk_neutral.cpp -o out/option_price_unit.out
	g++ -Wall -std=c++20 test/wang_transform_unit.cpp source/risk_neutral.cpp -o out/wang_transform_unit.out
	g++ -Wall -std=c++20 test/ROL_unit.cpp source/risk_neutral.cpp -o out/ROL_unit.out

pyinit:
	mkdir -p .venv
	pipenv install

pytests:
	g++ -Wall -shared -std=c++20 -fPIC $(python -m pybind11 --includes) source/cpy_credit.cpp source/risk_neutral.cpp -o cpy_credit.so
