# Commands to conduct builds of various tests and output libraries as needed.

tests: 
	g++ -Wall -std=c++20 test/option_price_unit.cpp source/risk_neutral.cpp -o out/option_price_unit.out
	g++ -Wall -std=c++20 test/wang_transform_unit.cpp source/risk_neutral.cpp -o out/wang_transform_unit.out
