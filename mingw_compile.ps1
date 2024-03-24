# Compile script for cpy_credit under Windows.
# Note: Change the versioning names as needed for your setup

# Binding dependencies
$home_dir = [System.Environment]::GetEnvironmentVariable("USERPROFILE")
$compiler_path = "$home_dir\llvm-mingw-20240320-msvcrt-x86_64\bin\g++.exe"

$python_includes = "$home_dir\AppData\Local\Programs\Python\Python312\Include"
$pybind11_includes = "$home_dir\AppData\Local\Programs\Python\Python312\Lib\site-packages\pybind11\include"
$python_lib_dll = "$home_dir\AppData\Local\Programs\Python\Python312"
$boost_lib_path = "$home_dir\boost_1_82_0"

# Run the following command to compile
Invoke-Expression "$compiler_path -O3 -Wall -shared -std=c++20 -fPIC -I$python_includes -I$pybind11_includes -I$boost_lib_path .\source\cpy_credit.cpp .\source\risk_neutral.cpp -o .\out\cpy_credit.pyd -L$python_lib_dll -lpython312 -static-libgcc -static-libstdc++"