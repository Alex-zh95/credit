/* @Filename:       cpy_credit.cpp
 * @Desrciption:    Implement exports for various functions within risk_neutral.h for Python interface.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // Needed for binding std::vector
namespace py = pybind11;

#include "risk_neutral.h"

PYBIND11_MODULE(cpy_credit, m) {
    m.doc() = "Module containing procedures for structural credit models and capital determination.";

    m.def(
        "get_vanilla_asset_volatility",
        &get_vanilla_asset_volatility,
        "Numerically derive asset volatility given debt and equity",
        py::arg("E"), py::arg("sigma_e"), py::arg("L"), py::arg("r"), py::arg("t") = 1, py::arg("n_iter") = 50
    );

    m.def(
        "get_vanilla_default_probability",
        &get_vanilla_default_probability,
        "Attain risk-neutral default probability",
        py::arg("a0"), py::arg("mu_a"), py::arg("sigma_a"), py::arg("L"), py::arg("t") = 1
    );

    m.def(
        "wang_transform",
        &wang_transform,
        "Apply Wang transform to provided probability vector",
        py::arg("P"), py::arg("sharpe_ratio"), py::arg("inverse") = false
    );

    m.def(
        "get_returns_with_put",
        &get_returns_with_put,
        "Calculate expected returns with protection from put option",
        py::arg("y"), py::arg("y_var"), py::arg("put"), py::arg("r")
    );

    m.def(
        "get_min_capital_ROL",
        &get_min_capital_ROL,
        "Obtain minimum rates on line. Set p=0 and i=risk-free rate if not using options",
        py::arg("y"), py::arg("p"), py::arg("i")
    );
}
