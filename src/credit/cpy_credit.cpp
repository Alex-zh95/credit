/* @Filename:       cpy_credit.cpp
 * @Description:    Implement exports for various functions within risk_neutral.h for Python interface.
 */

#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>           // Needed for binding STL containers
namespace py = pybind11;

#include "risk_neutral.hpp"
#include "stvol.hpp"

// Helper function for creating unique pointers to StVol's Underlying struct
std::unique_ptr<StVol::HestonCallMdl> HestonConstructorHelper(std::unique_ptr<StVol::Underlying> _uPtr, double _strike, double _maturity = 1.0)
{
    return std::make_unique<StVol::HestonCallMdl>(std::move(_uPtr), _strike, _maturity);
}

PYBIND11_MODULE(cpy_credit, m) {
    m.doc() = "Module containing procedures for structural credit models and capital determination.";

    // Bindings to risk_neutral.hpp
    m.def(
        "get_asset_volatility",
        &get_asset_volatility,
        "Numerically derive asset volatility given debt and equity",
        py::arg("E"), py::arg("sigma_e"), py::arg("L"), py::arg("r"), py::arg("t") = 1, py::arg("n_iter") = 50
    );

    m.def(
        "get_vanilla_default_probability",
        &get_vanilla_default_probability,
        "Attain risk-neutral default probability (Merton)",
        py::arg("a0"), py::arg("rf"), py::arg("sigma_a"), py::arg("L"), py::arg("t") = 1
    );

    m.def(
        "get_fpt_default_probability",
        &get_fpt_default_probability,
        "Attain risk-neutral default probability (First Passage Time).",
        py::arg("a0"), py::arg("rf"), py::arg("sigma_a"), py::arg("L"), py::arg("delta") = 0., py::arg("gamma") = 0., py::arg("t") = 1.
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

    // Bindings to stvol.hpp
    py::class_<StVol::Underlying>(m, "Underlying")
        .def(py::init<>())
        .def_readwrite("S0", &StVol::Underlying::S0)
        .def_readwrite("v0", &StVol::Underlying::v0)
        .def_readwrite("alpha", &StVol::Underlying::alpha)
        .def_readwrite("vSig", &StVol::Underlying::vSig)
        .def_readwrite("rho", &StVol::Underlying::rho)
        .def_readwrite("vTheta", &StVol::Underlying::vTheta)
        .def_readwrite("vLambda", &StVol::Underlying::vLambda)
        .def_readwrite("rf", &StVol::Underlying::rf);

    // Note: We will have to redefine how to expose this class or its calculations via helper C++ functions - i.e. probably separate handling of StVol::Underlying objects and wrappers for HestonMdlCall to avoid exposing unique pointers
    // py::class_<StVol::HestonCallMdl>(m, "HestonCallMdl")
    //     // Will only register the custom constructor - other constructors will be registered automatically
    //     .def(
    //         // ERROR: Illegal to use uniqe ptrs as function arguments for pybind
    //         py::init<std::unique_ptr<StVol::Underlying>, double, double>(),
    //         py::arg("_underlying"), py::arg("_K"), py::arg("_t") = 1.0
    //     )

    //     // Getters and setters
    //     .def_property("t", &StVol::HestonCallMdl::get_maturity, &StVol::HestonCallMdl::set_maturity)
    //     .def_property("P", &StVol::HestonCallMdl::get_option_price, nullptr)
    //     .def_property("K", &StVol::HestonCallMdl::get_strike, &StVol::HestonCallMdl::set_strike)

    //     // Public member functions
    //     .def("get_rn_exercise_probability", &StVol::HestonCallMdl::get_rn_exercise_probability)
    //     .def("get_delta", &StVol::HestonCallMdl::get_delta)
    //     .def("calc_option_price", &StVol::HestonCallMdl::calc_option_price);
}
