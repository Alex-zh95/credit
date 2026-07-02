/* @Filename:       cpy_credit.cpp
 * @Description:    Implement exports for various functions within risk_neutral.h for Python
 * interface.
 */

#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for binding STL containers
#include <tuple>
namespace py = pybind11;

#include "params.hpp"
#include "risk_neutral.hpp"
#include "stvol.hpp"

/* C++ Interfaces for StVol::HestonCallMdl
 * This is due to not being able to pass std::unique_ptr<...> as arguments
 * as Python does not implement such concepts.
 */

double get_Heston_call_price(StVol::HestonUnderlying U0, double strike, double t = 1) {
    auto _U = std::make_unique<StVol::HestonUnderlying>();
    _U->S0 = U0.S0;
    _U->v0 = U0.v0;
    _U->alpha = U0.alpha;
    _U->vTheta = U0.vTheta;
    _U->vSig = U0.vSig;
    _U->vLambda = U0.vLambda;
    _U->rho = U0.rho;
    _U->rf = U0.rf;

    StVol::HestonCallMdl mdl(std::move(_U), strike);
    mdl.calc_option_price();

    return mdl.get_option_price();
}

std::tuple<double, StVol::HestonUnderlying>
get_Heston_default_probability(StVol::HestonUnderlying U_fitted, double asset, double debt,
                               double maturity = 1.0) {
    // Build a Heston structural model representing the equity characteristics
    auto _U = std::make_unique<StVol::HestonUnderlying>();
    _U->S0 = U_fitted.S0;
    _U->v0 = U_fitted.v0;
    _U->alpha = U_fitted.alpha;
    _U->vTheta = U_fitted.vTheta;
    _U->vSig = U_fitted.vSig;
    _U->vLambda = U_fitted.vLambda;
    _U->rho = U_fitted.rho;
    _U->rf = U_fitted.rf;

    StVol::HestonCallMdl call(std::move(_U), debt);

    // Convert structural model for equity into structural model for asset
    auto _V = StVol::HestonAssetVolatilityImplied(call, asset, debt, maturity);

    StVol::HestonCallMdl structure(std::move(_V), debt);

    return std::make_tuple(1. - structure.get_rn_exercise_probability(),
                           structure.get_underlying());
}

/* Exposing definitions to Python. */
PYBIND11_MODULE(cpy_credit, m) {
    m.doc() =
        "Module containing procedures for structural credit models and capital determination.";

    // Bindings to risk_neutral.hpp
    m.def("get_asset_volatility", &get_asset_volatility,
          "Numerically derive asset volatility given debt and equity", py::arg("E"),
          py::arg("sigma_e"), py::arg("L"), py::arg("r"), py::arg("t") = 1, py::arg("n_iter") = 50);

    m.def(
        "get_vanilla_default_probability",
        [](double a0, double rf, double sigma_a, double L, double t) {
            return get_vanilla_default_probability(AssetDefaultParams(a0, rf, sigma_a, L, t));
        },
        "Attain risk-neutral default probability (Merton)", py::arg("a0"), py::arg("rf"),
        py::arg("sigma_a"), py::arg("L"), py::arg("t") = 1);

    m.def(
        "get_fpt_default_probability",
        [](double a0, double rf, double sigma_a, double L, double q, double gamma, double t) {
            return get_fpt_default_probability(AssetDefaultParams(a0, rf, sigma_a, L, t, q), gamma);
        },
        "Attain risk-neutral default probability (First Passage Time).", py::arg("a0"),
        py::arg("rf"), py::arg("sigma_a"), py::arg("L"), py::arg("delta") = 0.,
        py::arg("gamma") = 0., py::arg("t") = 1.);

    m.def("wang_transform", &wang_transform, "Apply Wang transform to provided probability vector",
          py::arg("P"), py::arg("sharpe_ratio"), py::arg("inverse") = false);

    m.def("get_returns_with_put", &get_returns_with_put,
          "Calculate expected returns with protection from put option", py::arg("y"),
          py::arg("y_var"), py::arg("put"), py::arg("r"));

    m.def("get_min_capital_ROL", &get_min_capital_ROL,
          "Obtain minimum rates on line. Set p=0 and i=risk-free rate if not using options",
          py::arg("y"), py::arg("p"), py::arg("i"));

    // Bindings to stvol.hpp
    py::class_<StVol::HestonUnderlying>(m, "Underlying")
        .def(py::init<>())
        .def_readwrite("S0", &StVol::HestonUnderlying::S0)
        .def_readwrite("v0", &StVol::HestonUnderlying::v0)
        .def_readwrite("alpha", &StVol::HestonUnderlying::alpha)
        .def_readwrite("vSig", &StVol::HestonUnderlying::vSig)
        .def_readwrite("rho", &StVol::HestonUnderlying::rho)
        .def_readwrite("vTheta", &StVol::HestonUnderlying::vTheta)
        .def_readwrite("vLambda", &StVol::HestonUnderlying::vLambda)
        .def_readwrite("rf", &StVol::HestonUnderlying::rf);

    m.def("fit_Heston", &StVol::fitHeston, "Fit a Heston model to available call options data",
          py::arg("S0"), py::arg("Ks"), py::arg("rfs"), py::arg("Ts"), py::arg("Ps"),
          py::arg("Volumes"));

    m.def("get_Heston_default_probability", &get_Heston_default_probability,
          "Attain risk-neutral default probability (Heston). Also returns a copy of derived asset "
          "vol model.",
          py::arg("U_fitted"), py::arg("asset"), py::arg("debt"), py::arg("maturity") = 1.0);

    m.def("get_Heston_call_price", &get_Heston_call_price,
          "Calculate price of a European vanilla call option under Heston model.", py::arg("U0"),
          py::arg("strike"), py::arg("t") = 1.0);
}
