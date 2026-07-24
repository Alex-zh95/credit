/* @Filename:       cpy_credit.cpp
 * @Description:    Implement exports for various functions within risk_neutral.h for Python
 * interface.
 */

#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>  // Needed for binding std::tuple
#include <nanobind/stl/vector.h> // Needed for binding std::vector containers
#include <tuple>
#include <utility>
namespace nb = nanobind;
using namespace nb::literals;

#include "params.hpp"
#include "risk_neutral.hpp"
#include "stvol.hpp"

/* C++ Interfaces for StVol::HestonCallMdl - the model takes its underlying by
 * value, so nanobind's copies pass straight through.
 */

double get_Heston_call_price(StVol::HestonUnderlying U0, double strike, double t = 1) {
    StVol::HestonCallMdl mdl(std::move(U0), strike, t);
    mdl.calc_option_price();

    return mdl.get_option_price();
}

std::tuple<double, StVol::HestonUnderlying>
get_Heston_default_probability(StVol::HestonUnderlying U_fitted, double asset, double debt,
                               double maturity = 1.0) {
    // Build a Heston structural model representing the equity characteristics
    StVol::HestonCallMdl call(std::move(U_fitted), debt);

    // Convert structural model for equity into structural model for asset
    auto V = StVol::HestonAssetVolatilityImplied(call, asset, debt, maturity);

    StVol::HestonCallMdl structure(std::move(V), debt);

    return std::make_tuple(1. - structure.get_rn_exercise_probability(),
                           structure.get_underlying());
}

/* Exposing definitions to Python. */
NB_MODULE(cpy_credit, m) {
    m.doc() =
        "Module containing procedures for structural credit models and capital determination.";

    // Bindings to risk_neutral.hpp
    m.def("get_asset_volatility", &get_asset_volatility,
          "Numerically derive asset volatility given debt and equity", "E"_a, "sigma_e"_a, "L"_a,
          "r"_a, "t"_a = 1, "n_iter"_a = 50);

    m.def(
        "get_vanilla_default_probability",
        [](double a0, double rf, double sigma_a, double L, double t) {
            return get_vanilla_default_probability(AssetDefaultParams(a0, rf, sigma_a, L, t));
        },
        "Attain risk-neutral default probability (Merton)", "a0"_a, "rf"_a, "sigma_a"_a, "L"_a,
        "t"_a = 1);

    m.def(
        "get_fpt_default_probability",
        [](double a0, double rf, double sigma_a, double L, double q, double gamma, double t) {
            return get_fpt_default_probability(AssetDefaultParams(a0, rf, sigma_a, L, t, q), gamma);
        },
        "Attain risk-neutral default probability (First Passage Time).", "a0"_a, "rf"_a,
        "sigma_a"_a, "L"_a, "delta"_a = 0., "gamma"_a = 0., "t"_a = 1.);

    m.def("wang_transform", &wang_transform, "Apply Wang transform to provided probability vector",
          "P"_a, "sharpe_ratio"_a, "inverse"_a = false);

    m.def("get_returns_with_put", &get_returns_with_put,
          "Calculate expected returns with protection from put option", "y"_a, "y_var"_a, "put"_a,
          "r"_a);

    m.def("get_min_capital_ROL", &get_min_capital_ROL,
          "Obtain minimum rates on line. Set p=0 and i=risk-free rate if not using options", "y"_a,
          "p"_a, "i"_a);

    // Bindings to params.hpp: base class holding parameters common to all underlyings
    nb::class_<StandardUnderlying>(m, "StandardUnderlying")
        .def(nb::init<double, double, double, double, double, bool, double>(), "S0"_a, "K"_a, "r"_a,
             "sigma"_a, "t"_a = 1., "call"_a = true, "q"_a = 0.)
        .def(nb::init<>())
        .def_rw("S0", &StandardUnderlying::S0)
        .def_rw("K", &StandardUnderlying::K)
        .def_rw("r", &StandardUnderlying::r)
        .def_rw("sigma", &StandardUnderlying::sigma)
        .def_rw("t", &StandardUnderlying::t)
        .def_rw("call", &StandardUnderlying::call)
        .def_rw("q", &StandardUnderlying::q);

    // Bindings to stvol.hpp: HestonUnderlying extends StandardUnderlying
    nb::class_<StVol::HestonUnderlying, StandardUnderlying>(m, "Underlying")
        .def(nb::init<>())
        .def_rw("v0", &StVol::HestonUnderlying::v0)
        .def_rw("alpha", &StVol::HestonUnderlying::alpha)
        .def_rw("vSig", &StVol::HestonUnderlying::vSig)
        .def_rw("rho", &StVol::HestonUnderlying::rho)
        .def_rw("vTheta", &StVol::HestonUnderlying::vTheta)
        .def_rw("vLambda", &StVol::HestonUnderlying::vLambda)
        // Backwards-compatible alias: "rf" maps onto the inherited risk-free rate "r"
        .def_prop_rw(
            "rf", [](const StVol::HestonUnderlying& u) { return u.r; },
            [](StVol::HestonUnderlying& u, double rf) { u.r = rf; });

    m.def("fit_Heston", &StVol::fitHeston, "Fit a Heston model to available call options data",
          "S0"_a, "Ks"_a, "rfs"_a, "Ts"_a, "Ps"_a, "Volumes"_a);

    m.def("get_Heston_default_probability", &get_Heston_default_probability,
          "Attain risk-neutral default probability (Heston). Also returns a copy of derived asset "
          "vol model.",
          "U_fitted"_a, "asset"_a, "debt"_a, "maturity"_a = 1.0);

    m.def("get_Heston_call_price", &get_Heston_call_price,
          "Calculate price of a European vanilla call option under Heston model.", "U0"_a,
          "strike"_a, "t"_a = 1.0);
}
