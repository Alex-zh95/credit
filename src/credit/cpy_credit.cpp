/* @Filename:       cpy_credit.cpp
 * @Description:    Implement exports for various functions within risk_neutral.h for Python interface.
 */

#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>           // Needed for binding STL containers
namespace py = pybind11;

#include "risk_neutral.hpp"
#include "stvol.hpp"

/* C++ Interfaces for StVol::HestonCallMdl 
 * This is due to not being able to pass std::unique_ptr<...> as arguments as Python does not implement such concepts.
 */

// Simplified accessor to get a default probability via training of Heston model on market data
StVol::Underlying fit_Hst(
    double S0,
    std::vector<double> Ks,
    std::vector<double> Ts,
    std::vector<double> Ps
)
{
    std::vector<StVol::HestonCallMdl> mdls;

    std::vector<double> initial_guess = {
        0.1,    // v0
        3.0,    // alpha
        0.05,   // vTheta
        0.03,   // vSig
        0.03,   // vLambda
        -0.8,   // vRho
        0.03,   // rf
    };

    auto nOptions = Ks.size();

    // Set up the vector of models for fit
    for (decltype(nOptions) j = 0; j < nOptions; ++j)
    {
        auto _U = std::make_unique<StVol::Underlying>();
        _U->S0 = S0;
        _U->v0 = initial_guess[0];
        _U->alpha = initial_guess[1];
        _U->vTheta = initial_guess[2];
        _U->vSig = initial_guess[3];
        _U->vLambda = initial_guess[4];
        _U->rho = initial_guess[5];
        _U->rf = initial_guess[6];

        StVol::HestonCallMdl mdl(std::move(_U), Ks[j]);
        mdl.calc_option_price();

        mdls.push_back(std::move(mdl));
    }

    // Fitted parameters to be stored in another StVol::Underlying
    auto V = StVol::market_calibration(mdls, Ps, initial_guess);

    return *V;
}

double get_Hst_call_price(
    StVol::Underlying U0,
    double strike,
    double t = 1
)
{
    auto _U = std::make_unique<StVol::Underlying>();
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

double get_Hst_default_probability(
    StVol::Underlying U_fitted,
    double strike
)
{
    auto _U = std::make_unique<StVol::Underlying>();
    _U->S0 = U_fitted.S0;
    _U->v0 = U_fitted.v0;
    _U->alpha = U_fitted.alpha;
    _U->vTheta = U_fitted.vTheta;
    _U->vSig = U_fitted.vSig;
    _U->vLambda = U_fitted.vLambda;
    _U->rho = U_fitted.rho;
    _U->rf = U_fitted.rf;

    StVol::HestonCallMdl mdl(std::move(_U), strike);
    mdl.calc_option_price();

    return (1. - mdl.get_rn_exercise_probability());
}

/* Exposing definitions to Python. */
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

    m.def(
        "fit_Hst",
        &fit_Hst,
        "Fit a Heston model to available call options data",
        py::arg("S0"), py::arg("Ks"), py::arg("Ts"), py::arg("Ps")
    );

    m.def(
        "get_Hst_default_probability",
        &get_Hst_default_probability,
        "Attain risk-neutral default probability (Heston).",
        py::arg("U_fitted"), py::arg("strike")
    );

    m.def(
        "get_Hst_call_price",
        &get_Hst_call_price,
        "Calculate price of a European vanilla call option under Heston model.",
        py::arg("U0"), py::arg("strike"), py::arg("t") = 1.0
    );
}
