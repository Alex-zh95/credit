#ifndef PARAMS_HPP
#define PARAMS_HPP

struct OptionParams {
    double S0;
    double K;
    double r;
    double sigma;
    double t;
    bool call;
    double q;

    OptionParams(double S0_, double K_, double r_, double sigma_,
                 double t_ = 1., bool call_ = true, double q_ = 0.)
        : S0(S0_), K(K_), r(r_), sigma(sigma_), t(t_), call(call_), q(q_) {}
};

struct AssetDefaultParams : OptionParams {
    // a0 maps to S0, rf maps to r, sigma_a maps to sigma, L maps to K
    AssetDefaultParams(double a0_, double rf_, double sigma_a_, double L_,
                       double t_ = 1., double q_ = 0.)
        : OptionParams(a0_, L_, rf_, sigma_a_, t_, true, q_) {}
};

#endif
