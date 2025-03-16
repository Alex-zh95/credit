#ifndef UTILS_H
#define UTILS_H

 // Useful consts such as tolerances
 #ifndef TOL
 #define TOL 1e-3
 #endif

 #ifndef N_ITER
 #define N_ITER 50
 #endif

/* Description
 * -----------
 * Root solver via the secant ion root
 *  T tol = 1e-3:       Tolerance level - stop if reached
 *  int n_iter = 50:    Maximum number of iterations before finishing
 *
 * Returns:
 * --------        
 * T x:                Root solution
 */
template <typename Func, typename T>
T secant_root(Func &&f, T x0, T tol = TOL, int n_iter = N_ITER)
{
    // Define a neighborhood around the initial value x0
    auto x1 = x0 / 2;
    auto x2 = x0 * 2;

    for (auto iter = 0; iter < n_iter; ++iter)
    {
        x0 = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
        auto c = f(x0);

        if ((c < tol) && (c > -tol))
            break;

        if (c < -tol)
        {
            x1 = x2;
            x2 = x0;
        }
        else
        {
            x2 = x1;
            x1 = x0;
        }
    }

    return x0;
}

#endif
