#ifndef TEMPLATE_UTILS_H
#define TEMPLATE_UTILS_H
/* @Filename:        template_utils.hpp
 * @Description:     Generic templated functions/solvers.
 */

/* @Description:    Root solver via the secant method.
 *
 * @Params:         Func& f:            Equation to solve for roots
 *                  T x0:               Initial guess for the function root
 *                  T tol = 1e-8:       Tolerance level - stop if reached
 *                  int n_iter = 50:    Maximum number of iterations before finishing
 *
 * @Returns:        T x:                Root solution
 */
template <typename Func, typename T>
T secant_root(Func &&f, T x0, T tol = 1e-8, int n_iter = 50)
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
