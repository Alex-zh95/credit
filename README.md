# Credit

Repository of code primarily relating to modeling credit risk as a for of insurance. The ideas presented here can also be extended to apply to reinsurance in an analogous vein: the P&L of an underwriting firm can be used to evaluate likelihood of a stop loss reinsurance (or approximately a very high excess layer) policy responding.

## Applications

1. Direct modeling of individual credit risks via an internal structural models.
2. Generating minimum premiums for high layer reinsurance from a capital point of view.

## Premium generation

### Default risk via structural model (loss cost)

- Use Merton model to evaluate probabilities of default.
- Under this framework, assume that assets, $A$, compriese purely of equities, $E$, and (short-term) liabilities, $L$.
    - After period of consideration (e.g. one year), value to shareholders is $\max(A-L, 0)$, which can be viewed as holding a call option on assets with strike price as the face value of liabilities.
    - Call option can be evaluated via solving Black-Scholes equation.
- Challenge is asset volatility, which is not observable.
    1. Start with some initial guess of asset volatility.
    2. Using this asset volatility, imply the asset price by solving the Merton model backwards. We perform this numerically.
    3. With a new vector of asset prices, update our guess of asset volatility. Repeat steps 2 and 3 until convergence.
- The probability of default is equivalent to the probability of the call option not being exercised (i.e. expiring out of the money and so equity is worthless).
    - This can be viewed as insurance rate on line post discounting to present value.
    - The probability derived is in risk-neutral world. We can apply risk-adjustment using the Wang transform with the appropriate market price of risk.

Idea for this method of iteratively solving for implied asset volatility, see [1](#references).

Default can also be defined as the point at which a reinsurance contract would pay out.

### Asymptotically minimum premiums (cost of capital)

- Consider premium to be a sum of risk load and discounted expected loss. The focus would be on cost of capital.
- One way to build up a cost of capital is to select an investment and buy a put option to protect against the case that the investment does not meet the risk-free rate.
    - The put option price contributes to the cost of capital, however, there are offsets in the form form increase in return and decrease in variance of selected investment.
    - Risk load is given by:

    $$R = S\frac{(1+y)(1+p) - (1+i)}{(1+y)(1+r)} - L \left( \frac{1}{1+y} - \frac{1}{1+r} \right)$$

    where $S$ is the safety level of loss (e.g. the limit of loss), $L$ the expected loss on portfolio, $y$ is the return on equity target, $p$ the unit price of put option (i.e. put price per dollar of selected investment), $i$ the mean investment return with option protection (simulate) and $r$ is the risk free rate.

    - Rate on line is obtained by dividing both sides of the above by $S$.
    - Taking $L \rightarrow 0$ yields a minimum rate on line (applicable if we are looking at very high excesses).
- Incorporating the put option will allow for investment market sentiment in affecting risk load more dynamically.
    - We can start off by using underwriter-led rates on line, looking at available options for the types of investments desired evaluate and adjust back the risk load.
    - 

## Volatility

Volatility is hard to observe generally. Classic Black-Scholes implies constant volatility over the period of observation, which may not hold. To allow for volatility to be a stochastic quantity, use the Heston model or time series. A common time series implementation is the GARCH(1,1) model.

## References

1. [Idea for implied asset volatility](https://www.bradfordlynch.com/blog/2017/05/20/ProbabilityOfDefault.html)
2. Hull, J. (2021): ‘Options, Futures, and Other Derivatives, Global Edition’.
3. Cummins, J.D. (1990): ‘Asset Pricing Models and Insurance Ratemaking’, ASTIN Bulletin, 20(2), pp. 125–166. doi:10.2143/AST.20.2.2005438
4. [Basket call options](https://quant.stackexchange.com/questions/57235/do-basket-options-have-a-closed-form-valuation-formula)
5. Mikhailov S., Nögel U. (2003): ‘Heston’s Stochastic Volatility Model Implementation, Calibration and Some Extensions’.
