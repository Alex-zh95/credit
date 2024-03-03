# Credit

Repository of code relating to modeling credit risk.

## Applications

1. Direct modeling of individual credit risks via an internal structural models.
2. Generating minimum premiums for high layer reinsurance from a capital point of view.

## Credit risk via structural model

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

Idea for this method of iteratively solving for implied asset volatility [here](https://www.bradfordlynch.com/blog/2017/05/20/ProbabilityOfDefault.html)

## Asymptotically minimum premiums

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
