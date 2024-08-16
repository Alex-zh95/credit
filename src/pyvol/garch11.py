"""
@Filename:      garch11.py
@Description:   Script used for fitting a GARCH(1,1) model to imply volatility. This method is more useful if there are no implied volatility data available.
"""

import numpy as np
import arch
from scipy import stats

from typing import Tuple


class GarchVol:
    '''
    Encapsulate instrument ticker information and accompanying modeling statistics pertaining to volatility.

    Note attributes are in % form (i.e. sigma=3.5 means 3.5%)

    @attributes:    sigma:  float:              empirical or garch-estimated volatility of return
                    mu:     float:              empirical or garch-estimated mean return
    '''

    def __init__(self, X: np.ndarray):
        '''
        Initialize the model.

        @params:        x:      np.ndarray:         vector of ticker information to store
        '''
        self._X = X
        self._dlogX = np.zeros(X.shape)
        self._mu = 0.0
        self._sigma = 0.0

        # Store volatility model in this param
        self._sigmdl = None

        # Residual variances - we run diagnostics here
        self._epsilons = np.zeros(X.shape[0] - 1)

    def _reforecast(self, horizon: int = 1):  # -> arch.ARCHModelForecast:
        '''
        Private fn: Perform forecast using _sigmdl but with custom horizon.

        Returns a forecast object for further use.
        '''
        if self._sigmdl is None:
            raise ValueError('_sigmdl not yet defined. Run method fit_garch at least once to initiate.')
        return self._sigmdl.forecast(horizon=horizon)

    def _log_rts(self) -> None:
        '''
        Private fn: Calculate the log returns for the underlying values and store the mean and volatilities.
        '''

        # arch documentation recommends scaling to int percentages
        self._dlogX = np.diff(np.log(self._X)) * 100
        self._mu = np.mean(self._dlogX)
        self._sigma = np.std(self._dlogX)

    def fit_garch(self, start: int = 0, end: int = 0, p: int = 1, q: int = 1, h: int = 252) -> None:
        '''
        Fit the GARCH(p,q) model to the underlying return data.

        @Params:    start:  int = 0:      Clip the start of the vector (no start clip by default)
                    end:    int = 0:      Clip the end of the vector (no end clip by default)
                    p:      int = 1:      ARCH p parameter (i.e. the AR(p) for underlying variance)
                    q:      int = 1:      GARCH q parameter (i.e. the MA(q) for underlying variance)
                    h:      int = 252:    Horizon for projection (trading year default)

        @Returns:   None
        '''
        dlogX = self._dlogX[start:len(self._dlogX)] if end == 0 else self._dlogX[start:end]
        am = arch.arch_model(dlogX, vol='GARCH', p=1, q=1, dist="normal")

        # Store handle to model
        self._sigmdl = am.fit()

        # Store results (for whole trading year)
        forecasts = self._reforecast(horizon=h)
        self._sigma = np.sqrt(forecasts.variance.values[-1, :])
        self._mu = forecasts.mean.values[-1, :]

    def diagnostics(self) -> Tuple[float]:
        '''
        Null hypothesis under KS test: Residuals are normal-distributed.

        @Params:     None
        @Returns:    (float, float):         KS-statistic and p-value
        '''
        forecasts = self._reforecast(horizon=len(self._dlogX))
        self._epsilons = forecasts.mean.values[-1, :] - self._dlogX

        res_mdl = stats.norm.fit(self._epsilons)
        test_result = stats.kstest(self._epsilons, "norm", args=res_mdl)

        return (test_result.statistic, test_result.pvalue)

    @property
    def X(self) -> np.ndarray:
        '''
        Underlying values data.
        '''
        return self._X

    @property
    def mu(self) -> float:
        '''
        Return the mean of the underlying returns. If this has not been calculated at elast once, then return the empirical mean.
        '''
        if np.max(self._mu) == 0.0:
            self._log_rts()

        return self._mu

    @property
    def sigma(self) -> float:
        '''
        Return the volatility of underlying returns. If this has not been calculated at lease once, then use empirical volatility.
        '''
        if np.max(self._sigma) == 0.0:
            self._log_rts()

        return self._sigma
