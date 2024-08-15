"""
@Filename:      garch11.py
@Description:   Script used for fitting a GARCH(1,1) model to imply volatility.
"""

import numpy as np
from arch import arch_model


class GarchVol:
    def __init__(self, X: np.ndarray):
        '''
        Encapsulate instrument ticker information and accompanying modeling statistics pertaining to volatility.

        @Params:    X: np.ndarray:      Vector of ticker information to store
        @Returns:   None
        '''
        self._X = X
        self._dlogX = np.zeros(X.shape)
        self._mu = 0.0
        self._sigma = 0.0

        # Parameters for volatility model
        self._sigmdl = None
        self._sigvar = np.zeros(X.shape - 1)
        self._nu = 0.0

    def _log_rts(self) -> None:
        '''
        Private fn: Calculate the log returns for the underlying values and store the mean and volatilities.
        '''
        self._dX = np.diff(np.log(self._X))
        self._mu = np.mean(self._dX)
        self._sigma = np.std(self._dX)

    def fit_garch(self, start: int = 0, end: int = 0, p: int = 1, q: int = 1) -> None:
        '''
        Fit the GARCH(p,q) model to the underlying return data.

        @Params:    start:  int = 0:      Clip the start of the vector (no start clip by default)
                    end:    int = -1:     Clip the end of the vector (no end clip by default)
                    p:      int = 1:      ARCH p parameter (i.e. the AR(p) for underlying variance)
                    q:      int = 1:      GARCH q parameter (i.e. the MA(q) for underlying variance)

        @Returns:   None
        '''
        dlogX = self._dlogX[start:end]
        am = arch_model(dlogX, vol='GARCH', p=1, q=1, dist="normal")

        # Store handle to model
        self._sigmdl = am.fit(update_freq=5)  # Limit the stdout to every 5 iterations

        # Store results (1-step ahead)
        forecasts = self._sigmdl.forecast(horizon=1)
        self._sigma = forecasts.variance.values[-1, :]
        self._nu = forecasts.mean.values[-1, :]

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
        if self._mu == 0.0:
            self._log_rts()

        return self._mu

    @property
    def sigma(self) -> float:
        '''
        Return the volatility of underlying returns. If this has not been calculated at lease once, then use empirical volatility.
        '''
        if self._sigma == 0.0:
            self._log_rts()

        return self._sigma
