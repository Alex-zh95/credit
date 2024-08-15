"""
@Filename:      garch11.py
@Description:   Script used for fitting a GARCH(1,1) model to imply volatility.

NOTE: Add a residuals checker as a diagnostic for GoF (e.g. resiudals should follow a N(0,1) distribution - use various tests of normality to verify this).
"""

import numpy as np
import arch


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

    def _reforecast(self, horizon: int = 1) -> arch.ARCHModelForecast:
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
        am = arch.arch_model(dlogX, vol='GARCH', p=1, q=1, dist="normal")

        # Store handle to model
        self._sigmdl = am.fit(update_freq=5)  # Limit the stdout to every 5 iterations

        # Store results (1-step ahead)
        forecasts = self._reforecast(horizon=1)
        self._nu = forecasts.variance.values[-1, :]
        self._vol_mean = forecasts.mean.values[-1, :]

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

    @property
    def nu(self) -> float:
        '''
        Return the volatility of the volatility (i.e. volatility of sigma).
        '''
        if self.nu == 0.0:
            raise Warning('Nu = 0. Possibly undefined.')

        return self.nu
