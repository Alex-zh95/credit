'''
Model layer: market data acquisition, Heston fitting and surface construction.

No tkinter/matplotlib imports live here. The heavy functions at module level
run inside a worker *process* (the nanobind calls hold the GIL, so a thread
would freeze the GUI); `VolSurfaceModel` owns the pool and returns futures.
'''

from __future__ import annotations

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import griddata
from scipy.special import ndtr

# Grid resolutions. The market surface is pure numpy (cheap, so fine);
# the model surfaces cost one C++ Heston integration per node (~0.7 ms each).
MARKET_GRID = (60, 40)          # (strikes, maturities)
MODEL_GRID = (36, 18)
ASSET_MONEYNESS = (0.50, 1.50)  # strike / asset-value range for the asset surface
MIN_QUOTES = 20


# --------------------------------------------------------------------------
# Data structures shared between model, controller and view
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class HestonParams:
    '''Plain (picklable) mirror of cpy_credit.Underlying.'''
    S0: float
    r: float
    q: float
    v0: float
    alpha: float
    vTheta: float
    vSig: float
    vLambda: float
    rho: float

    @classmethod
    def from_underlying(cls, u) -> 'HestonParams':
        return cls(S0=u.S0, r=u.r, q=u.q, v0=u.v0, alpha=u.alpha,
                   vTheta=u.vTheta, vSig=u.vSig, vLambda=u.vLambda, rho=u.rho)

    def to_underlying(self):
        from pyvol import cpy_credit as cc

        u = cc.Underlying()
        u.S0, u.r, u.q = self.S0, self.r, self.q
        u.v0, u.alpha, u.vTheta = self.v0, self.alpha, self.vTheta
        u.vSig, u.vLambda, u.rho = self.vSig, self.vLambda, self.rho
        return u


@dataclass(frozen=True)
class MarketSnapshot:
    '''Cleaned call-option chain for one symbol.'''
    symbol: str
    spot: float
    strikes: np.ndarray
    maturities: np.ndarray
    prices: np.ndarray
    rates: np.ndarray
    volumes: np.ndarray
    ivs: np.ndarray          # Black-Scholes IV of each mid quote (NaN when no-arb fails)
    synthetic: bool = False

    @property
    def n_quotes(self) -> int:
        return int(self.strikes.size)


@dataclass(frozen=True)
class AssetView:
    '''Structural-model output for one debt/assets ratio.'''
    debt_ratio: float
    horizon: float
    p_default: float
    asset: HestonParams


@dataclass(frozen=True)
class Surface:
    '''A gridded IV surface ready to plot. `ivs` is (len(maturities), len(strikes)).'''
    kind: str                # 'market' | 'equity' | 'asset'
    title: str
    xlabel: str
    xkey: str                # short label for tooltips, e.g. 'K' or 'K/A'
    strikes: np.ndarray
    maturities: np.ndarray
    ivs: np.ndarray
    points: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None  # raw (K, T, iv)


# --------------------------------------------------------------------------
# Vectorised Black-Scholes pricing / implied-vol inversion
# --------------------------------------------------------------------------

def bs_call_price(S, K, r, q, T, sigma):
    sigma = np.maximum(sigma, 1e-12)
    srt = sigma * np.sqrt(T)
    d1 = (np.log(S / K) + (r - q) * T) / srt + 0.5 * srt
    d2 = d1 - srt
    return S * np.exp(-q * T) * ndtr(d1) - K * np.exp(-r * T) * ndtr(d2)


def bs_implied_vol(price, S, K, r, q, T, lo=1e-4, hi=6.0, iters=48):
    '''
    Invert Black-Scholes for calls by bisection, vectorised over whole arrays.
    Quotes violating no-arbitrage bounds come back NaN, as do quotes with
    negligible vega (deep ITM/OTM), where the implied vol is not identifiable
    and inversion just amplifies price noise into spikes.
    '''
    price, K, r, T = np.broadcast_arrays(
        *(np.asarray(a, dtype=float) for a in (price, K, r, T)))
    with np.errstate(all='ignore'):
        intrinsic = np.maximum(S * np.exp(-q * T) - K * np.exp(-r * T), 0.0)
        upper = S * np.exp(-q * np.maximum(T, 0.0))
        valid = (price > intrinsic + 1e-10) & (price < upper - 1e-10) & (T > 0)

        lo_a = np.full(price.shape, lo)
        hi_a = np.full(price.shape, hi)
        for _ in range(iters):
            mid = 0.5 * (lo_a + hi_a)
            over = bs_call_price(S, K, r, q, T, mid) > price
            hi_a = np.where(over, mid, hi_a)
            lo_a = np.where(over, lo_a, mid)
        iv = 0.5 * (lo_a + hi_a)

        valid &= iv < hi * 0.99  # pinned at the bracket top => quote is junk
        # Identifiability: deep ITM a call is almost pure intrinsic value, so
        # price noise (or Heston quadrature error) inverts to huge spurious
        # IVs. Require a minimum of genuine time value relative to spot.
        valid &= (price - intrinsic) >= 1e-3 * upper
    return np.where(valid, iv, np.nan)


# --------------------------------------------------------------------------
# Worker-process jobs (top level so they pickle by reference under spawn)
# --------------------------------------------------------------------------

def worker_warmup(offline: bool) -> bool:
    '''Pre-import the heavy modules so the first real job is snappy.'''
    from pyvol import cpy_credit  # noqa: F401

    if not offline:
        from pyvol import market_data  # noqa: F401
    return True


def fetch_market(symbol: str) -> MarketSnapshot:
    '''Live option chain via yfinance + US Treasury curve (see heston_fit.py).'''
    from pyvol.market_data import get_call_information

    try:
        S0, df = get_call_information(symbol)
    except Exception as exc:
        raise RuntimeError(f'Could not fetch market data for {symbol!r}: {exc}') from exc

    df = df.dropna(subset=['strike', 'price', 'maturity', 'rf'])
    df = df[(df['price'] > 0) & (df['maturity'] > 1.5 / 365)]
    df = df[(df['strike'] > 0.3 * S0) & (df['strike'] < 3.0 * S0)]

    # Zero-volume quotes carry zero weight in the fit but still cost pricing
    # time, and their stale mids pollute the surface - drop them when we can.
    liquid = df[df['volume'].fillna(0) > 0]
    if len(liquid) >= max(MIN_QUOTES, 30):
        df = liquid
    if len(df) < MIN_QUOTES:
        raise RuntimeError(f'Only {len(df)} usable call quotes for {symbol!r}')

    K, T, P, RF = (np.ascontiguousarray(df[c].to_numpy(float))
                   for c in ('strike', 'maturity', 'price', 'rf'))
    V = df['volume'].fillna(0).to_numpy(float)
    if V.sum() <= 0:
        V = np.ones_like(V)

    ivs = bs_implied_vol(P, float(S0), K, RF, 0.0, T)
    return MarketSnapshot(symbol=symbol.upper(), spot=float(S0), strikes=K,
                          maturities=T, prices=P, rates=RF, volumes=V, ivs=ivs)


def synthetic_market(symbol: str) -> MarketSnapshot:
    '''Offline stand-in: a chain priced under a preset Heston model + 1% noise.'''
    from pyvol import cpy_credit as cc

    u = cc.Underlying()
    u.S0, u.q = 60.0, 0.0
    u.v0, u.alpha, u.vTheta = 0.16, 1.8, 0.12
    u.vSig, u.vLambda, u.rho = 0.45, 0.0, -0.55

    Ts = np.array([2 / 52, 1 / 12, 2 / 12, 4 / 12, 7 / 12, 10 / 12, 1.25, 1.75, 2.25])
    Ks = np.linspace(33.0, 96.0, 22)
    rng = np.random.default_rng(20260713)

    strikes, mats, prices, rates, vols = [], [], [], [], []
    for T in Ts:
        rf = 0.036 + 0.009 * (1 - np.exp(-1.2 * T))
        u.r = rf
        for K in Ks:
            p = cc.get_Heston_call_price(u, float(K), float(T))
            # Perturb the *time value*, not the price: quote noise lives in
            # the vol, and intrinsic value is never up for negotiation.
            intrinsic = max(u.S0 - K * np.exp(-rf * T), 0.0)
            tv = max(p - intrinsic, 0.0) * (1 + rng.normal(0, 0.03))
            strikes.append(K)
            mats.append(T)
            rates.append(rf)
            prices.append(max(intrinsic + tv, 0.005))
            vols.append(5 + round(400 * np.exp(-((K - u.S0) / (0.35 * u.S0)) ** 2)))

    K, T, P, RF, V = (np.asarray(a, dtype=float)
                      for a in (strikes, mats, prices, rates, vols))
    ivs = bs_implied_vol(P, u.S0, K, RF, 0.0, T)
    return MarketSnapshot(symbol=symbol.upper(), spot=u.S0, strikes=K, maturities=T,
                          prices=P, rates=RF, volumes=V, ivs=ivs, synthetic=True)


def fit_heston(snap: MarketSnapshot) -> HestonParams:
    from pyvol import cpy_credit as cc

    u = cc.fit_Heston(snap.spot, snap.strikes, snap.rates, snap.maturities,
                      snap.prices, snap.volumes)
    return HestonParams.from_underlying(u)


def derive_asset(equity: HestonParams, debt_ratio: float, horizon: float = 1.0) -> AssetView:
    from pyvol import cpy_credit as cc

    p, av = cc.get_Heston_default_probability(
        equity.to_underlying(), 1.0, float(debt_ratio), float(horizon))
    return AssetView(debt_ratio=float(debt_ratio), horizon=float(horizon),
                     p_default=float(p), asset=HestonParams.from_underlying(av))


def _valid_quote_mask(snap: MarketSnapshot) -> np.ndarray:
    return np.isfinite(snap.ivs) & (snap.ivs > 0.01) & (snap.ivs < 4.0)


def _rate_curve(snap: MarketSnapshot, grid_T: np.ndarray) -> np.ndarray:
    '''Interpolate the snapshot's per-quote risk-free rates onto grid maturities.'''
    uT, inv = np.unique(np.round(snap.maturities, 6), return_inverse=True)
    uR = np.bincount(inv, weights=snap.rates) / np.bincount(inv)
    return np.interp(grid_T, uT, uR)


def _model_grid(snap: MarketSnapshot):
    m = _valid_quote_mask(snap)
    K, T = snap.strikes[m], snap.maturities[m]
    nK, nT = MODEL_GRID
    gK = np.linspace(K.min(), K.max(), nK)
    gT = np.linspace(max(T.min(), 0.02), T.max(), nT)
    return gK, gT


def market_surface(snap: MarketSnapshot) -> Surface:
    '''Interpolate scattered quote IVs onto a regular grid (NaN outside the hull).'''
    m = _valid_quote_mask(snap)
    if m.sum() < 12:
        raise RuntimeError('Too few arbitrage-consistent quotes to build a market surface')

    K, T, iv = snap.strikes[m], snap.maturities[m], snap.ivs[m]
    nK, nT = MARKET_GRID
    gK = np.linspace(K.min(), K.max(), nK)
    gT = np.linspace(T.min(), T.max(), nT)
    KK, TT = np.meshgrid(gK, gT)
    Z = griddata((K, T), iv, (KK, TT), method='linear', rescale=True)

    return Surface(kind='market', title='Market implied volatility (mid quotes)',
                   xlabel='Strike ($)', xkey='K', strikes=gK, maturities=gT,
                   ivs=Z, points=(K, T, iv))


def heston_iv_grid(params: HestonParams, strikes: np.ndarray, maturities: np.ndarray,
                   rates: np.ndarray | None = None) -> np.ndarray:
    '''Price a strike x maturity grid under Heston (C++), invert to BS IVs.'''
    from pyvol import cpy_credit as cc

    if rates is None:
        rates = np.full(maturities.shape, params.r)

    u = params.to_underlying()
    prices = np.empty((maturities.size, strikes.size))
    for i, (T, rf) in enumerate(zip(maturities, rates)):
        u.r = rf  # get_Heston_call_price copies the underlying, mutation is safe
        for j, K in enumerate(strikes):
            prices[i, j] = cc.get_Heston_call_price(u, float(K), float(T))

    iv = bs_implied_vol(prices, params.S0, strikes[None, :],
                        rates[:, None], params.q, maturities[:, None])
    return np.where((iv > 0.005) & (iv < 4.0), iv, np.nan)


def equity_surface_job(snap: MarketSnapshot, equity: HestonParams) -> Surface:
    gK, gT = _model_grid(snap)
    Z = heston_iv_grid(equity, gK, gT, _rate_curve(snap, gT))
    return Surface(kind='equity', title='Heston fit - equity implied volatility',
                   xlabel='Strike ($)', xkey='K', strikes=gK, maturities=gT, ivs=Z)


def asset_surface_job(snap: MarketSnapshot, av: AssetView) -> Surface:
    _, gT = _model_grid(snap)
    nK = MODEL_GRID[0]
    gK = np.linspace(*ASSET_MONEYNESS, nK)  # strikes on the unit-normalised asset
    Z = heston_iv_grid(av.asset, gK, gT)
    return Surface(kind='asset',
                   title=f'Heston-implied asset volatility · debt/assets = {av.debt_ratio:.2f}',
                   xlabel='Strike / asset value', xkey='K/A',
                   strikes=gK, maturities=gT, ivs=Z)


# --------------------------------------------------------------------------
# The model facade owned by the GUI process
# --------------------------------------------------------------------------

class VolSurfaceModel:
    '''
    Submits compute jobs to a single spawn worker (jobs serialise in order,
    which the controller's pipeline relies on) and returns futures.
    '''

    def __init__(self, offline: bool = False):
        self.offline = offline
        self._pool = ProcessPoolExecutor(max_workers=1,
                                         mp_context=mp.get_context('spawn'))

    def warm_up(self):
        return self._pool.submit(worker_warmup, self.offline)

    def fetch(self, symbol: str):
        job = synthetic_market if self.offline else fetch_market
        return self._pool.submit(job, symbol)

    def fit(self, snap: MarketSnapshot):
        return self._pool.submit(fit_heston, snap)

    def derive_asset(self, equity: HestonParams, debt_ratio: float, horizon: float = 1.0):
        return self._pool.submit(derive_asset, equity, debt_ratio, horizon)

    def equity_iv_surface(self, snap: MarketSnapshot, equity: HestonParams):
        return self._pool.submit(equity_surface_job, snap, equity)

    def asset_iv_surface(self, snap: MarketSnapshot, av: AssetView):
        return self._pool.submit(asset_surface_job, snap, av)

    def market_surface(self, snap: MarketSnapshot) -> Surface:
        return market_surface(snap)  # pure numpy, fast enough for the GUI thread

    def shutdown(self):
        self._pool.shutdown(wait=False, cancel_futures=True)
