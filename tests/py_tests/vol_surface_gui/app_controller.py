'''
Controller layer: wires view events to model jobs and marshals worker
results back onto the Tk main thread.

Worker futures complete on an executor thread; their results are queued and
drained by a Tk `after` loop, so every view touch happens on the GUI thread.
A generation counter drops results that a newer "Load" has superseded, and
surfaces/asset derivations are cached per (kind, debt ratio) so toggling
between views is instant.
'''

from __future__ import annotations

import queue
import sys
import traceback
from dataclasses import dataclass, field

from app_model import AssetView, HestonParams, MarketSnapshot, Surface, VolSurfaceModel
from app_view import VolSurfaceView

READY_MSG = 'Ready — drag to rotate, scroll to zoom, hover for values'


@dataclass
class Session:
    symbol: str = ''
    snapshot: MarketSnapshot | None = None
    equity: HestonParams | None = None
    surfaces: dict = field(default_factory=dict)   # (kind, debt|None) -> Surface
    assets: dict = field(default_factory=dict)     # debt -> AssetView


class Controller:
    POLL_MS = 40

    def __init__(self, model: VolSurfaceModel, view: VolSurfaceView):
        self.model = model
        self.view = view
        self._queue: queue.Queue = queue.Queue()
        self._generation = 0
        self._pending: list = []
        self._inflight = 0
        self.session = Session()

        view.on_load = self.load_symbol
        view.on_kind = self.change_kind
        view.on_debt = self.change_debt
        view.on_close = self.shutdown

        self._chain(model.warm_up(), lambda _ok: None)
        view.after(self.POLL_MS, self._poll)

    # ----------------------------------------------------- future plumbing
    def _chain(self, future, callback):
        gen = self._generation
        self._pending.append(future)
        self._inflight += 1
        future.add_done_callback(lambda f: self._queue.put((gen, callback, f)))

    def _poll(self):
        processed, failed = 0, 0
        try:
            while True:
                gen, callback, future = self._queue.get_nowait()
                if gen != self._generation or future.cancelled():
                    continue
                self._inflight -= 1
                processed += 1
                try:
                    result = future.result()
                except Exception as exc:
                    traceback.print_exc(file=sys.stderr)
                    self.view.set_error(str(exc)[:160])
                    failed += 1
                    continue
                try:
                    callback(result)
                except Exception as exc:
                    traceback.print_exc(file=sys.stderr)
                    self.view.set_error(f'Internal error: {exc}'[:160])
                    failed += 1
            # not reached
        except queue.Empty:
            pass
        # Announce Ready only on the transition to idle, and never over an
        # error raised in this same drain - later ticks must not clobber
        # whatever the status bar currently says (errors, 'Saved ...', etc.)
        if (processed and not failed and self._inflight == 0
                and self.session.snapshot is not None):
            self.view.set_ready(READY_MSG)
        self.view.after(self.POLL_MS, self._poll)

    # ------------------------------------------------------------ pipeline
    def load_symbol(self, symbol: str):
        if not symbol:
            self.view.set_error('Enter a ticker symbol')
            return
        self._generation += 1
        for f in self._pending:
            f.cancel()
        self._pending.clear()
        self._inflight = 0
        self.session = Session(symbol=symbol)

        self.view.reset_params()
        for kind, _label in self.view.KINDS:
            self.view.set_kind_enabled(kind, False)
        self.view.kind = 'market'
        self.view.set_busy(f'Fetching option chain for {symbol}…', dim=True)
        try:
            self._chain(self.model.fetch(symbol), self._on_snapshot)
        except Exception as exc:  # e.g. BrokenProcessPool
            self.view.set_error(f'Compute worker unavailable: {exc}'[:160])

    def _on_snapshot(self, snap: MarketSnapshot):
        s = self.session
        s.snapshot = snap
        self.view.set_symbol_info(snap.symbol, snap.spot, snap.synthetic)

        surface = self.model.market_surface(snap)
        s.surfaces[('market', None)] = surface
        self.view.set_kind_enabled('market', True)
        if self.view.kind == 'market':
            self.view.show_surface(surface, snap.symbol)

        self.view.set_busy(
            f'Fitting Heston model to {snap.n_quotes} call quotes (C++)…')
        self._chain(self.model.fit(snap), self._on_fit)

    def _on_fit(self, equity: HestonParams):
        s = self.session
        s.equity = equity
        self.view.set_equity_params(equity)
        self.view.set_busy('Computing model IV surfaces…')
        self._chain(self.model.equity_iv_surface(s.snapshot, equity),
                    lambda surf: self._on_surface(('equity', None), surf))
        debt = self.view.debt
        if debt is not None:
            self._ensure_asset(debt)

    def _ensure_asset(self, debt: float):
        s = self.session
        if ('asset', debt) in s.surfaces:
            self._on_surface(('asset', debt), s.surfaces[('asset', debt)])
            return
        if debt in s.assets:
            self._on_asset(s.assets[debt])
            return
        self._chain(self.model.derive_asset(s.equity, debt), self._on_asset)

    def _on_asset(self, av: AssetView):
        s = self.session
        s.assets[av.debt_ratio] = av
        if av.debt_ratio == self.view.debt:
            self.view.set_asset_params(av)
        key = ('asset', av.debt_ratio)
        if key in s.surfaces:
            self._on_surface(key, s.surfaces[key])
        else:
            self._chain(self.model.asset_iv_surface(s.snapshot, av),
                        lambda surf, k=key: self._on_surface(k, surf))

    def _on_surface(self, key, surface: Surface):
        kind, debt = key
        s = self.session
        s.surfaces[key] = surface
        if debt is not None and debt != self.view.debt:
            return  # superseded debt ratio: cached for later, not shown
        self.view.set_kind_enabled(kind, True)
        if self.view.kind == kind:
            self.view.show_surface(surface, s.symbol.upper())

    # --------------------------------------------------------- view events
    def change_kind(self, kind: str):
        s = self.session
        key = (kind, self.view.debt if kind == 'asset' else None)
        surface = s.surfaces.get(key)
        if surface is not None:
            self.view.show_surface(surface, s.symbol.upper())
        elif kind == 'asset' and s.equity is not None and self.view.debt is not None:
            self.view.set_busy('Computing asset volatility surface…')
            self._ensure_asset(self.view.debt)

    def change_debt(self, debt: float):
        s = self.session
        if s.equity is None:
            return  # nothing fitted yet; the pipeline picks the ratio up later
        av = s.assets.get(debt)
        self.view.set_asset_params(av)
        if ('asset', debt) not in s.surfaces:
            self.view.set_kind_enabled('asset', False)
            self.view.set_busy(
                f'Deriving asset volatility for debt/assets = {debt:.2f}…')
        self._ensure_asset(debt)

    def shutdown(self):
        self._generation += 1  # drop anything still in flight
        self.model.shutdown()
