# Volatility-surface explorer (GUI demo)

Interactive 3D volatility surfaces built on the `heston_fit.py` concept: fetch a
live call-option chain, fit the Heston model in C++, and derive the firm's
**asset** volatility surface through the structural default model.

Pipeline: ticker → option chain (yfinance + US Treasury curve) → mid-quote
Black-Scholes IVs → `cc.fit_Heston` → `cc.get_Heston_default_probability` →
Heston-priced IV grids.

The UI is dark mode throughout, themed after
[Kanagawa](https://github.com/rebelot/kanagawa.nvim): sumiInk surfaces,
fuji/oldWhite ink, and a waveBlue→springBlue sequential ramp in which low
volatility recedes into the surface (contrast and ramp monotonicity checked
against the dark surfaces).

## Run

Requires the built `cpy_credit` extension (see repo `AGENTS.md`) and the pdm
virtualenv (which includes PySide6). From the repo root:

```sh
python tests/py_tests/vol_surface_gui                 # live data, loads CCJ
python tests/py_tests/vol_surface_gui --symbol AAPL   # pick the startup ticker
python tests/py_tests/vol_surface_gui --offline       # synthetic chain, no network
python tests/py_tests/vol_surface_gui --self-test     # headless-ish check + PNGs
```

## Controls

- **Ticker** — type any Yahoo Finance symbol (or pick a preset) and *Load & fit*.
- **Surface** — toggle between the market IV surface (gridded mid-quote IVs with
  the raw quotes scattered underneath), the fitted Heston equity IV surface, and
  the Heston-implied **asset** IV surface for the chosen debt/assets ratio.
- **Debt / assets** — leverage for the structural model; changing it re-derives
  the asset parameters, the 1-year default probability and the asset surface.
- **Plot** — left-drag to rotate (right-drag adjusts zoom distance), scroll
  wheel to zoom, toolbar for pan/zoom-to-rect/save, hover any grid node for a
  tooltip with the IV / strike / maturity at that point; the value is mirrored
  in the status bar. Fitted parameters sit in the right-hand panel.

## Architecture (MVC)

| File | Role |
|---|---|
| `app_model.py` | Data + compute: fetch, cleaning, vectorised Black-Scholes IV inversion, Heston fit/pricing jobs, worker-process pool. No GUI imports. |
| `app_view.py` | PySide6 window + matplotlib QtAgg canvas, hover layer, LOD swap during gestures, parameter panel. No business logic. |
| `app_controller.py` | Chains model futures, marshals results onto the GUI thread (via a toolkit-agnostic `after()` shim), caches per (surface kind, debt ratio), drops superseded generations. |
| `app_theme.py` | Kanagawa dark theme: chart chrome tokens, the sequential "wave" colormap, and the Qt palette/stylesheet applied at startup. |

## Why PySide6 / QtAgg

The pdm-managed interpreter is a python-build-standalone build whose
`_tkinter` is statically linked: the Tk C symbols are not exported, so
matplotlib's TkAgg bridge cannot load and a tkinter build must push frames
through a `tk.PhotoImage` (~18 ms/frame of pure transfer). The QtAgg backend
instead blits the Agg buffer via QImage at negligible cost, renders at the
retina device-pixel ratio, and provides mplot3d's native rotate plus the
standard toolbar. Measured on the same simulated drag: ~19 fps (Tk, 1×) →
**~36 fps at 2× resolution** (Qt).

One quirk: PySide6's shiboken import hook crashes on `six.moves`' fake lazy
modules, which the Qt backend reaches via `matplotlib.dates` → `dateutil`.
`app_view.py` pre-imports that chain (and defuses the hook if PySide6 was
imported first), so import order doesn't matter.

## Performance notes

- The pybind11 calls (`fit_Heston`, `get_Heston_call_price`) hold the GIL, so
  all C++ work runs in a **worker process** (spawned once, warmed at startup) —
  the UI never blocks; a background *thread* would still freeze the GUI. A
  live fit takes ~10–15 s on a few hundred quotes.
- Implied vols are inverted by **vectorised bisection** over the whole grid at
  once (`scipy.special.ndtr`), ~1 ms for a full surface. Nodes where the IV is
  not identifiable (negligible time value, e.g. deep ITM) become gaps rather
  than noise spikes.
- While rotating or zooming a **coarse, stroke-free twin** of the surface is
  swapped in (halving the Agg draw, which now bounds the frame rate) and full
  quality returns when the gesture ends — after a lull, for trackpad momentum.
- The mouse-over tooltip is **blitted**: the scene is cached per draw and only
  the marker/annotation are re-composited per mouse move; grid-node screen
  projections are cached against the view matrix.
- Surfaces and asset derivations are **memoised** per (kind, debt ratio);
  toggling views is instant, and a new Load cancels queued stale jobs.
