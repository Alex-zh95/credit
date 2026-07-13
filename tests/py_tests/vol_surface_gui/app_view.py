'''
View layer: PySide6 window with a matplotlib QtAgg canvas.

Qt is used because this project's pdm-managed interpreter is a
python-build-standalone build whose statically linked `_tkinter` cannot load
matplotlib's TkAgg C bridge. The QtAgg backend blits the Agg buffer straight
into the widget via QImage (no per-frame image encode/transfer), so
interaction runs at the cost of the Agg draw itself. mplot3d's native mouse
handling provides rotation (left-drag) and distance zoom (right-drag); the
scroll wheel zooms by scaling the axes limits. While a gesture is active a
coarse, stroke-free twin of the surface is swapped in to halve the draw time.

Holds no business logic - the controller injects data through `show_surface`
/ `set_*` and receives user intent through the `on_*` callbacks. `after()`
mirrors tkinter's scheduling API so the controller stays toolkit-agnostic.
'''

from __future__ import annotations

import os
import sys

os.environ.setdefault('QT_API', 'pyside6')


def _defuse_shiboken_hook():
    '''
    PySide6's shiboken "feature" import hook introspects every module imported
    after PySide6 and crashes on six.moves' fake lazy modules (PYSIDE-2226
    family). If PySide6 is already loaded, wrap the hook so it swallows
    introspection errors instead of aborting unrelated imports.
    '''
    feature = sys.modules.get('shibokensupport.feature')
    if feature is None or getattr(feature, '_viz_defused', False):
        return

    original = feature.feature_imported

    def safe_feature_imported(module, _original=original):
        try:
            _original(module)
        except Exception:
            pass

    feature.feature_imported = safe_feature_imported
    feature._viz_defused = True


_defuse_shiboken_hook()

import numpy as np
import matplotlib

import app_theme as theme

matplotlib.use('QtAgg')
matplotlib.rcParams.update({
    'font.size': 9.5,
    'font.sans-serif': ['Helvetica Neue', 'Arial', 'DejaVu Sans'],
    'grid.color': theme.GRID,
    'grid.linewidth': 0.7,
    'axes.labelcolor': theme.MUTED,
    'xtick.color': theme.MUTED,
    'ytick.color': theme.MUTED,
    'axes.edgecolor': theme.BASELINE,
    'text.color': theme.INK,
    'figure.facecolor': theme.SURFACE,
    'savefig.facecolor': theme.SURFACE,
})

# The Qt backend pulls in six.moves via matplotlib.dates -> dateutil. Load
# that chain before PySide6 installs its import hook (see _defuse_shiboken_hook
# for the case where PySide6 got imported first).
import matplotlib.dates  # noqa: F401  isort: skip

from matplotlib import colors as mcolors
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.mplot3d import proj3d

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QFontDatabase
from PySide6.QtWidgets import (
    QButtonGroup, QComboBox, QDoubleSpinBox, QFileDialog, QFrame, QHBoxLayout,
    QLabel, QMainWindow, QPushButton, QRadioButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget)

from app_model import AssetView, HestonParams, Surface

_PARAM_ROWS = (
    ('spot', 'Spot price'),
    ('v0', 'Spot variance v0'),
    ('alpha', 'Mean-rev rate α'),
    ('vTheta', 'Long-run var θ'),
    ('vSig', 'Vol of vol σ'),
    ('vLambda', 'Price of vol λ'),
    ('rho', 'Corr ρ'),
    ('pdef', 'P(default)'),
)

DEFAULT_ELEV, DEFAULT_AZIM = 22.0, -60.0


class HoverProbe:
    '''
    Pure hit-testing: grid nodes projected to display space, nearest-node
    lookup within a 24px target. The projection is cached against the view
    matrix + canvas size, so a mouse move costs one vectorised distance pass.
    '''

    RADIUS_PX = 24
    DEPTH_TIE_PX = 8

    def __init__(self):
        self.ax = None
        self._pts = None
        self._proj_key = None
        self._proj = None

    def attach(self, ax, X, Y, Z):
        self.ax = ax
        z = np.ma.filled(np.ma.masked_invalid(Z), np.nan).ravel()
        ok = np.isfinite(z)
        self._pts = np.column_stack([np.ravel(X)[ok], np.ravel(Y)[ok], z[ok]])
        self._proj_key = None

    def _projected(self, size):
        M = self.ax.get_proj()
        key = (M.tobytes(), size)
        if key != self._proj_key:
            x2, y2, z2 = proj3d.proj_transform(
                self._pts[:, 0], self._pts[:, 1], self._pts[:, 2], M)
            disp = self.ax.transData.transform(np.column_stack([x2, y2]))
            self._proj = (disp, np.asarray(z2))
            self._proj_key = key
        return self._proj

    def pick(self, x_disp, y_disp, size):
        '''Return (data_point, display_xy) for the node under the cursor, or None.'''
        if self.ax is None or self._pts is None or not len(self._pts):
            return None
        disp, depth = self._projected(size)
        d2 = (disp[:, 0] - x_disp) ** 2 + (disp[:, 1] - y_disp) ** 2
        nearest = int(np.argmin(d2))
        if d2[nearest] > self.RADIUS_PX ** 2:
            return None
        # Among near-tied candidates prefer the one facing the camera
        # (mplot3d draws large projected z first, i.e. farthest away).
        near = np.flatnonzero(d2 <= (np.sqrt(d2[nearest]) + self.DEPTH_TIE_PX) ** 2)
        idx = int(near[np.argmin(depth[near])])
        return self._pts[idx], disp[idx]


class HoverLayer:
    '''
    Mouse-over tooltip rendered with blitting: the scene is cached once per
    draw and only the marker + annotation are re-composited per mouse move,
    so hovering never triggers a full matplotlib redraw.
    '''

    def __init__(self, canvas, probe: HoverProbe, is_busy, on_point=None):
        self.canvas = canvas
        self.probe = probe
        self.is_busy = is_busy          # suppress hover during drags/toolbar modes
        self.on_point = on_point        # mirrors the tooltip into the status bar
        self.ax = None
        self._fmt = None
        self._bg = None
        self._ann = None
        self._dot = None
        canvas.mpl_connect('draw_event', self._on_draw)
        canvas.mpl_connect('motion_notify_event', self._on_motion)
        canvas.mpl_connect('figure_leave_event', lambda _e: self.hide())
        canvas.mpl_connect('axes_leave_event', lambda _e: self.hide())

    def attach(self, ax, fmt):
        self.ax = ax
        self._fmt = fmt
        self._bg = None
        self._dot, = ax.plot([0.0], [0.0], [0.0], marker='o', ms=8, ls='',
                             mfc=theme.ACCENT, mec=theme.SURFACE, mew=1.6,
                             animated=True, zorder=20, visible=False)
        self._ann = ax.annotate(
            '', xy=(0, 0), xycoords='figure pixels',
            xytext=(14, 14), textcoords='offset points',
            fontsize=9, color=theme.INK, linespacing=1.45,
            bbox=dict(boxstyle='round,pad=0.45', fc=theme.TOOLTIP_BG,
                      ec=theme.BASELINE, lw=0.8, alpha=0.97),
            annotation_clip=False, animated=True, zorder=21, visible=False)

    def _on_draw(self, _event):
        self._bg = self.canvas.copy_from_bbox(self.canvas.figure.bbox)

    def _on_motion(self, event):
        if self.ax is None:
            return
        if event.inaxes is not self.ax or event.button is not None or self.is_busy():
            self.hide()
            return
        size = tuple(self.canvas.figure.bbox.size)
        hit = self.probe.pick(event.x, event.y, size)
        if hit is None:
            self.hide()
            return
        (x, y, z), (px, py) = hit
        self._dot.set_data_3d([x], [y], [z])
        self._ann.xy = (px, py)
        flip = px > self.canvas.figure.bbox.width - 170
        self._ann.set_horizontalalignment('right' if flip else 'left')
        self._ann.xyann = (-14, 14) if flip else (14, 14)
        self._ann.set_text(self._fmt(x, y, z))
        self._dot.set_visible(True)
        self._ann.set_visible(True)
        self._blit()
        if self.on_point:
            self.on_point((x, y, z))

    def hide(self):
        if self._ann is not None and (self._ann.get_visible() or self._dot.get_visible()):
            self._ann.set_visible(False)
            self._dot.set_visible(False)
            self._blit()
        if self.on_point:
            self.on_point(None)

    def _blit(self):
        if self._bg is None:
            self.canvas.draw()  # triggers _on_draw, which caches the background
            return
        self.canvas.restore_region(self._bg)
        if self._dot.get_visible():
            self.ax.draw_artist(self._dot)
            self.ax.draw_artist(self._ann)
        self.canvas.blit(self.canvas.figure.bbox)


class VolSurfaceView(QMainWindow):
    KINDS = (('market', 'Market IV'),
             ('equity', 'Heston fit IV'),
             ('asset', 'Asset IV (Heston)'))
    SYMBOL_PRESETS = ('CCJ', 'AAPL', 'MSFT', 'NVDA', 'TSLA', 'SPY')
    DEBT_DEFAULT = 0.75

    def __init__(self, offline: bool = False):
        super().__init__()
        self.setWindowTitle('Volatility surface explorer — pyvol / Heston demo')
        self.resize(1280, 860)

        # Controller-injected callbacks
        self.on_load = None
        self.on_kind = None
        self.on_debt = None
        self.on_close = None

        self._offline = offline
        self._surf_artist = None
        self._lod_artist = None            # coarse twin shown while interacting
        self._scatter_artist = None
        self._ax = None
        self._home = None                  # (xlim, ylim, zlim) at plot time
        self._view_angles = (DEFAULT_ELEV, DEFAULT_AZIM)
        self._interacting = False
        self._mouse_down = False
        self._readout_fmt = None
        self._tooltip_fmt = None

        self._refine_timer = QTimer(self)  # swap back to full quality after a lull
        self._refine_timer.setSingleShot(True)
        self._refine_timer.timeout.connect(self._end_interaction)
        self._debt_timer = QTimer(self)    # debounce: each recompute costs ~1 s
        self._debt_timer.setSingleShot(True)
        self._debt_timer.timeout.connect(self._fire_debt)

        self.probe = HoverProbe()
        self._build_ui()
        self.hover = HoverLayer(self.canvas, self.probe,
                                is_busy=self._pointer_busy,
                                on_point=self._mirror_readout)

        self.canvas.mpl_connect('scroll_event', self._on_scroll)
        self.canvas.mpl_connect('button_press_event', self._on_press)
        self.canvas.mpl_connect('button_release_event', self._on_release)

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        central = QWidget(self)
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(10, 8, 10, 4)

        bar = QHBoxLayout()
        root.addLayout(bar)

        bar.addWidget(QLabel('Ticker'))
        self._symbol_combo = QComboBox()
        self._symbol_combo.setEditable(True)
        self._symbol_combo.addItems(self.SYMBOL_PRESETS)
        self._symbol_combo.lineEdit().returnPressed.connect(self._fire_load)
        bar.addWidget(self._symbol_combo)

        load_btn = QPushButton('Load && fit')
        load_btn.clicked.connect(self._fire_load)
        bar.addWidget(load_btn)

        bar.addWidget(self._separator())
        bar.addWidget(QLabel('Surface'))
        self._kind_buttons = {}
        self._kind_group = QButtonGroup(self)
        for kind, label in self.KINDS:
            rb = QRadioButton(label)
            rb.setEnabled(False)
            # `clicked` fires on user action only, so programmatic setChecked
            # from the controller never loops back into on_kind
            rb.clicked.connect(lambda _=False, k=kind: self.on_kind and self.on_kind(k))
            self._kind_group.addButton(rb)
            bar.addWidget(rb)
            self._kind_buttons[kind] = rb
        self._kind_buttons['market'].setChecked(True)

        bar.addWidget(self._separator())
        bar.addWidget(QLabel('Debt / assets'))
        self._debt_spin = QDoubleSpinBox()
        self._debt_spin.setRange(0.05, 0.95)
        self._debt_spin.setSingleStep(0.05)
        self._debt_spin.setDecimals(2)
        self._debt_spin.setValue(self.DEBT_DEFAULT)
        self._debt_spin.setKeyboardTracking(False)
        self._debt_spin.valueChanged.connect(
            lambda _v: self._debt_timer.start(350))
        bar.addWidget(self._debt_spin)

        bar.addWidget(self._separator())
        reset_btn = QPushButton('Reset view')
        reset_btn.clicked.connect(self._reset_view)
        bar.addWidget(reset_btn)
        save_btn = QPushButton('Save PNG…')
        save_btn.clicked.connect(self._save_png)
        bar.addWidget(save_btn)

        bar.addStretch(1)
        if self._offline:
            offline_label = QLabel('OFFLINE — synthetic data')
            offline_label.setStyleSheet(f'color: {theme.INK_SECONDARY};')
            bar.addWidget(offline_label)

        body = QHBoxLayout()
        root.addLayout(body, stretch=1)

        plot_col = QVBoxLayout()
        body.addLayout(plot_col, stretch=1)
        self.figure = Figure(figsize=(9.2, 6.6), dpi=100, facecolor=theme.SURFACE)
        self.canvas = FigureCanvasQTAgg(self.figure)
        plot_col.addWidget(self.canvas, stretch=1)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        plot_col.addWidget(self.toolbar)

        side = QVBoxLayout()
        body.addLayout(side)
        heading = QLabel('Heston parameters')
        heading.setStyleSheet(f'color: {theme.INK_SECONDARY};')
        side.addWidget(heading)

        self._table = QTableWidget(len(_PARAM_ROWS), 3)
        self._table.setHorizontalHeaderLabels(['', 'Equity', 'Asset'])
        self._table.verticalHeader().setVisible(False)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self._table.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        for row, (_iid, label) in enumerate(_PARAM_ROWS):
            self._table.setItem(row, 0, QTableWidgetItem(label))
            for col in (1, 2):
                item = QTableWidgetItem('—')
                item.setTextAlignment(
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self._table.setItem(row, col, item)
        self._table.setColumnWidth(0, 138)
        self._table.setColumnWidth(1, 96)
        self._table.setColumnWidth(2, 96)
        self._table.setFixedWidth(348)
        self._table.setFixedHeight(
            self._table.horizontalHeader().height()
            + sum(self._table.rowHeight(r) for r in range(len(_PARAM_ROWS))) + 4)
        side.addWidget(self._table)

        hint = QLabel('Drag to rotate, scroll to zoom,\nhover a node for its value.\n'
                      'Asset columns use assets\nnormalised to 1.0.')
        hint.setStyleSheet(f'color: {theme.MUTED};')
        side.addWidget(hint)
        side.addStretch(1)

        self._status = QLabel('Starting compute worker…')
        self._status.setStyleSheet(f'color: {theme.INK_SECONDARY};')
        self.statusBar().addWidget(self._status)
        self._readout = QLabel('')
        self._readout.setFont(
            QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.statusBar().addPermanentWidget(self._readout)

    @staticmethod
    def _separator():
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.VLine)
        sep.setFrameShadow(QFrame.Shadow.Sunken)
        return sep

    # --------------------------------------------- tkinter-compatible shims
    def after(self, ms: int, callback):
        '''Scheduling shim so the controller stays toolkit-agnostic.'''
        QTimer.singleShot(ms, callback)

    def destroy(self):
        self.close()

    def mainloop(self):  # pragma: no cover - the app runs QApplication.exec()
        raise RuntimeError('run the Qt event loop via QApplication.exec()')

    def closeEvent(self, event):
        if self.on_close:
            self.on_close()
        event.accept()

    # ------------------------------------------------------- user intent
    @property
    def symbol(self) -> str:
        return self._symbol_combo.currentText().strip().upper()

    @property
    def kind(self) -> str:
        for kind, button in self._kind_buttons.items():
            if button.isChecked():
                return kind
        return 'market'

    @kind.setter
    def kind(self, value: str):
        self._kind_buttons[value].setChecked(True)

    @property
    def debt(self) -> float | None:
        return round(self._debt_spin.value(), 4)

    def _fire_load(self):
        if self.on_load:
            self.on_load(self.symbol)

    def _fire_debt(self):
        if self.on_debt and self.debt is not None:
            self.on_debt(self.debt)

    # ----------------------------------------------------- 3D interaction
    def _pointer_busy(self) -> bool:
        return bool(self._mouse_down or self._interacting or self.toolbar.mode)

    def _begin_interaction(self):
        '''Swap in the coarse surface so interaction frames render cheaply.'''
        self._refine_timer.stop()
        if self._interacting or self._lod_artist is None:
            return
        self._interacting = True
        self._surf_artist.set_visible(False)
        self._lod_artist.set_visible(True)
        if self._scatter_artist is not None:
            self._scatter_artist.set_visible(False)

    def _end_interaction(self):
        '''Swap back to full quality and redraw.'''
        self._refine_timer.stop()
        if not self._interacting:
            return
        self._interacting = False
        self._surf_artist.set_visible(True)
        self._lod_artist.set_visible(False)
        if self._scatter_artist is not None:
            self._scatter_artist.set_visible(True)
        self.canvas.draw_idle()

    def _on_press(self, event):
        # mplot3d handles the rotation itself; we only manage LOD + tooltip
        if event.inaxes is not self._ax or self.toolbar.mode:
            return
        self._mouse_down = True
        self.hover.hide()
        self._begin_interaction()

    def _on_release(self, _event):
        if not self._mouse_down:
            return
        self._mouse_down = False
        if self._ax is not None:
            self._view_angles = (self._ax.elev, self._ax.azim)
        self._end_interaction()

    def _on_scroll(self, event):
        if self._ax is None or event.inaxes is not self._ax:
            return
        steps = max(-5.0, min(5.0, event.step))
        if not steps:
            return
        self._begin_interaction()
        self._zoom(0.94 ** steps)
        # Trackpads keep emitting momentum events; refine once they stop
        self._refine_timer.start(180)

    def _zoom(self, factor):
        for get_lim, set_lim in ((self._ax.get_xlim3d, self._ax.set_xlim3d),
                                 (self._ax.get_ylim3d, self._ax.set_ylim3d),
                                 (self._ax.get_zlim3d, self._ax.set_zlim3d)):
            lo, hi = get_lim()
            centre, half = (lo + hi) / 2, (hi - lo) / 2 * factor
            set_lim(centre - half, centre + half)
        self.canvas.draw_idle()

    def _reset_view(self):
        if self._ax is None:
            return
        self._end_interaction()
        self._ax.view_init(elev=DEFAULT_ELEV, azim=DEFAULT_AZIM)
        self._view_angles = (DEFAULT_ELEV, DEFAULT_AZIM)
        if self._home is not None:
            self._ax.set_xlim3d(self._home[0])
            self._ax.set_ylim3d(self._home[1])
            self._ax.set_zlim3d(self._home[2])
        self.canvas.draw_idle()

    def _save_png(self):
        if self._ax is None:
            return
        path, _filter = QFileDialog.getSaveFileName(
            self, 'Save PNG', 'vol_surface.png', 'PNG image (*.png)')
        if path:
            self.figure.savefig(path, dpi=160)
            self.set_ready(f'Saved {path}')

    def _mirror_readout(self, point):
        if point is None or self._readout_fmt is None:
            self._readout.setText('')
        else:
            self._readout.setText(self._readout_fmt(*point))

    # --------------------------------------------------------- rendering
    def show_surface(self, surface: Surface, symbol: str):
        if self._ax is not None:
            self._view_angles = (self._ax.elev, self._ax.azim)
        self._refine_timer.stop()
        self._interacting = False

        self.figure.clf()
        ax = self.figure.add_subplot(projection='3d', facecolor=theme.SURFACE)
        self._ax = ax

        KK, TT = np.meshgrid(surface.strikes, surface.maturities)
        Z = np.ma.masked_invalid(surface.ivs)
        vmin, vmax = np.nanpercentile(surface.ivs, [2, 98])
        norm = mcolors.Normalize(vmin=vmin, vmax=max(vmax, vmin + 1e-6))
        cmap = theme.sequential_cmap()
        self._surf_artist = ax.plot_surface(
            KK, TT, Z, cmap=cmap, norm=norm, rcount=64, ccount=64,
            linewidth=0.3, edgecolor=mcolors.to_rgba(theme.SURFACE, 0.35),
            antialiased=True)
        # Coarse, stroke-free twin swapped in during rotation/zoom: it halves
        # the Agg draw time, which is what bounds the interactive frame rate
        self._lod_artist = ax.plot_surface(
            KK, TT, Z, cmap=cmap, norm=norm,
            rcount=28, ccount=28, linewidth=0, antialiased=False, visible=False)

        self._scatter_artist = None
        if surface.points is not None:
            pK, pT, pIV = surface.points
            self._scatter_artist = ax.scatter(
                pK, pT, pIV, s=7, c=theme.INK_SECONDARY, alpha=0.30,
                linewidths=0, depthshade=False)

        ax.set_xlabel(surface.xlabel, labelpad=8)
        ax.set_ylabel('Maturity (years)', labelpad=8)
        ax.set_zlabel('Implied vol', labelpad=6)
        ax.zaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=0))
        ax.set_title(f'{symbol} — {surface.title}', loc='left',
                     color=theme.INK_SECONDARY, fontsize=11, pad=14)
        for axis_name in ('x', 'y', 'z'):
            ax.tick_params(axis=axis_name, colors=theme.MUTED, labelsize=8.5)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_pane_color(mcolors.to_rgba(theme.PAGE, 0.85))
            axis.pane.set_edgecolor(theme.GRID)
        ax.set_box_aspect((1.35, 1.0, 0.72))
        ax.view_init(elev=self._view_angles[0], azim=self._view_angles[1])
        self.figure.subplots_adjust(left=0.0, right=0.98, bottom=0.02, top=0.92)
        self._home = (ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d())

        xkey = surface.xkey
        self.probe.attach(ax, KK, TT, Z)
        self.hover.attach(
            ax, fmt=lambda x, y, z: f'{z:.1%}\n{xkey} {x:,.2f} · {y:.2f}y')
        self._readout_fmt = lambda x, y, z: (
            f'IV {z:.1%}   {xkey} {x:,.2f}   T {y:.2f}y')

        self.toolbar.update()  # reset the home/back navigation stack
        self.canvas.draw_idle()

    # ------------------------------------------------- controller inputs
    def set_kind_enabled(self, kind: str, enabled: bool):
        self._kind_buttons[kind].setEnabled(enabled)

    def set_busy(self, message: str, dim: bool = False):
        self._status.setStyleSheet(f'color: {theme.INK_SECONDARY};')
        self._status.setText(message)
        if dim and self._surf_artist is not None:
            for artist in (self._surf_artist, self._lod_artist):
                artist.set_alpha(0.55)  # hold the stale frame, dimmed
            self.canvas.draw_idle()

    def set_ready(self, message: str = 'Ready'):
        self._status.setStyleSheet(f'color: {theme.INK_SECONDARY};')
        self._status.setText(message)
        if self._surf_artist is not None and self._surf_artist.get_alpha() is not None:
            for artist in (self._surf_artist, self._lod_artist):
                artist.set_alpha(None)
            self.canvas.draw_idle()

    def set_error(self, message: str):
        self._status.setStyleSheet(f'color: {theme.CRITICAL};')
        self._status.setText(message)

    def _set_cell(self, iid: str, column: str, text: str):
        row = next(i for i, (r, _l) in enumerate(_PARAM_ROWS) if r == iid)
        col = {'param': 0, 'equity': 1, 'asset': 2}[column]
        self._table.item(row, col).setText(text)

    def set_symbol_info(self, symbol: str, spot: float, synthetic: bool):
        # synthetic data is already flagged by the OFFLINE badge in the toolbar
        self._set_cell('spot', 'equity', f'${spot:,.2f}')
        self._set_cell('spot', 'asset', '1.00 (norm.)')

    def set_equity_params(self, p: HestonParams | None):
        for iid in ('v0', 'alpha', 'vTheta', 'vSig', 'vLambda', 'rho'):
            self._set_cell(iid, 'equity', '—' if p is None else f'{getattr(p, iid):,.5f}')

    def set_asset_params(self, av: AssetView | None):
        for iid in ('v0', 'alpha', 'vTheta', 'vSig', 'vLambda', 'rho'):
            self._set_cell(iid, 'asset',
                           '—' if av is None else f'{getattr(av.asset, iid):,.5f}')
        if av is None:
            self._set_cell('pdef', 'param', 'P(default)')
            self._set_cell('pdef', 'asset', '—')
        else:
            self._set_cell('pdef', 'param', f'P(default, {av.horizon:g}y)')
            self._set_cell('pdef', 'asset', f'{av.p_default:.2%}')

    def reset_params(self):
        for iid, label in _PARAM_ROWS:
            self._set_cell(iid, 'param', label)
            self._set_cell(iid, 'equity', '—')
            self._set_cell(iid, 'asset', '—')
