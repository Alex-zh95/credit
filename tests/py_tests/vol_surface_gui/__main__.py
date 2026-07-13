'''
Interactive volatility-surface explorer built on pyvol / cpy_credit.

Run from the repo root (the directory itself is the entry point):

    python tests/py_tests/vol_surface_gui                # live market data
    python tests/py_tests/vol_surface_gui --offline      # synthetic chain, no network
    python tests/py_tests/vol_surface_gui --self-test    # render snapshots and exit
'''

from __future__ import annotations

import argparse
import sys
import tempfile
import time
from pathlib import Path

from app_controller import Controller
from app_model import VolSurfaceModel
from app_view import VolSurfaceView

import app_theme as theme

from PySide6.QtWidgets import QApplication


class SelfTest:
    '''
    Drives the full offline pipeline, simulates hover/rotate/zoom on each
    surface and writes PNG snapshots, then closes the window. The process
    exits non-zero on timeout or a failed check.
    '''

    TIMEOUT_S = 240

    def __init__(self, view: VolSurfaceView, controller: Controller, out_dir: Path):
        self.view = view
        self.controller = controller
        self.out_dir = out_dir
        self.t0 = time.monotonic()
        self.ok = False
        view.after(500, self._wait_ready)

    def _wait_ready(self):
        session = self.controller.session
        done = {('market', None), ('equity', None),
                ('asset', self.view.debt)} <= set(session.surfaces)
        if done and self.controller._inflight == 0:
            self.view.after(100, self._run_checks)
        elif time.monotonic() - self.t0 > self.TIMEOUT_S:
            print('SELF-TEST FAILED: pipeline did not complete in time', file=sys.stderr)
            self.view.destroy()
        else:
            self.view.after(500, self._wait_ready)

    def _run_checks(self):
        try:
            for kind, _label in self.view.KINDS:
                self.view.kind = kind
                self.controller.change_kind(kind)
                self._snapshot(kind)
            self._check_rotate_zoom()
            window_png = self.out_dir / 'window.png'
            self.view.grab().save(str(window_png))
            print(f'snapshot: {window_png} (full window)')
            self.ok = True
            print(f'SELF-TEST OK ({time.monotonic() - self.t0:.1f}s)')
        except Exception as exc:
            print(f'SELF-TEST FAILED: {exc}', file=sys.stderr)
        finally:
            self.view.destroy()

    def _snapshot(self, kind: str):
        from matplotlib.backend_bases import MouseEvent

        view = self.view
        view.canvas.draw()

        # Simulate a mouse-over on a middle grid node to exercise the hover
        # hit-testing, then persist the tooltip into the snapshot.
        size = tuple(view.canvas.figure.bbox.size)
        disp, _depth = view.probe._projected(size)
        px, py = disp[len(disp) // 2]
        event = MouseEvent('motion_notify_event', view.canvas, px, py)
        view.canvas.callbacks.process('motion_notify_event', event)
        hover = view.hover
        if not hover._ann.get_visible():
            raise AssertionError(f'hover did not resolve a node on {kind!r}')
        tooltip = hover._ann.get_text().splitlines()[0]

        for artist in (hover._ann, hover._dot):
            artist.set_animated(False)
        view.canvas.draw()
        path = self.out_dir / f'{self.controller.session.symbol}_{kind}.png'
        # figure.dpi keeps the screen renderer's pixel space, so the
        # figure-pixel-anchored tooltip lands where it did on screen
        view.figure.savefig(path, dpi=view.figure.dpi)
        for artist in (hover._ann, hover._dot):
            artist.set_animated(True)
        hover.hide()
        print(f'snapshot: {path} (tooltip value {tooltip!r})')

    def _check_rotate_zoom(self):
        from matplotlib.backend_bases import MouseButton, MouseEvent

        view = self.view
        canvas = view.canvas
        ax = view._ax
        elev0, azim0 = ax.elev, ax.azim
        x0 = (ax.bbox.x0 + ax.bbox.x1) / 2
        y0 = (ax.bbox.y0 + ax.bbox.y1) / 2

        # mplot3d's native rotation, driven through synthetic mouse events
        for name, x, y in (('button_press_event', x0, y0),
                           ('motion_notify_event', x0 + 40, y0 + 15),
                           ('button_release_event', x0 + 40, y0 + 15)):
            event = MouseEvent(name, canvas, x, y, button=MouseButton.LEFT)
            canvas.callbacks.process(name, event)
            if name == 'button_press_event' and not view._interacting:
                raise AssertionError('press did not switch to the LOD surface')
        if (ax.elev, ax.azim) == (elev0, azim0):
            raise AssertionError('drag-rotation did not change the view angles')
        if view._interacting or not view._surf_artist.get_visible():
            raise AssertionError('full-quality surface not restored after drag')

        from types import SimpleNamespace

        xlim0 = ax.get_xlim3d()
        view._on_scroll(SimpleNamespace(inaxes=ax, step=1))
        if not (ax.get_xlim3d()[1] - ax.get_xlim3d()[0]) < (xlim0[1] - xlim0[0]):
            raise AssertionError('zoom did not shrink the axes limits')
        view._reset_view()
        print('rotate/zoom/reset interaction OK')


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--symbol', default='CCJ', help='ticker loaded on startup')
    parser.add_argument('--offline', action='store_true',
                        help='use a synthetic option chain instead of live data')
    parser.add_argument('--self-test', action='store_true',
                        help='run the offline pipeline, save snapshots, exit')
    parser.add_argument('--snapshot-dir', type=Path, default=None,
                        help='where --self-test writes its PNGs')
    args = parser.parse_args(argv)

    offline = args.offline or args.self_test
    app = QApplication(sys.argv[:1])
    theme.apply_qt_theme(app)
    model = VolSurfaceModel(offline=offline)
    view = VolSurfaceView(offline=offline)
    controller = Controller(model, view)
    view.show()

    self_test = None
    if args.self_test:
        out_dir = args.snapshot_dir or Path(tempfile.mkdtemp(prefix='vol_surface_'))
        out_dir.mkdir(parents=True, exist_ok=True)
        self_test = SelfTest(view, controller, out_dir)

    symbol = 'DEMO' if args.self_test else args.symbol
    view.after(200, lambda: controller.load_symbol(symbol))
    app.exec()
    model.shutdown()
    return 0 if (self_test is None or self_test.ok) else 1


if __name__ == '__main__':  # guard is required: the spawn worker re-imports us
    sys.exit(main())
